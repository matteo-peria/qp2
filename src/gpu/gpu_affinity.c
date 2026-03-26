
#include <sys/syscall.h>   /* SYS_getcpu — defined for every Linux arch */
#include <unistd.h>        /* syscall()                                  */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <limits.h>

/* ── PCI vendor IDs ────────────────────────────────────────────────── */
#define VENDOR_NVIDIA  0x10de
#define VENDOR_AMD     0x1002
#define VENDOR_INTEL   0x8086   /* Xe / Arc / Ponte Vecchio             */

/* ── PCI class codes (upper 2 bytes of the 3-byte class register) ─── */
#define CLASS_VGA      0x0300   /* VGA-compatible controller            */
#define CLASS_3D       0x0302   /* 3D controller (headless compute GPU) */
#define CLASS_DISPLAY  0x0380   /* Display controller                   */

#define MAX_GPUS  32
#define MAX_CPUS  4096

typedef struct {
    unsigned int domain, bus, slot, func;
    int          device_index;   /* 0-based, sorted by BDF               */
    int          numa_node;      /* -1 on non-NUMA systems               */
} gpu_info_t;

static gpu_info_t gpus[MAX_GPUS];
static int        ngpus = 0;
static int        cpu_to_gpu[MAX_CPUS];
static int        affinity_ready = 0;

/* ── tiny helpers ──────────────────────────────────────────────────── */

static int read_hex(const char *path, unsigned long *val)
{
    FILE *f = fopen(path, "r");
    if (!f) return -1;
    int ok = (fscanf(f, "%lx", val) == 1);
    fclose(f);
    return ok ? 0 : -1;
}

/* Returns the NUMA node that owns logical CPU core 'cpu',
   by parsing /sys/devices/system/node/nodeN/cpulist              */
static int cpu_numa_node(int cpu)
{
    char path[PATH_MAX], buf[8192];
    for (int n = 0; n < 256; n++) {
        snprintf(path, sizeof(path),
                 "/sys/devices/system/node/node%d/cpulist", n);
        FILE *f = fopen(path, "r");
        if (!f) break;
        if (fgets(buf, sizeof(buf), f) == NULL) {
          fclose(f);
          break;
        }
        fclose(f);

        /* cpulist is a range string: "0-23,48-71"  */
        char *tok = strtok(buf, ",\n");
        while (tok) {
            int lo, hi;
            if (sscanf(tok, "%d-%d", &lo, &hi) == 2) {
                if (cpu >= lo && cpu <= hi) return n;
            } else if (sscanf(tok, "%d", &lo) == 1) {
                if (cpu == lo) return n;
            }
            tok = strtok(NULL, ",\n");
        }
    }
    return -1;
}

static int bdf_cmp(const void *a, const void *b)
{
    const gpu_info_t *ga = a, *gb = b;
    if (ga->domain != gb->domain) return (int)(ga->domain - gb->domain);
    if (ga->bus    != gb->bus   ) return (int)(ga->bus    - gb->bus);
    if (ga->slot   != gb->slot  ) return (int)(ga->slot   - gb->slot);
    return (int)(ga->func - gb->func);
}

/* ── public API ────────────────────────────────────────────────────── */

/*
 * Scan /sys/bus/pci/devices/, collect GPU entries, sort by BDF,
 * then build the cpu -> gpu map.
 *
 * REQUIREMENT (CUDA): export CUDA_DEVICE_ORDER=PCI_BUS_ID
 * HIP already defaults to PCI bus order.
 *
 * Fortran visible as: build_gpu_affinity_map()
 */
void build_gpu_affinity_map_(void)
{
    DIR *dir = opendir("/sys/bus/pci/devices");
    if (!dir) {
        fprintf(stderr, "[gpu_affinity] cannot open /sys/bus/pci/devices\n");
        return;
    }

    ngpus = 0;
    struct dirent *e;
    while ((e = readdir(dir)) && ngpus < MAX_GPUS) {
        if (e->d_name[0] == '.') continue;

        char base[PATH_MAX], path[PATH_MAX];
        snprintf(base, sizeof(base), "/sys/bus/pci/devices/%s", e->d_name);

        /* Filter by vendor ─────────────────────────────────────────── */
        unsigned long vendor;
        path[0] = '\0';
        strncat(path, base, sizeof(path) - 1);
        strncat(path, "/vendor", sizeof(path) - strlen(path) - 1);

        if (read_hex(path, &vendor) < 0) continue;
        if (vendor != VENDOR_NVIDIA &&
            vendor != VENDOR_AMD    &&
            vendor != VENDOR_INTEL) continue;

        /* Filter by PCI class (drop the programming-interface byte) ─ */
        unsigned long pci_class;
        path[0] = '\0';
        strncat(path, base, sizeof(path) - 1);
        strncat(path, "/class", sizeof(path) - strlen(path) - 1);
        if (read_hex(path, &pci_class) < 0) continue;
        unsigned int cls = (unsigned int)(pci_class >> 8);
        if (cls != CLASS_VGA && cls != CLASS_3D && cls != CLASS_DISPLAY)
            continue;

        /* NUMA node (/sys value is -1 when not a NUMA system) ──────── */
        long numa = -1;
        path[0] = '\0';
        strncat(path, base, sizeof(path) - 1);
        strncat(path, "/numa_node", sizeof(path) - strlen(path) - 1);
        FILE *f = fopen(path, "r");
        if (f) {
          if (fscanf(f, "%ld", &numa) == EOF) numa = 0;
          fclose(f);
        }

        /* BDF from directory name:  dddd:bb:ss.f ───────────────────── */
        unsigned int dom, bus, slot, func;
        if (sscanf(e->d_name, "%x:%x:%x.%x",
                   &dom, &bus, &slot, &func) != 4) continue;

        gpus[ngpus++] = (gpu_info_t){
            .domain = dom, .bus = bus, .slot = slot, .func = func,
            .numa_node = (int)numa
        };
    }
    closedir(dir);

    /* Sort by BDF so device_index matches the runtime's PCI-order ─── */
    qsort(gpus, ngpus, sizeof(gpu_info_t), bdf_cmp);
    for (int g = 0; g < ngpus; g++) gpus[g].device_index = g;

    /* Build cpu -> closest gpu map ─────────────────────────────────── */
    memset(cpu_to_gpu, -1, sizeof(cpu_to_gpu));
    for (int c = 0; c < MAX_CPUS; c++) {
        int node = cpu_numa_node(c);
        if (node < 0) continue;
        /* First GPU on the same NUMA node wins.
           If multiple GPUs share one node, extend this to round-robin  */
        for (int g = 0; g < ngpus; g++) {
            if (gpus[g].numa_node == node) {
                cpu_to_gpu[c] = gpus[g].device_index;
                break;
            }
        }
    }

    affinity_ready = 1;

    printf("[gpu_affinity] %d GPU(s) found (sorted by PCI BDF):\n", ngpus);
    for (int g = 0; g < ngpus; g++)
        printf("  device %d  %04x:%02x:%02x.%x  NUMA node %d\n",
               gpus[g].device_index,
               gpus[g].domain, gpus[g].bus, gpus[g].slot, gpus[g].func,
               gpus[g].numa_node);
}


static int get_current_cpu(void)
{
    unsigned int cpu, node;
    if (syscall(SYS_getcpu, &cpu, &node, NULL) == 0)
        return (int)cpu;
    return 0;   /* fallback: should never happen on Linux >= 2.6.19 */
}

/*
 * Returns the GPU device index closest to the calling thread.
 * Must be called from within the OMP parallel region.
 */
int get_closest_gpu_(void)
{
    if (!affinity_ready) build_gpu_affinity_map_();
    int cpu = get_current_cpu();
    if (cpu < 0 || cpu >= MAX_CPUS) return 0;
    int gpu = cpu_to_gpu[cpu];
    return (gpu < 0) ? 0 : gpu;
}
