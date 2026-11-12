# Dependency validation helpers for configure_managed.

function fail() {
    echo "You can try to install it using the -i option."
    echo "Please refer to INSTALL.rst to install the missing dependencies."
    exit -1
}

function success() {
    echo ""
    echo "Configuration successful."
    exit 0
}

function not_found() {
    echo 'not_found'
}

function find_exe() {
    command -v "$1" 2> /dev/null || not_found
}

function find_python_lib() {
    python3 -c "import $1" &> /dev/null && echo "$1" || not_found
}

function check_c_program() {
    eval "${CONFIG_CC} ${CONFIG_CFLAGS} \
\"${QP_ROOT}\"/external/tmp.c -o \"${QP_ROOT}\"/external/tmp.exe \
${CONFIG_LIB} $*" 2> /dev/null
    if [[ -x "${QP_ROOT}"/external/tmp.exe ]] && \
           "${QP_ROOT}"/external/tmp.exe 2> /dev/null ; then
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.c
        echo "found"
    else
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.c
        not_found
    fi
}

function check_zeromq() {
    cat > "${QP_ROOT}"/external/tmp.c <<'EOF'
#include <zmq.h>
int main(void) {
  int major, minor, patch;
  zmq_version(&major, &minor, &patch);
  return major < 0;
}
EOF
    check_c_program -lzmq
}

function check_gmp() {
    cat > "${QP_ROOT}"/external/tmp.c <<'EOF'
#include <gmp.h>
int main(void) {
  mpz_t value;
  mpz_init_set_ui(value, 42ul);
  int result = mpz_cmp_ui(value, 42ul);
  mpz_clear(value);
  return result != 0;
}
EOF
    check_c_program -lgmp
}

function check_zlib() {
    cat > "${QP_ROOT}"/external/tmp.c <<'EOF'
#include <zlib.h>
int main(void) {
  return zlibVersion() == 0;
}
EOF
    check_c_program -lz
}

function check_simple_gpu() {
    cat > "${QP_ROOT}"/external/tmp.c <<'EOF'
#include <stddef.h>
#include <simple_gpu.h>
int main(void) {
  return gpu_ndevices() < 0;
}
EOF
    check_c_program -lsimple_gpu-cpu
}

function check_trexio() {
    cat > "${QP_ROOT}"/external/tmp.f90 <<'EOF'
program main
  use trexio
  if (.not. trexio_has_back_end(TREXIO_TEXT)) stop 1
end program main
EOF
    eval "${CONFIG_FC} ${CONFIG_FCFLAGS} \
\"${QP_ROOT}\"/src/trexio/trexio_module.F90 \
\"${QP_ROOT}\"/external/tmp.f90 -o \"${QP_ROOT}\"/external/tmp.exe \
${CONFIG_LIB} -ltrexio" 2> /dev/null
    if [[ -x "${QP_ROOT}"/external/tmp.exe ]] && \
           "${QP_ROOT}"/external/tmp.exe 2> /dev/null ; then
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90 trexio.mod
        echo "found"
    else
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90 trexio.mod
        not_found
    fi
}

function check_qmckl() {
    cat > "${QP_ROOT}"/external/tmp.f90 <<'EOF'
program main
  use qmckl
  integer(qmckl_context) :: context
  integer(qmckl_exit_code) :: rc
  context = qmckl_context_create()
  if (context == QMCKL_NULL_CONTEXT) stop 1
  rc = qmckl_context_destroy(context)
  if (rc /= QMCKL_SUCCESS) stop 2
end program main
EOF
    eval "${CONFIG_FC} ${CONFIG_FCFLAGS} \
\"${QP_ROOT}\"/plugins/local/qmckl/qmckl.F90 \
\"${QP_ROOT}\"/external/tmp.f90 -o \"${QP_ROOT}\"/external/tmp.exe \
${CONFIG_LIB} ${BLAS_LIBS:-${CONFIG_LAPACK_LIB}} \
-lqmckl -ltrexio" 2> /dev/null
    if [[ -x "${QP_ROOT}"/external/tmp.exe ]] && \
           "${QP_ROOT}"/external/tmp.exe 2> /dev/null ; then
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90 qmckl.mod \
           qmckl_constants.mod qmckl_verificarlo_f.mod
        echo "found"
    else
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90 qmckl.mod \
           qmckl_constants.mod qmckl_verificarlo_f.mod
        not_found
    fi
}

function check_f77zmq() {
    cat > "${QP_ROOT}"/external/tmp.f90 <<'EOF'
program main
  use f77_zmq
  integer :: major, minor, patch, rc
  rc = f77_zmq_version(major, minor, patch)
  if (rc /= 0) stop 1
end program main
EOF
    eval "${CONFIG_FC} ${CONFIG_FCFLAGS} -I\"${QP_ROOT}/src/zmq\" \
\"${QP_ROOT}/src/zmq/f77_zmq_module.F90\" \
\"${QP_ROOT}\"/external/tmp.f90 -o \"${QP_ROOT}\"/external/tmp.exe \
${CONFIG_LIB} -lf77zmq -lzmq -lpthread" 2> /dev/null
    if [[ -x "${QP_ROOT}"/external/tmp.exe ]] && \
           "${QP_ROOT}"/external/tmp.exe 2> /dev/null ; then
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90 f77_zmq.mod
        echo "found"
    else
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90 f77_zmq.mod
        not_found
    fi
}

function check_blas_lapack() {
    local linear_algebra_libs
    linear_algebra_libs="${LAPACK_LIBS:-${CONFIG_LAPACK_LIB}}"
    if [[ -n ${BLAS_LIBS:-} && \
          ${BLAS_LIBS} != "${linear_algebra_libs}" ]] ; then
        linear_algebra_libs="${linear_algebra_libs} ${BLAS_LIBS}"
    fi
    cat > "${QP_ROOT}"/external/tmp.f90 <<'EOF'
program main
  implicit none
  double precision :: a(2,2), b(2,2), c(2,2), work(8), w(2)
  integer :: info
  a = reshape((/ 2.d0, 0.d0, 0.d0, 3.d0 /), shape(a))
  b = 0.d0
  b(1,1) = 1.d0
  b(2,2) = 1.d0
  c = 0.d0
  call dgemm('N', 'N', 2, 2, 2, 1.d0, a, 2, b, 2, 0.d0, c, 2)
  if (maxval(abs(c - a)) > 1.d-12) stop 1
  call dsyev('N', 'U', 2, a, 2, w, work, size(work), info)
  if (info /= 0) stop 2
  if (maxval(abs(w - (/ 2.d0, 3.d0 /))) > 1.d-12) stop 3
end program main
EOF
    eval "${CONFIG_FC} ${CONFIG_FCFLAGS} \
\"${QP_ROOT}\"/external/tmp.f90 -o \"${QP_ROOT}\"/external/tmp.exe \
${CONFIG_LIB} ${linear_algebra_libs}" 2> /dev/null
    if [[ -x "${QP_ROOT}"/external/tmp.exe ]] && \
           "${QP_ROOT}"/external/tmp.exe 2> /dev/null ; then
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90
        echo "found"
    else
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90
        not_found
    fi
}

function check_hdf5() {
    local h5fc h5fc_show h5fc_flags
    h5fc=$(find_hdf5_fortran_wrapper)
    if [[ ${h5fc} = $(not_found) ]] ; then
        not_found
        return
    fi
    h5fc_show="$(${h5fc} -show 2> /dev/null)"
    if [[ -z ${h5fc_show} || ${h5fc_show} != *" "* ]] ; then
        not_found
        return
    fi
    h5fc_flags="${h5fc_show#* }"
    cat > "${QP_ROOT}"/external/tmp.f90 <<'EOF'
program main
  use hdf5
  integer :: ierr
  call h5open_f(ierr)
  if (ierr /= 0) stop 1
  call h5close_f(ierr)
  if (ierr /= 0) stop 2
end program main
EOF
    eval "${CONFIG_FC} ${CONFIG_FCFLAGS} \
\"${QP_ROOT}\"/external/tmp.f90 -o \"${QP_ROOT}\"/external/tmp.exe \
${h5fc_flags}" 2> /dev/null
    if [[ -x "${QP_ROOT}"/external/tmp.exe ]] && \
           "${QP_ROOT}"/external/tmp.exe 2> /dev/null ; then
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90
        echo "found"
    else
        rm -f "${QP_ROOT}"/external/tmp.exe \
           "${QP_ROOT}"/external/tmp.f90
        not_found
    fi
}

function find_dir() {
    if [[ -d $1 ]] ; then
        echo "$1"
    else
        not_found
    fi
}
