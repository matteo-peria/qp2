# Package installation recipes for configure_managed.

function install_requested_packages() {
    for PACKAGE in ${PACKAGES} ; do

        if [[ ${PACKAGE} = openblas && -n ${EXTERNAL_OPENBLAS_PREFIX} ]] ; then
            echo "Using external OpenBLAS: ${EXTERNAL_OPENBLAS_PREFIX}"
            continue
        elif [[ ${PACKAGE} = hdf5 && -n ${EXTERNAL_HDF5_PREFIX} ]] ; then
            echo "Using external HDF5: ${EXTERNAL_HDF5_PREFIX}"
            continue
        elif [[ ( ${PACKAGE} = trexio || ${PACKAGE} = trexio-nohdf5 ) && \
                    -n ${EXTERNAL_TREXIO_PREFIX} ]] ; then
            echo "Using external TREXIO: ${EXTERNAL_TREXIO_PREFIX}"
            continue
        elif [[ ${PACKAGE} = qmckl && -n ${EXTERNAL_QMCKL_PREFIX} ]] ; then
            echo "Using external QMCkl: ${EXTERNAL_QMCKL_PREFIX}"
            continue
        fi

        if [[ ${PACKAGE} = ninja ]] ; then

            execute << EOF
              rm -f "\${QP_ROOT}"/bin/ninja
              tar -zxvf "\${QP_ROOT}"/external/qp2-dependencies/${ARCHITECTURE}/ninja.tar.gz
              mv ninja "\${QP_ROOT}"/bin/
EOF

        elif [[ ${PACKAGE} = cmake ]] ; then

            VERSION=3.31.10
            execute << EOF
              cd "\${QP_ROOT}"/external
              rm -f cmake-${VERSION}-linux-${ARCHITECTURE}.sh
              wget https://github.com/Kitware/CMake/releases/download/v${VERSION}/cmake-${VERSION}-linux-${ARCHITECTURE}.sh
              chmod +x cmake-${VERSION}-linux-${ARCHITECTURE}.sh
              ./cmake-${VERSION}-linux-${ARCHITECTURE}.sh --prefix=$QP_ROOT --skip-license
EOF

        elif [[ ${PACKAGE} = openblas ]] ; then

            VERSION=$OPENBLAS_VERSION
            PREFIX="$(package_prefix openblas "${VERSION}")"
            CC_ARGS="${CONFIG_CC#${CONFIG_CC%% *}} ${CONFIG_CFLAGS}"
            FC_ARGS="${CONFIG_FC#${CONFIG_FC%% *}} ${CONFIG_FCFLAGS}"
            execute << EOF
              cd "\${QP_ROOT}"/external
              wget https://github.com/OpenMathLib/OpenBLAS/releases/download/v${VERSION}/OpenBLAS-${VERSION}.tar.gz
              rm -rf OpenBLAS-${VERSION}
              tar -zxf OpenBLAS-${VERSION}.tar.gz && rm OpenBLAS-${VERSION}.tar.gz
              cd OpenBLAS-${VERSION}
              make -j 8 CC="${CONFIG_CC%% *}" FC="${CONFIG_FC%% *}" CFLAGS="${CC_ARGS}" FFLAGS="${FC_ARGS}" USE_OPENMP=1
              make PREFIX="${PREFIX}" install
EOF
            activate_dependency_prefixes

        elif [[ ${PACKAGE} = hdf5 ]] ; then

            VERSION=$HDF5_VERSION
            PREFIX="$(package_prefix hdf5 "${VERSION}")"
            execute << EOF
              cd "\${QP_ROOT}"/external
              wget https://github.com/HDFGroup/hdf5/releases/download/hdf5_${VERSION}/hdf5-${VERSION}.tar.gz
              rm -rf hdf5-${VERSION}
              tar -zxf hdf5-${VERSION}.tar.gz && rm hdf5-${VERSION}.tar.gz
              cd hdf5-${VERSION}
              ${CONFIGURE_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ./configure --prefix="${PREFIX}" --enable-fortran --enable-shared || exit 1
              (make -j 8 || make) && make check && make install
EOF
            activate_dependency_prefixes
            HDF5_PREFIX="${PREFIX}"
            H5CC="${PREFIX}/bin/h5cc"
            CONFIG_H5CC="${H5CC}"
            CONFIG_HDF5_PREFIX="${PREFIX}"
            TREXIO_H5CC_ENV="H5CC=\"${H5CC}\""
            TREXIO_HDF5_WITH_ARG="--with-hdf5=\"${PREFIX}\""

        elif [[ ${PACKAGE} = trexio-nohdf5 ]] ; then

            PREFIX="${MANAGED_TREXIO_NOHDF5_PREFIX}"
            if [[ ${TREXIO_SOURCE} = git ]] ; then
                require_git_source_tools TREXIO
                SOURCE_DIR="trexio-git-${TREXIO_RESOLVED_COMMIT:0:12}"
                execute << EOF
              cd "\${QP_ROOT}"/external
              rm -rf "${SOURCE_DIR}"
              git clone --no-checkout https://github.com/TREX-CoE/trexio.git "${SOURCE_DIR}"
              cd "${SOURCE_DIR}"
              git checkout --detach "${TREXIO_RESOLVED_COMMIT}"
              ./autogen.sh
              ${CONFIGURE_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ./configure --prefix="${PREFIX}" --without-hdf5 || exit 1
              (make -j 8 || make) && make check && make -j 8 install
EOF
            else
                VERSION=$TREXIO_VERSION
                execute << EOF
              cd "\${QP_ROOT}"/external
              wget https://github.com/TREX-CoE/trexio/releases/download/v${VERSION}/trexio-${VERSION}.tar.gz
              rm -rf trexio-${VERSION}
              tar -zxf trexio-${VERSION}.tar.gz && rm trexio-${VERSION}.tar.gz
              cd trexio-${VERSION}
              ${CONFIGURE_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ./configure --prefix="${PREFIX}" --without-hdf5 || exit 1
              (make -j 8 || make) && make check && make -j 8 install
              tar -zxvf "\${QP_ROOT}"/external/qp2-dependencies/${ARCHITECTURE}/ninja.tar.gz
              mv ninja "\${QP_ROOT}"/bin/
EOF
            fi
            activate_dependency_prefixes
        elif [[ ${PACKAGE} = trexio ]] ; then

            PREFIX="${MANAGED_TREXIO_PREFIX}"
            if [[ ${TREXIO_SOURCE} = git ]] ; then
                require_git_source_tools TREXIO
                SOURCE_DIR="trexio-git-${TREXIO_RESOLVED_COMMIT:0:12}"
                execute << EOF
              cd "\${QP_ROOT}"/external
              rm -rf "${SOURCE_DIR}"
              git clone --no-checkout https://github.com/TREX-CoE/trexio.git "${SOURCE_DIR}"
              cd "${SOURCE_DIR}"
              git checkout --detach "${TREXIO_RESOLVED_COMMIT}"
              ./autogen.sh
              ${TREXIO_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ${TREXIO_H5CC_ENV} ${TREXIO_HDF5_ENV} ./configure --prefix="${PREFIX}" ${TREXIO_HDF5_WITH_ARG} || exit 1
              (make -j 8 || make) && make check && make -j 8 install
EOF
            else
                VERSION=$TREXIO_VERSION
                execute << EOF
              cd "\${QP_ROOT}"/external
              wget https://github.com/TREX-CoE/trexio/releases/download/v${VERSION}/trexio-${VERSION}.tar.gz
              rm -rf trexio-${VERSION}
              tar -zxf trexio-${VERSION}.tar.gz && rm trexio-${VERSION}.tar.gz
              cd trexio-${VERSION}
              ${TREXIO_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ${TREXIO_H5CC_ENV} ${TREXIO_HDF5_ENV} ./configure --prefix="${PREFIX}" ${TREXIO_HDF5_WITH_ARG} || exit 1
              (make -j 8 || make) && make check && make -j 8 install
EOF
            fi
            activate_dependency_prefixes
        elif [[ ${PACKAGE} = qmckl ]] ; then

            PREFIX="${MANAGED_QMCKL_PREFIX}"
            if [[ ${QMCKL_SOURCE} = git ]] ; then
                require_git_source_tools QMCkl
                SOURCE_DIR="qmckl-git-${QMCKL_RESOLVED_COMMIT:0:12}"
                execute << EOF
              cd "\${QP_ROOT}"/external
              rm -rf "${SOURCE_DIR}"
              git clone --no-checkout https://github.com/TREX-CoE/qmckl.git "${SOURCE_DIR}"
              cd "${SOURCE_DIR}"
              git checkout --detach "${QMCKL_RESOLVED_COMMIT}"
              ./autogen.sh
              ${CONFIGURE_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_F77_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ${CONFIGURE_FFLAGS_ENV} ${CONFIGURE_BLAS_ENV} ${CONFIGURE_LAPACK_ENV} ./configure --prefix="${PREFIX}" --enable-hpc --disable-doc || exit 1
              (make -j 8 || make) && make check && make install
EOF
            else
                VERSION=$QMCKL_VERSION
                execute << EOF
              cd "\${QP_ROOT}"/external
              wget https://github.com/TREX-CoE/qmckl/releases/download/v${VERSION}/qmckl-${VERSION}.tar.gz
              rm -rf qmckl-${VERSION}
              tar -zxf qmckl-${VERSION}.tar.gz && rm qmckl-${VERSION}.tar.gz
              cd qmckl-${VERSION}
              ${CONFIGURE_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_F77_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ${CONFIGURE_FFLAGS_ENV} ${CONFIGURE_BLAS_ENV} ${CONFIGURE_LAPACK_ENV} ./configure --prefix="${PREFIX}" --enable-hpc --disable-doc || exit 1
              (make -j 8 || make) && make check && make install
EOF
            fi
            activate_dependency_prefixes
        elif [[ ${PACKAGE} = gmp ]] ; then

            VERSION=6.1.2
            PREFIX="$(package_prefix gmp "${VERSION}")"
            GMP_TARGET=""
            [[ ${ARCHITECTURE} = arm64Darwin ]] && GMP_TARGET="--build=aarch64-apple-darwin --host=aarch64-apple-darwin --disable-assembly"
            execute << EOF
              cd "\${QP_ROOT}"/external
              rm -rf gmp-${VERSION}
              tar --bzip2 --extract --file qp2-dependencies/gmp-6.1.2.tar.bz2
              cd gmp-${VERSION}
              ${CONFIGURE_CC_ENV} ${CONFIGURE_CFLAGS_ENV} ./configure --prefix="${PREFIX}" ${GMP_TARGET} && make -j 8
              make -j 8 install
EOF
            activate_dependency_prefixes

        elif [[ ${PACKAGE} = zeromq ]] ; then

            [[ -n ${CONFIG_CXX} ]] || error \
                                          "CXX is not defined in ${CONFIG:-the selected compiler environment}." \
                                          "Add CXX to the selected configuration file."
            VERSION=$ZEROMQ_VERSION
            PREFIX="$(package_prefix zeromq "${VERSION}")"
            execute << EOF
              cd "\${QP_ROOT}"/external
              wget https://github.com/zeromq/libzmq/releases/download/v${VERSION}/zeromq-${VERSION}.tar.gz
              rm -rf zeromq-${VERSION}
              tar -zxf zeromq-${VERSION}.tar.gz && rm zeromq-${VERSION}.tar.gz
              cd zeromq-${VERSION}
              [[ "${SYSTEM}" = "Darwin" ]] && ./autogen.sh
              ${CONFIGURE_CC_ENV} CXX="${CONFIG_CXX}" ${CONFIGURE_CFLAGS_ENV} ./configure --prefix="${PREFIX}" --without-libsodium --enable-libunwind=no
              make -j 8
              make install
EOF
            activate_dependency_prefixes


        elif [[ ${PACKAGE} = f77zmq ]] ; then

            VERSION=$F77ZMQ_VERSION
            PREFIX="$(package_prefix f77zmq "${VERSION}")"
            ZMQ_PREFIX="$(package_prefix zeromq "${ZEROMQ_VERSION}")"
            execute << EOF
              cd "\${QP_ROOT}"/external
              rm -rf f77-zmq-${VERSION}
              tar --gunzip --extract --file qp2-dependencies/f77-zmq-${VERSION}.tar.gz
              cd f77-zmq-${VERSION}
              ${CONFIGURE_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_F77_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ${CONFIGURE_FFLAGS_ENV} ./configure --prefix="${PREFIX}"
              export ZMQ_H="${ZMQ_PREFIX}/include/zmq.h"
              make && make check && make install
EOF
            activate_dependency_prefixes


        elif [[ ${PACKAGE} = ocaml ]] ; then

            execute <<EOF
            rm -rf "${QP_ROOT}"/external/opampack
            cd "${QP_ROOT}"/external/
            tar --gunzip --extract --file qp2-dependencies/${ARCHITECTURE}/opampack.tar.gz
            cd "${QP_ROOT}"/external/opampack
            ./install.sh
            export OPAMROOT="${QP_ROOT}"/external/opampack/opamroot
            eval \$("${QP_ROOT}"/external/opampack/opam env)
EOF
            source "${QP_ROOT}"/etc/ocaml.rc

        elif [[ ${PACKAGE} = bse ]] ; then

            execute << EOF
              cd "\${QP_ROOT}"/external
              tar --gunzip --extract --file qp2-dependencies/bse-v0.8.11.tar.gz
              python3 -m pip install -e basis_set_exchange-*
EOF

        elif [[ ${PACKAGE} = zlib ]] ; then

            VERSION=$ZLIB_VERSION
            PREFIX="$(package_prefix zlib "${VERSION}")"
            execute << EOF
              cd "\${QP_ROOT}"/external
              wget https://zlib.net/fossils/zlib-${VERSION}.tar.gz
              rm -rf zlib-${VERSION}
              tar -zxf zlib-${VERSION}.tar.gz && rm zlib-${VERSION}.tar.gz
              cd zlib-${VERSION}/
              ${CONFIGURE_CC_ENV} ${CONFIGURE_CFLAGS_ENV} ./configure --prefix="${PREFIX}" && \
              make && make install
EOF
            activate_dependency_prefixes


        elif [[ ${PACKAGE} = docopt ]] ; then

            execute << EOF
              cd "\${QP_ROOT}"/external
              tar --gunzip --extract --file qp2-dependencies/docopt-0.6.2.tar.gz
              mv docopt-*/docopt.py "\${QP_ROOT}/external/Python"
EOF


        elif [[ ${PACKAGE} = resultsFile ]] ; then

            execute << EOF
              cd "\${QP_ROOT}"/external
              tar --gunzip --extract --file qp2-dependencies/resultsFile-v2.3.tar.gz
              rm -rf "\${QP_ROOT}/external/Python/resultsFile"
              mv resultsFile-*/resultsFile "\${QP_ROOT}/external/Python/"
EOF

        elif [[ ${PACKAGE} = bats ]] ; then

            execute << EOF
              cd "\${QP_ROOT}"/external
              tar -zxf qp2-dependencies/bats-v1.7.0.tar.gz
              ( cd bats-core-1.7.0/ ; ./install.sh \${QP_ROOT})
EOF

        elif [[ ${PACKAGE} = armadillo ]] ; then

            execute << EOF
              cd "\${QP_ROOT}"/external
              cd armadillo-code
              #TODO
EOF

        elif [[ ${PACKAGE} = openorb ]] ; then

            execute << EOF
              cd "\${QP_ROOT}"/external
              cd OpenOrbitalOptimizer
              #TODO
EOF
        elif [[ ${PACKAGE} = simple_gpu ]] ; then

            VERSION=$SIMPLE_GPU_VERSION
            PREFIX="$(package_prefix simple_gpu "${VERSION}")"
            execute << EOF
              cd "\${QP_ROOT}"/external
              wget https://github.com/scemama/simple_gpu/releases/download/v${VERSION}/simple_gpu-${VERSION}.tar.gz
              rm -rf simple_gpu-${VERSION}
              tar -zxf simple_gpu-${VERSION}.tar.gz && rm simple_gpu-${VERSION}.tar.gz
              cd simple_gpu-${VERSION}
              ${CONFIGURE_CC_ENV} ${CONFIGURE_FC_ENV} ${CONFIGURE_F77_ENV} ${CONFIGURE_CFLAGS_ENV} ${CONFIGURE_FCFLAGS_ENV} ${CONFIGURE_FFLAGS_ENV} ${CONFIGURE_BLAS_ENV} ${CONFIGURE_LAPACK_ENV} ./configure --prefix="${PREFIX}" || exit 1
              (make -j 8 || make) && make check && make install
EOF
            activate_dependency_prefixes
        else
            error "${PACKAGE} unknown."
            fail
        fi


    done
}
