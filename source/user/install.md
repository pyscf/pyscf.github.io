# How to install PySCF 

:::{warning}

PySCF is not supported natively on Windows. You must use the Windows
Subsystem for Linux.

:::

## Install with pip

This is the recommended way to install PySCF for non-developers,
```bash
pip install --prefer-binary pyscf
```

The pip package provides a precompiled PySCF code (python wheel) that
works on almost all Linux systems, most macOS systems, and
the Windows Subsystem for Linux.

If you already have installed
PySCF via pip, you can upgrade it to the new version,
```bash
pip install --upgrade pyscf
```

:::{seealso}
Newly introduced features are available in the
[pyscf-forge](https://github.com/pyscf/pyscf-forge) package, which can be
installed with `pip`,
```bash
pip install pyscf-forge
```
Some other features are only maintained as [Extensions](user/extensions) to PySCF. 
:::

## Build from source with pip

If you're interested in a new feature that's not included in the latest release or you simply 
want the latest and greatest PySCF, you can build from source using pip,
```bash
pip install git+https://github.com/pyscf/pyscf
```

To install the features developed on a particular branch,
```bash
pip install git+https://github.com/pyscf/pyscf@<branch_name>
```

This install method compiles and links C extensions against the
libraries in your system. See [Build from source](#build-from-source) for a full 
list of prerequisites. If you would like to tune the CMake compilation 
parameters, you can set them with the environment variable `CMAKE_CONFIGURE_ARGS`, 
for example,
```bash
export CMAKE_CONFIGURE_ARGS="-DBUILD_MARCH_NATIVE=ON -DBLA_VENDOR=Intel10_64lp_seq"
```

See [CMake options](#cmake-options) for more details about CMake configuration.


## Build from source

You can manually install PySCF by building it from source.
Prerequisites for manual installation are
- C compiler
- C++ compiler (optional, but required for XCFun and some extensions)
- CMake >= 3.10
- Python >= 3.7
- Numpy >= 1.13
- Scipy >= 1.3
- h5py >= 2.7

You can download the latest version of PySCF (or the development
branch) from GitHub,
```bash
git clone https://github.com/pyscf/pyscf.git
cd pyscf
```
Next, you need to build the C extensions in `pyscf/lib`
```bash
cd pyscf/lib
mkdir build
cd build
cmake ..
make
```
This will automatically download the analytical GTO integral library
[libcint](https://github.com/sunqm/libcint.git) and the DFT exchange
correlation functional libraries [Libxc](https://libxc.gitlab.io/) 
and [XCFun](https://github.com/dftlibs/xcfun.git). Finally, to allow Python to
find the `pyscf` package, add the top-level `pyscf` directory (not the
`pyscf/pyscf` subdirectory) to `PYTHONPATH`.  For example, if `pyscf`
is installed in ``/opt``, you should update `PYTHONPATH` with something
like,
```bash
export PYTHONPATH=/opt/pyscf:$PYTHONPATH
```

To ensure that the installation was successful, you can use python to 
try to import and pring the PySCF version,
```bash
python -c "import pyscf; print(pyscf.__version__)"
```

See [CMake options](#cmake-options) for details about CMake configuration.


## Install with conda

If you have a [Conda](https://conda.io/docs/) (or
[Anaconda](https://www.continuum.io/downloads) environment, PySCF
package can be installed from the Conda cloud (for Linux and macOS
systems),
```bash
conda install -c pyscf pyscf
```

<!-- 
Extension modules are not available on the Conda cloud. They should be
installed either with pip or through the environment variable
`PYSCF_EXT_PATH` (see [Installing extensions](extensions)).
-->

## Install on Fedora

If you are running Fedora Linux, you can install PySCF as a
distribution package,
```bash
dnf install python3-pyscf
```

If you are running on an x86-64 platform, dnf should automatically
install the optimized integral library, qcint, instead of the
cross-platform libcint library.

<!--
Extension modules are not available in the Fedora package.
-->

## Docker image

The following command starts a container with a jupyter notebook
server that listens for HTTP connections on port 8888,
```
docker run -it -p 8888:8888 pyscf/pyscf:latest
```

Now you can visit ``https://localhost:8888`` with your browser to use
PySCF in the notebook.

Another way to use PySCF in a docker container is to start an Ipython
shell,
```bash
docker run -it pyscf/pyscf:latest start.sh ipython
```

## Advanced build options

### CMake options

A full build of PySCF may take a long time to finish, and the CMake options
listed below can be used to speed up compilation or omit packages that fail to compile
(e.g., `XCFun` may fail to build if a proper C++ compiler is not available).

:::{warning}

If both `-DENABLE_LIBXC=OFF` and `-DENABLE_XCFUN=OFF` are set, importing the
dft module will lead to an `ImportError`.

:::

| Flags                | Default | Comments
| -------------------- | ------- | ----------------------------------------------------------------- | 
| `ENABLE_LIBXC`       | ON      | Whether to use `Libxc` library in PySCF. If `-DENABLE_LIBXC=OFF` is appended to cmake command, `Libxc` will not be compiled. |
| `ENABLE_XCFUN`       | ON      | Whether to use `XCFun` library in PySCF. If `-DENABLE_XCFUN=OFF` is appended to cmake command, `XCFun` will not be compiled. |
| `BUILD_LIBXC`        |  ON     | Set it to `OFF` to skip compiling `Libxc`. The dft module still calls `Libxc` library by default. The dft module will be linked against the `Libxc` library from an earlier build. |
| `BUILD_XCFUN`        | ON      | Set it to `OFF` to skip compiling `XCFun`. The dft module will be linked against the `XCFun` library from an earlier build. |
| `BUILD_LIBCINT`      | ON      | Set it to `OFF` to skip compiling `libcint`. The integral library from an earlier build will be used. |
| `WITH_F12`           | ON      | Whether to compile the F12 integrals. |
| `DISABLE_DFT`        | OFF     | Set this flag to skip the entire dft module. Neither `Libxc` nor `XCFun` will be compiled. |
| `BUILD_MARCH_NATIVE` | OFF     | Whether to let the compiler optimize the code against CPU architecture. |

#### CMake config file

CMake options can be saved in a configuration file `pyscf/lib/cmake.arch.inc`.
The settings in this file will be automatically loaded and overwrite the
default CMake options during compilation.  For example, you can set
`CMAKE_C_FLAGS` in this file to include advanced compiler optimiztion flags,
```bash
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math -mtune=native -march=native")
```

Other settings, variables, and flags can also be set in this file,
```bash
  set(ENABLE_XCFUN Off)
  set(WITH_F12 Off)
```
Some examples of platform-specific configurations can be found in the
[cmake_user_inc_examples](https://github.com/pyscf/pyscf/tree/master/pyscf/lib/cmake_user_inc_examples)
directory.


### Environment variables and global configurations

| Env variable        | Comment                                              |
| ------------------- | ---------------------------------------------------- |
| `PYSCF_MAX_MEMORY`  | Maximum memory to use in MB                          |
| `PYSCF_TMPDIR`      | Directory for temporary files                        |
| `PYSCF_CONFIG_FILE` | File where various PySCF default settings are stored |
| `PYSCF_EXT_PATH`    | Path for finding external extensions                 |

`PYSCF_MAX_MEMORY` sets the default maximum memory in MB when creating
`Mole` (or `Cell`) object. It corresponds to the attribute
`max_memory``of Mole` (or `Cell`) object.

The environment variable `PYSCF_TMPDIR` controls which directory is
used to store intermediate files and temporary data when PySCF is run;
it is also commonly known as the scratch directory. If this
environment variable is not set, the system-wide temporary directory
`TMPDIR` will be used as the temp directory, instead. It is highly
recommended to set this variable to a directory with enough disk
space, as many quantum chemistry methods may consume a huge amount of
temporary storage space. It is equally important that the scratch
directory has fast i/o: for instance, using a network-mounted scratch
disk is often much slower than local disks.

`PYSCF_CONFIG_FILE` is a Python file that can be used to predefine and
override several default parameters in the program: you may already
have noticed statements like `getattr(__config__, "FOOBAR")` many
places in the source code. These global parameters are defined in
`PYSCF_CONFIG_FILE` and loaded when the pyscf module is imported.  By
default, this environment variable points to `~/.pyscf_conf.py`.

`PYSCF_EXT_PATH` allows PySCF to find any possible extension
packages. This is documented in detail in :ref:`installing_extproj`.

### Install without internet

In typical installations, all external libraries (libcint, Libxc, XCFun) are
downloaded and installed when the C extensions are compiled, thus requiring
internet access. In this section, we show how to install the external libraries
without internet access, assuming you have downloaded the libcint, Libxc, and
XCFun libraries,
```bash
git clone https://github.com/sunqm/libcint.git
tar czf libcint.tar.gz libcint
wget https://gitlab.com/libxc/libxc/-/archive/6.0.0/libxc-6.0.0.tar.gz
wget -O xcfun.tar.gz https://github.com/fishjojo/xcfun/archive/refs/tags/cmake-3.5.tar.gz
```

Assuming `/opt` is the place where these libraries will be installed, they
can be compiled,
```bash
tar xvzf libcint.tar.gz
cd libcint
mkdir build && cd build
cmake -DWITH_F12=1 -DWITH_RANGE_COULOMB=1 -DWITH_COULOMB_ERF=1 \
  -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
make && make install

tar xvzf libxc-6.0.0.tar.gz
cd libxc-6.0.0
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=1 \
  -DENABLE_FORTRAN=0 -DDISABLE_KXC=0 -DDISABLE_LXC=1 \
  -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
make && make install

tar xvzf xcfun.tar.gz
cd xcfun-cmake-3.5
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=1 -DXCFUN_MAX_ORDER=3 -DXCFUN_ENABLE_TESTS=0 \
  -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
make && make install
```

Next, compile PySCF,
```bash
cd pyscf/pyscf/lib
mkdir build && cd build
cmake -DBUILD_LIBCINT=0 -DBUILD_LIBXC=0 -DBUILD_XCFUN=0 -DCMAKE_INSTALL_PREFIX:PATH=/opt ..
make
```

Finally, update the ``PYTHONPATH`` environment variable so that the Python
interpreter can find your installation of PySCF.


### Using optimized BLAS

The default installation tries to find the BLAS libraries
automatically. This automated setup script may end up linking the code
to slow versions of BLAS libraries, like the reference NETLIB
implementation. Using an optimized linear algebra library like ATLAS,
BLIS, or OpenBLAS may speed up certain parts of PySCF by factors of 10 to 1000.

You can compile PySCF against BLAS libraries from other vendors to
improve performance. For example, the Intel Math Kernel Library (MKL)
can provide a 10x speedup in many modules,
```bash
cd pyscf/lib/build
cmake -DBLA_VENDOR=Intel10_64lp_seq ..
make
```

When linking the program to MKL, CMake may have problems finding the
correct MKL libraries for some versions of MKL.  Setting
``LD_LIBRARY_PATH`` to include the MKL dynamic libraries can sometimes
help, for example,
```bash
export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64:$LD_LIBRARY_PATH
```

If you are using Anaconda, you can link PySCF to the MKL library shipped with Anaconda,
```bash
export MKLROOT=/path/to/anaconda2
export LD_LIBRARY_PATH=$MKLROOT/lib:$LD_LIBRARY_PATH
cd pyscf/lib/build
cmake -DBLA_VENDOR=Intel10_64lp_seq ..
make
```

You can also link to other BLAS libraries by setting `BLA_VENDOR`,
e.g., `BLA_VENDOR=ATLAS`, `BLA_VENDOR=IBMESSL`,
`BLA_VENDOR=OpenBLAS` (requiring cmake-3.6).  Please refer to the
[CMake manual](http://www.cmake.org/cmake/help/v3.6/module/FindBLAS.html>) 
for more details on the use of the `FindBLAS` macro.

If setting the CMake `BLA_VENDOR` variable does not result in the
right BLAS library being chosen, you can specify the BLAS libraries to
use by hand by setting the `BLAS_LIBRARIES` CMake argument,
```bash
cmake -DBLAS_LIBRARIES=-lopenblaso ..
```

You can also hardcode the libraries you want to use in the file `lib/CMakeLists.txt`,
```bash
set(BLAS_LIBRARIES "${BLAS_LIBRARIES};/path/to/mkl/lib/intel64/libmkl_intel_lp64.so")
set(BLAS_LIBRARIES "${BLAS_LIBRARIES};/path/to/mkl/lib/intel64/libmkl_sequential.so")
set(BLAS_LIBRARIES "${BLAS_LIBRARIES};/path/to/mkl/lib/intel64/libmkl_core.so")
set(BLAS_LIBRARIES "${BLAS_LIBRARIES};/path/to/mkl/lib/intel64/libmkl_avx.so")
```

:::{note}
The MKL library may lead to an OSError at runtime,
```bash
OSError: ... mkl/lib/intel64/libmkl_avx.so: undefined symbol: ownLastTriangle_64fc
```
or
```bash
MKL FATAL ERROR: Cannot load libmkl_avx.so or libmkl_def.so
```
This can be solved by preloading the MKL core library,
```bash
export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_avx.so:$MKLROOT/lib/intel64/libmkl_core.so
```
:::

### Using the Qcint optimized integral library

The default integral library used by PySCF is 
[libcint](https://github.com/sunqm/libcint), which is implemented within a
model that maximizes its compatibility with various high-performance
computing systems. However, on x86-64 platforms, libcint has a more
efficient counterpart, [Qcint](https://github.com/sunqm/qcint), which is
heavily optimized with x86 SIMD instructions (AVX-512/AVX2/AVX/SSE3).
To replace the default libcint library with the qcint library, edit the
URL of the integral library in `lib/CMakeLists.txt` file,
```{code-block} bash
:caption: lib/CMakeLists.txt

ExternalProject_Add(libcint
    GIT_REPOSITORY
    https://github.com/sunqm/qcint.git
     ...
```
