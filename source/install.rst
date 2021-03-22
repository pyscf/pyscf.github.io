.. _installing:

Installation
************

Installation with pip
=====================

This is the recommended way to install PySCF::

  $ pip install pyscf

The pip package provides a precompiled PySCF code (python wheel) which
works on almost all Linux systems, and most of Mac OS X systems, and
the Ubuntu subsystems on Windows 10. If you already have installed
PySCF via pip, you can upgrade it to the new version with::

  $ pip install --upgrade pyscf

Since PySCF version 2.0, some modules are developed independently; see
:ref:`installing_extproj`. Individual extension modules (for example
the geometry optimization module) can be installed using pip's extra
dependency mechanism::

  $ pip install pyscf[geomopt]

All extension modules can be installed with::

  $ pip install pyscf[all]

The extension modules can be found in `https://github.com/pyscf` (see
also :ref:`installing_extproj`).

Installing the latest code on GitHub with pip
---------------------------------------------
The latest code on github can be installed with::

  $ pip install git+https://github.com/pyscf/pyscf

To install the features developed on a particular branch use::

  $ pip install git+https://github.com/pyscf/pyscf@<branch_name>

This install method compiles and links C extensions against the
libraries in your system. It requires CMake, a BLAS library and the
GCC compiler (more details of the prerequisites can be found in
:ref:`compile_c_extensions`). The C extensions are compiled with the
default settings specified in the `CMakeLists.txt` file. If you would
like to tune the CMake compilation parameters, you can set them with
the environment variable `CMAKE_CONFIGURE_ARGS`. The contents of this
environment variable will be passed in full to CMake. For example, if
you have multiple BLAS libraries available in your system, and MKL is
the one you would like to use, you can accomplish this by specifying
the environment variable (see also :ref:`installing_blas`) as::

  $ export CMAKE_CONFIGURE_ARGS="-DBLA_VENDOR=Intel10_64lp_seq"

To install the latest versions of the extension modules from GitHub,
you can place the GitHub repo url with a `git+` prefix in the argument
list of the pip command::

  $ pip install git+https://github.com/pyscf/geomopt


Installation on Fedora
======================

If you are running Fedora Linux, you can install PySCF as a
distribution package::

  # dnf install python3-pyscf

If you are running on an X86-64 platform, dnf should automatically
install the optimized integral library, qcint, instead of the
cross-platform libcint library.

Extension modules are not available in the Fedora package.

Installation with conda
=======================

If you have a `Conda <https://conda.io/docs/>`_ (or `Anaconda
<https://www.continuum.io/downloads#linux>`_) environment, PySCF
package can be installed from the Conda cloud::

  $ conda install -c pyscf pyscf

Extension modules are not available on the Conda cloud. They should be
installed either with pip, or through the environment variable
`PYSCF_EXT_PATH` (see the section :ref:`installing_extproj`).


PySCF docker image
==================

The following command starts a container with the jupyter notebook
server that listens for HTTP connections on port 8888::

  $ docker run -it -p 8888:8888 pyscf/pyscf:latest

Now, you can visit ``https://localhost:8888`` with your browser to use
PySCF in the notebook.

Another way to use PySCF in a docker container is to start an Ipython
shell::

  $ docker run -it pyscf/pyscf:latest start.sh ipython


.. _compile_c_extensions:
Compiling from source code
==========================

Prerequisites for manual install are

* CMake >= 3.10
* Python >= 3.6
* Numpy >= 1.13
* Scipy >= 0.19
* h5py >= 2.7

You can download the latest version of PySCF (or the development
branch) from github::

  $ git clone https://github.com/pyscf/pyscf.git
  $ cd pyscf
  $ git checkout dev  # optional if you'd like to try out the development branch

Next, you need to build the C extensions in :file:`pyscf/lib`::

  $ cd pyscf/lib
  $ mkdir build
  $ cd build
  $ cmake ..
  $ make

This will automatically download the analytical GTO integral library
`libcint <https://github.com/sunqm/libcint.git>`_ and the DFT exchange
correlation functional libraries `Libxc
<http://www.tddft.org/programs/Libxc>`_ and `XCFun
<https://github.com/dftlibs/xcfun.git>`_.  Finally, to allow Python to
find the `pyscf` package, add the top-level `pyscf` directory (not the
`pyscf/pyscf` subdirectory) to `PYTHONPATH`.  For example, if `pyscf`
is installed in ``/opt``, you adjust `PYTHONPATH` with something
like::

  export PYTHONPATH=/opt/pyscf:$PYTHONPATH

To ensure that the installation was successful, you can start a Python
shell, and type::

  >>> import pyscf

For Mac OS X/macOS, you may get an import error if your OS X/macOS
version is 10.11 or newer::

    OSError: dlopen(xxx/pyscf/pyscf/lib/libcgto.dylib, 6): Library not loaded: libcint.3.0.dylib
    Referenced from: xxx/pyscf/pyscf/lib/libcgto.dylib
    Reason: unsafe use of relative rpath libcint.3.0.dylib in xxx/pyscf/pyscf/lib/libcgto.dylib with restricted binary

This is caused by the incorrect RPATH.  The script
``pyscf/lib/_runme_to_fix_dylib_osx10.11.sh`` in the ``pyscf/lib``
directory can be used to fix this problem::

    cd pyscf/lib
    sh _runme_to_fix_dylib_osx10.11.sh

.. note::

  RPATH has been built in the dynamic library.  This may cause library
  loading error on some systems.  You can run
  ``pyscf/lib/_runme_to_remove_rpath.sh`` to remove the rpath code
  from the library head.  Another workaround is to set
  ``-DCMAKE_SKIP_RPATH=1`` and ``-DCMAKE_MACOSX_RPATH=0`` in the CMake
  command line.  When the RPATH was removed, you need to add
  ``pyscf/lib`` and ``pyscf/lib/deps/lib`` in ``LD_LIBRARY_PATH``.



Environment variables and global configures
===========================================

----------------------- ---------------------------------------------------------
Env variable            Comments
----------------------- ---------------------------------------------------------
`PYSCF_MAX_MEMORY`      Maximum memory to use in MB
`PYSCF_TMPDIR`          Directory for temporary files
`PYSCF_CONFIG_FILE`     File where various PySCF default settings are stored
`PYSCF_EXT_PATH`        Path for finding external extensions
----------------------- ---------------------------------------------------------

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


.. _installing_wo_network:
Installation without network
============================

In the usual case, all external libraries (libcint, libxc, xcfun) are
downloaded and installed when the C extensions are compiled, thus
requiring network access. In this section, we show how to install the
external libraries without accessing to network. First, you need to
install the libcint, Libxc, and XCFun libraries::

    $ git clone https://github.com/sunqm/libcint.git
    $ tar czf libcint.tar.gz libcint

    $ wget https://gitlab.com/libxc/libxc/-/archive/4.3.4/libxc-4.3.4.tar.gz

    $ git clone https://github.com/sunqm/xcfun.git
    $ tar czf xcfun.tar.gz xcfun

Assuming ``/opt`` is the place where these libraries will be installed, these
packages should be compiled with the flags::

    $ tar xvzf libcint.tar.gz
    $ cd libcint
    $ mkdir build && cd build
    $ cmake -DWITH_F12=1 -DWITH_RANGE_COULOMB=1 -DWITH_COULOMB_ERF=1 \
        -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
    $ make && make install

    $ tar xvzf libxc-4.3.4.tar.gz
    $ cd libxc-4.3.4
    $ mkdir build && cd build
    $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=1 \
        -DENABLE_FORTRAN=0 -DDISABLE_KXC=0 -DDISABLE_LXC=1 \
        -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
    $ make && make install

    $ tar xvzf xcfun.tar.gz
    $ cd xcfun
    $ mkdir build && cd build
    $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=1 -DXC_MAX_ORDER=3 -DXCFUN_ENABLE_TESTS=0 \
        -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
    $ make && make install

Next, you can compile PySCF::

    $ cd pyscf/pyscf/lib
    $ mkdir build && cd build
    $ cmake -DBUILD_LIBCINT=0 -DBUILD_LIBXC=0 -DBUILD_XCFUN=0 -DCMAKE_INSTALL_PREFIX:PATH=/opt ..
    $ make

Finally, you should update the ``PYTHONPATH`` environment variable so
that the Python interpreter can find your installation of PySCF.


.. _installing_blas:
Using optimized BLAS
====================

The default installation tries to find the BLAS libraries
automatically. This automated setup script may end up linking the code
to slow versions of BLAS libraries, like the reference NETLIB
implementation. Using an optimized linear algebra library like ATLAS,
BLIS or OpenBLAS may, however, speed up certain parts of PySCF by
orders of magnitudes; speedups by a factor of 1000x over the reference
implementation are not uncommon.

You can compile PySCF against BLAS libraries from other vendors to
improve performance. For example, the Intel Math Kernel Library (MKL)
can provide a 10x speedup in many modules::

  $ cd pyscf/lib/build
  $ cmake -DBLA_VENDOR=Intel10_64lp_seq ..
  $ make

When linking the program to MKL, CMake may have problems to find the
correct MKL libraries for some versions of MKL.  Setting
``LD_LIBRARY_PATH`` to include the MKL dynamic libraries can sometimes
help, e.g.::

  export LD_LIBRARY_PATH=/opt/intel/compilers_and_libraries_2018/linux/mkl/lib/intel64:$LD_LIBRARY_PATH

If you are using Anaconda as your Python-side platform, you can link
PySCF to the MKL library shipped with Anaconda::

  $ export MKLROOT=/path/to/anaconda2
  $ export LD_LIBRARY_PATH=$MKLROOT/lib:$LD_LIBRARY_PATH
  $ cd pyscf/lib/build
  $ cmake -DBLA_VENDOR=Intel10_64lp_seq ..
  $ make

You can also link to other BLAS libraries by setting ``BLA_VENDOR``,
eg ``BLA_VENDOR=ATLAS``, ``BLA_VENDOR=IBMESSL``,
``BLA_VENDOR=OpenBLAS`` (requiring cmake-3.6).  Please refer to the
`cmake manual
<http://www.cmake.org/cmake/help/v3.6/module/FindBLAS.html>`_ for more
details on the use of the ``FindBLAS`` macro.

If setting the CMake ``BLA_VENDOR`` variable does not result in the
right BLAS library being chosen, you can specify the BLAS libraries to
use by hand by setting the ``BLAS_LIBRARIES`` CMake argument::
  $ cmake -DBLAS_LIBRARIES=-lopenblaso ..

You can also hardcode the libraries you want to use in
:file:`lib/CMakeLists.txt`::

  set(BLAS_LIBRARIES "${BLAS_LIBRARIES};/path/to/mkl/lib/intel64/libmkl_intel_lp64.so")
  set(BLAS_LIBRARIES "${BLAS_LIBRARIES};/path/to/mkl/lib/intel64/libmkl_sequential.so")
  set(BLAS_LIBRARIES "${BLAS_LIBRARIES};/path/to/mkl/lib/intel64/libmkl_core.so")
  set(BLAS_LIBRARIES "${BLAS_LIBRARIES};/path/to/mkl/lib/intel64/libmkl_avx.so")

.. note::
  MKL library may lead to an OSError at runtime:
  ``OSError: ... mkl/lib/intel64/libmkl_avx.so: undefined symbol: ownLastTriangle_64fc``
  or ``MKL FATAL ERROR: Cannot load libmkl_avx.so or libmkl_def.so.``.
  It can be solved by preloading MKL core library with:
  ``export LD_PRELOAD=$MKLROOT/lib/intel64/libmkl_avx.so:$MKLROOT/lib/intel64/libmkl_core.so``


.. _installing_qcint:
Using optimized integral library
================================

The default integral library used by PySCF is libcint
(https://github.com/sunqm/libcint), which is implemented within a
model that maximizes its compatibility with various high performance
computer systems. On X86-64 platforms, however, libcint has a more
efficient counterpart, Qcint (https://github.com/sunqm/qcint) which is
heavily optimized with X86 SIMD instructions (AVX-512/AVX2/AVX/SSE3).
To replace the default libcint library with qcint library, edit the
URL of the integral library in lib/CMakeLists.txt file::

  ExternalProject_Add(libcint
     GIT_REPOSITORY
     https://github.com/sunqm/qcint.git
     ...


Cmake config file
=================

CMake options can be saved in a configuration file
``pyscf/lib/cmake.arch.inc``.  The settings in this file will be
automatically loaded and overwrite the default CMake options during
compilation.  For example, you can set ``CMAKE_C_FLAGS`` in this file
to include advanced compiler optimization flags::

  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -ffast-math -mtune=native -march=native")

Other settings, variables, and flags can also be set in this file::

  set(ENABLE_XCFUN Off)
  set(WITH_F12 Off)

Some examples of platform-specific configurations can be found in
directory ``pyscf/lib/cmake_arch_config``.


.. _installing_extproj:
Extension modules
=================

As of PySCF-2.0, some modules have been moved from the main code trunk
to extension projects hosted at https://github.com/pyscf.

------------------- ---------------------------------------------------------
Project             URL
------------------- ---------------------------------------------------------
cornell_shci        https://github.com/pyscf/cornell_shci
dftd3               https://github.com/pyscf/dftd3
dmrgscf             https://github.com/pyscf/dmrgscf
doci                https://github.com/pyscf/doci
fciqmcscf           https://github.com/pyscf/fciqmcscf
icmpspt             https://github.com/pyscf/icmpspt
mbd                 https://github.com/pyscf/mbd
naive_hci           https://github.com/pyscf/naive_hci
nao                 https://github.com/pyscf/nao
rt                  https://github.com/pyscf/rt
semiempirical       https://github.com/pyscf/semiempirical
shciscf             https://github.com/pyscf/shciscf
zquatev             https://github.com/sunqm/zquatev
tblis               https://github.com/pyscf/pyscf-tblis
------------------- ---------------------------------------------------------

Based on the technique of namespace packages specified in `PEP 420
<https://www.python.org/dev/peps/pep-0420/>`, PySCF has developed a
couple of methods to install the extension modules.

* Pypi command. For pypi version newer than 19.0, projects hosted on
  GitHub can be installed on the command line::

    $ pip install git+https://github.com/pyscf/semiemprical

  A particular release on github can be installed with the release URL
  you can look up on GitHub::

    $ pip install https://github.com/pyscf/semiemprical/archive/v0.1.0.tar.gz

* Pypi command for local paths. If you wish to load an extension
  module developed in a local directory, you can use the local install
  mode of pip. Use of a Python virtual environment is recommended to
  avoid polluting the system default Python runtime environment; for
  example::

    $ python -m venv /home/abc/pyscf-local-env
    $ source /home/abc/pyscf-local-env/bin/activate
    $ git clone https://github.com/pyscf/semiemprical /home/abc/semiemprical
    $ pip install -e /home/abc/semiemprical

* Environment variable `PYSCF_EXT_PATH`. You can place the location of
  each extension module (or a file that contains these locations) in
  this environment variable. The PySCF library will parse the paths
  defined in this environment variable, and load the relevent
  submodules. For example::

    $ git clone https://github.com/pyscf/semiempirical /home/abc/semiempirical
    $ git clone https://github.com/pyscf/doci /home/abc/doci
    $ git clone https://github.com/pyscf/dftd3 /home/abc/dftd3
    $ echo /home/abc/doci > /home/abc/.pyscf_ext_modules
    $ echo /home/abc/dftd3 >> /home/abc/.pyscf_ext_modules
    $ export PYSCF_EXT_PATH=/home/abc/semiempirical:/home/abc/.pyscf_ext_modules

  Using this definition of `PYSCF_EXT_PATH`, the three extension
  submodules (semiempirical, doci, dftd3) are loaded when PySCF is
  imported, and you don't have to use a Python virtual environment.

Once the extension modules have been correctly installed (with any of
the methods shown above), you can use them as regular submodules
developed inside the pyscf main project::

    >>> import pyscf
    >>> from pyscf.semiempirical import MINDO
    >>> mol = pyscf.M(atom='N 0 0 0; N 0 0 1')
    >>> MINDO(mol).run()


NAO
---

The :mod:`nao` module includes basic functions for numerical atomic
orbitals (NAO) and NAO-based TDDFT methods.  This module was
contributed by Marc Barbry and Peter Koval. More details of :mod:`nao`
can be found in
https://github.com/pyscf/nao/blob/master/README.md. This module can be
installed with::

    $ pip install https://github.com/pyscf/nao


DMRG solver
-----------

Density matrix renormalization group (DMRG) theory is a powerful
method for solving ab initio quantum chemistry problems. PySCF can be
used with two implementations of DMRG: Block
(https://sanshar.github.io/Block) and CheMPS2
(http://sebwouters.github.io/CheMPS2/index.html).  `Installing Block
<https://sanshar.github.io/Block/build.html>`_ requires a C++11
compiler.  If C++11 is not supported by your compiler, you can
register and download the precompiled Block binary from
https://sanshar.github.io/Block/build.html.  Before using Block or
CheMPS2, you need create a configuration file
``pyscf/dmrgscf/settings.py`` (as shown by settings.py.example) to
store the path where the DMRG solver was installed.


TBLIS
-----

`TBLIS <https://github.com/devinamatthews/tblis>`_ provides a native
algorithm for performing tensor contraction for arbitrarily
high-dimensional tensors. The native algorithm in TBLIS does not need
to transform tensors into matrices by permutations, then call BLAS for
the the matrix contraction, and back-permute the results. This means
that tensor transposes and data moves are largely avoided by TBLIS.
The interface to TBLIS offers an efficient implementation for
:func:`numpy.einsum` style tensor contraction.  The tlibs-einsum
plugin can be enabled with::

  $ pip install pyscf-tblis
