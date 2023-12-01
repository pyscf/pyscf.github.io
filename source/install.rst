.. _installing:

Install PySCF
*************

1) Install with `pip` (easiest method)
======================================
This is the recommended way to install PySCF for non-developers::

  $ pip install --prefer-binary pyscf

The pip package provides a precompiled PySCF code (python wheel) which
works on almost all Linux systems, and most of Mac OS X systems, and
the Ubuntu subsystems on Windows 10. If you already have installed
PySCF via pip, you can upgrade it to the new version with::

  $ pip install --upgrade pyscf

.. note::
   Since PySCF version 2.1, the Linux wheels require manylinux2010 (for x86_64) or manylinux2014 (for aarch64). So the pip version should >= 19.3 for installing on Linux.


2) Build from source with `pip`
===============================

If you're interested in a new feature, that's not included in the latest release or you simply 
want the latest and greatest PySCF you can build from source using pip.::

  $ pip install git+https://github.com/pyscf/pyscf

To install the features developed on a particular branch use::

  $ pip install git+https://github.com/pyscf/pyscf@<branch_name>

This install method compiles and links C extensions against the
libraries in your system. See :ref:`compile_c_extensions` for a full 
list of prerequisites. If you would like to tune the CMake compilation 
parameters, you can set them with the environment variable `CMAKE_CONFIGURE_ARGS`, 
for example:: 

  $ export CMAKE_CONFIGURE_ARGS="-DBUILD_MARCH_NATIVE=ON"

See :ref:`cmake_options` for more details about CMake configuration.

.. _compile_c_extensions:

3) Build from source
====================

Prerequisites for manual install are

.. note::

  * C compiler
  * C++ compiler (optional, but required for XCFun and some extensions)
  * CMake >= 3.10
  * Python >= 3.7
  * Numpy >= 1.13
  * Scipy >= 0.19
  * h5py >= 2.7

You can download the latest version of PySCF (or the development
branch) from github::

  $ git clone https://github.com/pyscf/pyscf.git
  $ cd pyscf

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

See :ref:`cmake_options` for details about CMake configuration.



4) Installation with conda
==========================

If you have a `Conda <https://conda.io/docs/>`_ (or `Anaconda
<https://www.continuum.io/downloads#linux>`_) environment, PySCF
package can be installed from the Conda cloud
(for Linux and Mac OS X systems)::

  $ conda install -c pyscf pyscf

Extension modules are not available on the Conda cloud. They should be
installed either with pip, or through the environment variable
`PYSCF_EXT_PATH` (see the section :ref:`installing_extproj`).


5) Installation on Fedora
==========================

If you are running Fedora Linux, you can install PySCF as a
distribution package::

  # dnf install python3-pyscf

If you are running on an X86-64 platform, dnf should automatically
install the optimized integral library, qcint, instead of the
cross-platform libcint library.

Extension modules are not available in the Fedora package.



6) PySCF docker image
=====================

The following command starts a container with the jupyter notebook
server that listens for HTTP connections on port 8888::

  $ docker run -it -p 8888:8888 pyscf/pyscf:latest

Now, you can visit ``https://localhost:8888`` with your browser to use
PySCF in the notebook.

Another way to use PySCF in a docker container is to start an Ipython
shell::

  $ docker run -it pyscf/pyscf:latest start.sh ipython


Advanced build options
**********************

.. _cmake_options:

CMake options
=============

A full build of PySCF may take a long time to finish.
`XCFun` may fail to build if a proper C++ compiler is not available, such as on certain old operating systems.
The CMake options listed below can be used to speed up compilation or omit extensions that fail to compile.
Note:  If both `-DENABLE_LIBXC=OFF` and `-DENABLE_XCFUN=OFF` are set, importing the dft module will lead to an `ImportError`.

==================== ======= =================================================================
Flags                Default Comments
==================== ======= =================================================================
`ENABLE_LIBXC`       ON      Whether to use `LibXC` library in PySCF. If `-DENABLE_LIBXC=OFF`
                             is appended to cmake command, `LibXC` will not be compiled.
`ENABLE_XCFUN`       ON      Whether to use `XCFun` library in PySCF. If `-DENABLE_XCFUN=OFF`
                             is appended to cmake command, `XCFun` will not be compiled.
`BUILD_LIBXC`        ON      Set it to `OFF` to skip compiling `Libxc`. The dft module
                             still calls `LibXC` library by default. The dft module will be
                             linked against the `LibXC` library from an earlier build.
`BUILD_XCFUN`        ON      Set it to `OFF` to skip compiling `XCFun`. The dft module
                             will be linked against the `XCFun` library from an earlier build.
`BUILD_LIBCINT`      ON      Set it to `OFF` to skip compiling `libcint`. The integral
                             library from an earlier build will be used.
`WITH_F12`           ON      Whether to compile the F12 integrals.
`DISABLE_DFT`        OFF     Set this flag to skip the entire dft module. Neither `LibXC`
                             nor `XCFun` will be compiled.
`BUILD_MARCH_NATIVE` OFF     Whether to let the compiler optimize the code against CPU
                             architecture
==================== ======= =================================================================

CMake config file
-----------------

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
directory ``pyscf/lib/cmake_user_inc_examples``.


Environment variables and global configures
===========================================

======================= =========================================================
Env variable            Comments
======================= =========================================================
`PYSCF_MAX_MEMORY`      Maximum memory to use in MB
`PYSCF_TMPDIR`          Directory for temporary files
`PYSCF_CONFIG_FILE`     File where various PySCF default settings are stored
`PYSCF_EXT_PATH`        Path for finding external extensions
======================= =========================================================

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
download the libcint, Libxc, and XCFun libraries::

    $ git clone https://github.com/sunqm/libcint.git
    $ tar czf libcint.tar.gz libcint

    $ wget https://gitlab.com/libxc/libxc/-/archive/6.0.0/libxc-6.0.0.tar.gz

    $ wget -O xcfun.tar.gz https://github.com/fishjojo/xcfun/archive/refs/tags/cmake-3.5.tar.gz

Assuming ``/opt`` is the place where these libraries will be installed, these
packages should be compiled with the flags::

    $ tar xvzf libcint.tar.gz
    $ cd libcint
    $ mkdir build && cd build
    $ cmake -DWITH_F12=1 -DWITH_RANGE_COULOMB=1 -DWITH_COULOMB_ERF=1 \
        -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
    $ make && make install

    $ tar xvzf libxc-6.0.0.tar.gz
    $ cd libxc-6.0.0
    $ mkdir build && cd build
    $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=1 \
        -DENABLE_FORTRAN=0 -DDISABLE_KXC=0 -DDISABLE_LXC=1 \
        -DCMAKE_INSTALL_PREFIX:PATH=/opt -DCMAKE_INSTALL_LIBDIR:PATH=lib ..
    $ make && make install

    $ tar xvzf xcfun.tar.gz
    $ cd xcfun-cmake-3.5
    $ mkdir build && cd build
    $ cmake -DCMAKE_BUILD_TYPE=RELEASE -DBUILD_SHARED_LIBS=1 -DXCFUN_MAX_ORDER=3 -DXCFUN_ENABLE_TESTS=0 \
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


.. _installing_extproj:

Install PySCF extensions
************************

Extension modules
=================

As of PySCF-2.0, some modules have been moved from the main code trunk
to extension projects hosted at https://github.com/pyscf.

=================== =========================================================
Project             URL
=================== =========================================================
cas_ac0             https://github.com/CQCL/pyscf-ac0
cornell-shci        https://github.com/pyscf/cornell-shci
dftd3               https://github.com/pyscf/dftd3
dmrgscf             https://github.com/pyscf/dmrgscf
doci                https://github.com/pyscf/doci
fciqmc              https://github.com/pyscf/fciqmc
icmpspt             https://github.com/pyscf/icmpspt
mbd                 https://github.com/pyscf/mbd
naive-hci           https://github.com/pyscf/naive-hci
nao                 https://github.com/pyscf/nao
qsdopt              https://github.com/pyscf/qsdopt
rt                  https://github.com/pyscf/rt
semiempirical       https://github.com/pyscf/semiempirical
shciscf             https://github.com/pyscf/shciscf
zquatev             https://github.com/sunqm/zquatev
tblis               https://github.com/pyscf/pyscf-tblis
=================== =========================================================

Install extensions
==================

Since PySCF version 2.0, some modules are developed independently; see
:ref:`installing_extproj`. Individual extension modules (for example
the geometry optimization module) can be installed using pip's extra
dependency mechanism::

  $ pip install pyscf[geomopt]

All extension modules can be installed with::

  $ pip install pyscf[all]

The extension modules can be found in `https://github.com/pyscf` (see
also :ref:`installing_extproj`).

Install extensions (advanced)
=============================

Based on the technique of namespace packages specified in `PEP 420
<https://www.python.org/dev/peps/pep-0420/>`_, PySCF has developed a
couple of methods to install the extension modules.

* Pip command. For pip version newer than 19.0, projects hosted on
  GitHub can be installed on the command line::

    $ pip install git+https://github.com/pyscf/semiempirical

  A particular release on github can be installed with the release URL
  you can look up on GitHub::

    $ pip install https://github.com/pyscf/semiempirical/archive/v0.1.0.tar.gz

* Pip command for local paths. If you wish to load an extension
  module developed in a local directory, you can use the local install
  mode of pip. Use of a Python virtual environment is recommended to
  avoid polluting the system default Python runtime environment; for
  example::

    $ python -m venv /home/abc/pyscf-local-env
    $ source /home/abc/pyscf-local-env/bin/activate
    $ git clone https://github.com/pyscf/semiempirical /home/abc/semiempirical
    $ pip install -e /home/abc/semiempirical

* Environment variable `PYSCF_EXT_PATH`. You can place the location of
  each extension module (or a file that contains these locations) in
  this environment variable. The PySCF library will parse the paths
  defined in this environment variable, and load the relevant
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
    >>> from pyscf.semiempirical import MINDO3
    >>> mol = pyscf.M(atom='N 0 0 0; N 0 0 1')
    >>> MINDO(mol).run()

Common examples
===============
... NAO
... ---
... The :mod:`nao` module includes basic functions for numerical atomic
orbitals (NAO) and NAO-based TDDFT methods.  This module was
contributed by Marc Barbry and Peter Koval. More details of :mod:`nao`
can be found in
https://github.com/pyscf/nao/blob/master/README.md. This module can be
installed with::
...    $ pip install https://github.com/pyscf/nao


DMRG solvers
------------

Density matrix renormalization group (DMRG) theory is a powerful
method for solving ab initio quantum chemistry problems. PySCF can be
used with three implementations of DMRG: Block
(https://sanshar.github.io/Block), block2
(https://block2.readthedocs.io/en/latest), and CheMPS2
(http://sebwouters.github.io/CheMPS2/index.html).

`Installing Block <https://sanshar.github.io/Block/build.html>`_ requires a C++11
compiler.  If C++11 is not supported by your compiler, you can
download the precompiled Block binary from https://sanshar.github.io/Block/build.html.

``block2`` can be easily installed via ``pip install block2`` or ``pip install block2-mpi``,
or `building from source <https://block2.readthedocs.io/en/latest/user/installation.html>`_.

Before using Block or CheMPS2, you need create a configuration file
``pyscf/dmrgscf/settings.py`` (as shown by settings.py.example) to
store the path where the DMRG solver was installed.


TBLIS
-----

`TBLIS <https://github.com/devinamatthews/tblis>`_ provides a native
algorithm for performing tensor contraction for arbitrarily
high-dimensional tensors. The native algorithm in TBLIS does not need
to transform tensors into matrices by permutations, then call BLAS for
the the matrix contraction, and back-permute the results. This means
that tensor transposes and data moves are largely avoided by TBLIS. This
leads to speedups in many correlated quantum chemistry methods in PySCF, such as
the coupled cluster methods.
The interface to TBLIS offers an efficient implementation for
:func:`numpy.einsum` style tensor contraction.  The tblis-einsum
plugin can be enabled with::

  $ pip install pyscf-tblis

Troubleshooting
***************

`error: command 'cmake' failed`
===============================

In some cases, users who install PySCF with `pip install pyscf` may see an error like the following::

  Building wheels for collected packages: pyscf
    Building wheel for pyscf (setup.py) ... error
    error: subprocess-exited-with-error
    × python setup.py bdist_wheel did not run successfully.
    │ exit code: 1
    ╰─> [7 lines of output]
        scipy>1.1.0 may crash when calling scipy.linalg.eigh. (Issues https://github.com/scipy/scipy/issues/15362 https://github.com/scipy/scipy/issues/16151)
        running bdist_wheel
        running build
        running build_ext
        Configuring extensions
        cmake -S/Users/<user>/personal/codes/chemistry/pyscf/pyscf/lib -Bbuild/temp.macosx-12-x86_64-cpython-310
        error: command 'cmake' failed: No such file or directory
        [end of output]

Here, ``pip`` chose not to install a binary wheel and is trying to build from source. 
If that's not your intention, you should install with the command ``pip install --prefer-binary pyscf``.
On the other hand, if you are intentionally trying to build from source, you're missing the required `cmake` program.
See the docs for building from source above and issue `1684 <https://github.com/pyscf/pyscf/issues/1684>`_ for more details.

MacOS: `Library not loaded`
===========================

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
