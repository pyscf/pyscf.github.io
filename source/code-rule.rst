.. _code_stand:

Code standard
*************

General considerations
======================

PySCF program code base is designed to provide an environment convenient for
method developing, quick testing and calculations for systems of moderate size. We emphasize first the
simplicity, next the generality, finally the efficiency. We favor the
implementations which have clear structure with optimization at Python level.
When Python performance becomes a major bottleneck, C code can be implemented to
improve efficiency. The following is a set of
guidelines for contributing to package. These are just guidelines, not rules.
Feel free to propose changes.

* 90/10 functional/OOP, unless performance critical, functions are pure.

* 90/10 Python/C, only computational hot spots were written in C.

* To extend python function with C:

  - Except complex numbers and variable length array, following C89 (gnu89) standard for C code.
    http://flash-gordon.me.uk/ansi.c.txt

  - Following C89 (gnu89) standard for C code;

  - Using ctypes to call C functions

* Conservative on advanced language feature.

* Minimal dependence principle

  - Minimal requirements on 3rd party program or libraries.

  - Loose-coupling between modules so that the failure of one module can
    have minimal effects on other modules.

  - For 3rd party Python library, implementing either a back up plan or
    an error/exception handler to avoid breaking the import chain

* External C/Fortran libraries. These are libraries to be compiled and linked in
  the PySCF C extensions. Their compiling/linking flags are resolved in the cmake config system.
  - BLAS, FFTW: Yes.
  - LAPACK: Yes but not recommended.  These libraries can be used in the
    PySCF C level library. But we recommend to restructure your code and move
    the relevant linear algebra and sparse matrix operations to Python code.
  - MPI and other parallel libraries: No.
    The MPI communications should be implemented at python level through MPI4py
    library.

* Code format.
  It is recommended to write code compliance with the [PEP8](https://www.python.org/dev/peps/pep-0008/) style.


Name conventions
================

* The prefix or suffix underscore in the function names have special meanings

  - functions with prefix-underscore like ``_fn`` are private functions.
    They are typically not documented, and not recommended to use.

  - functions with suffix-underscore like ``fn_`` means that they have side
    effects.  The side effects include the change of the input arguments,
    the runtime modification of the class definitions (attributes or
    members), or module definitions (global variables or functions) etc.

  - regular (pure) functions do not have underscore as the prefix or suffix.

API conventions
===============

* :class:`gto.Mole` (or :class:`gto.Cell` for PBC calculations) holds all global
  parameters, like the log level, the max memory usage etc.  They are used as
  the default values for all other classes.

* Class for quantum chemistry models or algorithms

  - Most QC method classes (like HF, CASSCF, FCI, ...) have three attributes
    ``verbose``, ``stdout`` and ``max_memory`` which are copied directly from
    :class:`gto.Mole` (or :class:`gto.Cell`.  Overwriting these attributes only affects the behavior of the
    local instance for that method class.  In the following example,
    ``mf.verbose`` mutes all messages produced by :class:`RHF`
    method, and the output of :class:`MP2` is written in the log file
    ``example.log``::

    >>> from pyscf import gto, scf, mp
    >>> mol = gto.M(atom='H 0 0 0; H 0 0 1', verbose=5)
    >>> mf = scf.RHF(mol)
    >>> mf.verbose = 0
    >>> mf.kernel()
    >>> mp2 = mp.MP2(mf)
    >>> mp2.stdout = open('example.log', 'w')
    >>> mp2.kernel()

  - Method class are only to hold the options or environments (like
    convergence threshold, max iterations, ...) to control the
    behavior/convergence of the method. Intermediate status at runtime are
    **not** supposed to be saved in the method class (in contrast to the object
    oriented paradigm).  However, the final results or outputs can be kept in
    the method object so that they can be easily accessed in the subsequent steps.
    We need to assume the attributes for results
    will be used as default inputs or environments for other objects in the rest
    parts of the program.
    The results attributes should be immutable, once they were generated
    and stored (after calling the :func:`kernel()` method) in a particular object.

  - In __init__ function, initialize/define the problem size. The
    problem size parameters (like num_orbitals etc) can be considered as
    environments. They should be immutable.

  - Kernel functions:
    Classes for QC models should provide a method :func:`kernel` as the entrance/main function.
    The :func:`kernel` function then call other code to finish the calculation.
    Although not required, it is recommended to let the kernel function return certain key results.
    If your class is inherited from the :class:`pyscf.lib.StreamObject`,
    the class has a method :func:`run` which will call the :func:`kernel`
    function and return the object itself. One can simply call the
    :func:`kernel` method or :func:`run` method to start the flow of a QC method.

* Function arguments

  - The first argument is a handler. The handler is one of :class:`gto.Mole`
    object, a mean-field object, or a post-Hartree-Fock object.

* Return value.
  Create returns for all functions whenever possible. For methods
  defined in class, return self instead of None if the method does not have
  particular return values.


Unit Tests and Example Scripts
==============================

* Examples to run modules should be placed in the appropriate directory inside
  the /examples directory.  While the examples should be easy enough to run on a
  modest personal computer; however, should not be trivial and instead showcase
  the functionality of the module.  The format for naming examples is::

    /examples/name_of_module/XX-function_name.py

  where XX is a two-digit numeric string.

* Test cases are placed in the /test/name_of_module directory and performed with
  nosetest (https://nose.readthedocs.io/en/latest/). These tests are to ensure
  the robustness of both simple functions and more complex drivers between
  version changes.

.. include:: design.rst
