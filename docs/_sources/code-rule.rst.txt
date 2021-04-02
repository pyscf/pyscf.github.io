.. _code_stand:

Code standard
*************

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


Name convention
===============

* The prefix or suffix underscore in the function names have special meanings

  - functions with prefix-underscore like ``_fn`` are private functions.
    They are typically not documented, and not recommended to use.

  - functions with suffix-underscore like ``fn_`` means that they have side
    effects.  The side effects include the change of the input arguments,
    the runtime modification of the class definitions (attributes or
    members), or module definitions (global variables or functions) etc.

  - regular (pure) functions do not have underscore as the prefix or suffix.

API convention
==============

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

  - In __init__ function, initialize/define the problem size.  The
    problem size parameters (like num_orbitals etc) can be considered as
    environments.  They should be immutable.

  - Kernel functions:
    Classes for QC models should provide a method :func:`kernel` as the entrance/main function.
    The :func:`kernel` function then call other code to finish the calculation.
    Although not required, it is recommended to let the kernel function return certain key results.
    If your class is inherited from the :class:`pyscf.lib.StreamObject`,
    the class has a method :func:`run` which will call the :func:`kernel`
    function and return the object itself. One can simply call the
    :func:`kernel` method or :func:`run` method to start the flow of a QC method.

* Function arguments

  - The first argument is a handler.  The handler is one of :class:`gto.Mole`
    object, a mean-field object, or a post-Hartree-Fock object.

.. include:: design.rst
