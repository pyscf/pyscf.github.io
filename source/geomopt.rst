geomopt --- Geometry optimization
**********************************

The :mod:`geomopt` module implements geometry optimization via
interfaces to geomeTRIC (https://github.com/leeping/geomeTRIC) and PyBerny (https://github.com/azag0/pyberny).

Examples
========

Relevant examples
:file:`examples/geomopt/01-geomeTRIC.py`
:file:`examples/geomopt/01-pyberny.py`
:file:`examples/geomopt/02-as_pyscf_method.py`
:file:`examples/geomopt/10-with_qmmm.py`
:file:`examples/geomopt/11-with_ghost_atom.py`
:file:`examples/geomopt/12-excited_states.py`
:file:`examples/geomopt/12-mcscf_excited_states.py`
:file:`examples/geomopt/13-ccsd_t.py`


Program reference
=================

.. automodule:: pyscf.geomopt
 
pyberny
------

.. automodule:: pyscf.geomopt.berny_solver
   :members:

geomeTRIC
-------

.. automodule:: pyscf.geomopt.geometric_solver
   :members:


addons
------

.. automodule:: pyscf.geomopt.addons
   :members:

