geomopt --- Geometry optimization
**********************************

The :mod:`geomopt` module implements geometry optimization via
interfaces to geomeTRIC (https://github.com/leeping/geomeTRIC) and PyBerny (https://github.com/azag0/pyberny).
For example, optimize the structure of N2 molecule using PyBerny::

    from pyscf import gto, scf
    from pyscf.geomopt.berny_solver import optimize
    mol = gto.M(atom='N 0 0 0; N 0 0 1.2', basis='ccpvdz')
    mf = scf.RHF(mol)
    mol_eq = optimize(mf)
    print(mol_eq.atom_coords())

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

