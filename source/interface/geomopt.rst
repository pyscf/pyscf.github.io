.. _geomopt:

:mod:`geomopt` --- Geometry optimization
****************************************

The :mod:`geomopt` module implements geometry optimization via
interfaces to `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ 
and `PyBerny <https://github.com/azag0/pyberny>`_.
The following example shows how to optimize the structure of N\ :sub:`2`\  molecule with PyBerny::

    from pyscf import gto, scf
    from pyscf.geomopt.berny_solver import optimize
    mol = gto.M(atom='N 0 0 0; N 0 0 1.2', basis='ccpvdz')
    mf = scf.RHF(mol)
    mol_eq = optimize(mf)
    print(mol_eq.atom_coords())

Examples
========

:download:`examples/geomopt/01-geomeTRIC.py </../examples/geomopt/01-geomeTRIC.py>`
:download:`examples/geomopt/01-pyberny.py </../examples/geomopt/01-pyberny.py>`
:download:`examples/geomopt/02-as_pyscf_method.py </../examples/geomopt/02-as_pyscf_method.py>`
:download:`examples/geomopt/10-with_qmmm.py </../examples/geomopt/10-with_qmmm.py>`
:download:`examples/geomopt/12-excited_states.py </../examples/geomopt/12-excited_states.py>`
:download:`examples/geomopt/12-mcscf_excited_states.py </../examples/geomopt/12-mcscf_excited_states.py>`
:download:`examples/geomopt/13-ccsd_t.py </../examples/geomopt/13-ccsd_t.py>`
:download:`examples/geomopt/14-with_solvent.py </../examples/geomopt/14-with_solvent.py>`
:download:`examples/geomopt/20-callback.py </../examples/geomopt/20-callback.py>`


Program reference
=================

.. automodule:: pyscf.geomopt.berny_solver
   :members:

.. automodule:: pyscf.geomopt.geometric_solver
   :members:

.. automodule:: pyscf.geomopt.addons
   :members:
