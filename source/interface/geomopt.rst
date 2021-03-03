.. _geomopt:

:mod:`geomopt` --- Geometry optimization
****************************************

The :mod:`geomopt` module implements geometry optimization via
interfaces to `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ 
and `PyBerny <https://github.com/jhrmnn/pyberny>`_.
The following example shows how to optimize the structure of N\ :sub:`2`\  molecule with PyBerny::

    from pyscf import gto, scf
    from pyscf.geomopt.berny_solver import optimize
    mol = gto.M(atom='N 0 0 0; N 0 0 1.2', basis='ccpvdz')
    mf = scf.RHF(mol)
    mol_eq = optimize(mf)
    print(mol_eq.atom_coords())

Examples
========

* :source:`examples/geomopt/01-geomeTRIC.py`
* :source:`examples/geomopt/01-pyberny.py`
* :source:`examples/geomopt/02-as_pyscf_method.py`
* :source:`examples/geomopt/10-with_qmmm.py`
* :source:`examples/geomopt/11-with_ghost_atom.py`
* :source:`examples/geomopt/12-excited_states.py`
* :source:`examples/geomopt/12-mcscf_excited_states.py`
* :source:`examples/geomopt/13-ccsd_t.py`
* :source:`examples/geomopt/14-with_solvent.py`
* :source:`examples/geomopt/15-tddft_with_solvent.py`
* :source:`examples/geomopt/20-callback.py`

Program reference
=================

.. automodule:: pyscf.geomopt.berny_solver
   :members:

.. automodule:: pyscf.geomopt.geometric_solver
   :members:

.. automodule:: pyscf.geomopt.addons
   :members:
