.. _solvent:

:mod:`solvent` --- Solvent methods
**********************************
.. module:: solvent
   :synopsis: Solvation via domain-decomposition COSMO for mean-field and correlated methods

The :mod:`solvent` module implements the domain-decomposition COSMO solvent model
for mean-field and correlated methods.

To run a mean-field calculation with the implicit solvent model, one can try the
following example::

    from pyscf import gto, scf, dft
    from pyscf import solvent
    mol = gto.M(atom='''
    C        0.000000    0.000000             -0.542500
    O        0.000000    0.000000              0.677500
    H        0.000000    0.9353074360871938   -1.082500
    H        0.000000   -0.9353074360871938   -1.082500
                ''',
                verbose = 4)
    mf = scf.RHF(mol)
    solvent.ddCOSMO(mf).run()


Examples
========

* :source:`examples/solvent/00-scf_with_ddcosmo.py`
* :source:`examples/solvent/01-casscf_with_ddcosmo.py`
* :source:`examples/solvent/02-casci_with_ddcosmo.py`
* :source:`examples/solvent/03-ccsd_with_ddcosmo.py`
* :source:`examples/solvent/04-pe_potfile_from_pyframe.py`
* :source:`examples/solvent/04-scf_with_pe.py`
* :source:`examples/solvent/20-state_specific_casci.py`
* :source:`examples/solvent/21-tddft_equilibrium_solvation.py`
* :source:`examples/solvent/21-tddft_geomopt.py`
* :source:`examples/solvent/22-with_qmmm.py`

Program reference
=================

.. automodule:: pyscf.solvent
 
domain decomposition PCM
------------------------

.. automodule:: pyscf.solvent.ddpcm
   :members:

domain decomposition COSMO
--------------------------

.. automodule:: pyscf.solvent.ddcosmo
   :members:

.. automodule:: pyscf.solvent.ddcosmo_grad
      :members:

