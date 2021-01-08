.. _tdscf:

:mod:`tdscf` --- TDHF and TDDFT
*******************************

The :mod:`tdscf` module implements the time-dependent Hartree-Fock and 
time-dependent density functional theory.

An example to run a TDDFT calculation::
    
    from pyscf import gto, scf, dft, tddft
    mol = gto.Mole()
    mol.build(
        atom = 'H 0 0 0; F 0 0 1.1',
        basis = '631g',
        symmetry = True,
    )
    mf = dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.kernel()
    mytd = tddft.TDDFT(mf)
    mytd.kernel()
    mytd.analyze()

One can perform NTO analysis for TDDFT as::
    
    weights_1, nto_1 = mytd.get_nto(state=1, verbose=4)
    weights_2, nto_2 = mytd.get_nto(state=2, verbose=4)
    weights_3, nto_3 = mytd.get_nto(state=3, verbose=4)


Examples
========

* :source:`examples/tddft/00-simple_tddft.py`
* :source:`examples/tddft/01-nto_analysis.py`
* :source:`examples/tddft/02-tddft_for_camb3lyp.py`
* :source:`examples/tddft/21-matrix_A_B.py`
* :source:`examples/tddft/22-density.py`
* :source:`examples/tddft/30-change_xc_grids.py`
* :source:`examples/tddft/31-energy_transfer_coupling_matrix.py`

Program reference
=================

.. automodule:: pyscf.tdscf
   :members:

.. automodule:: pyscf.tdscf.rhf
   :members:

.. automodule:: pyscf.tdscf.uhf
   :members:

.. automodule:: pyscf.tdscf.rks
   :members:

.. automodule:: pyscf.tdscf.uks
   :members:

.. automodule:: pyscf.tdscf.common_slow
   :members:

.. automodule:: pyscf.tdscf.rhf_slow
   :members:

.. automodule:: pyscf.tdscf.proxy
   :members:

