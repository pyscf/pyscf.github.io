.. _mrpt:

:mod:`mrpt` --- Multi-reference perturbation theory
***************************************************

.. automodule:: pyscf.mrpt

The :mod:`mrpt` module implements the N-electron valence
(multi-reference) perturbation theory.


DMRG-NEVPT2
===========
In the calculations of DMRG-CASSCF followed by DMRG-NEVPT2,
:attr:`canonicalization` and :attr:`natorb` of DMRG-CASSCF should be set.
These settings ensure that DMRG-CASSCF orbitals are eigenvectors of general Fock
matrix. Current NEVPT2 implementation is based on this assumption. (See also
discussions in github issue https://github.com/pyscf/pyscf/issues/698)


Program reference
=================

N-electron valance perturbation theory (NEVPT2)
-----------------------------------------------

.. automodule:: pyscf.mrpt.nevpt2
   :members:
