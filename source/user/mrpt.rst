.. _user_mrpt:

Multi-Reference Perturbation Theory (MRPT)
******************************************

*Modules*: :mod:`mrpt`

Introduction
============

The second-order perturbative energy correction can be added on top of the multireference wave function.
The MRPT module in PySCF supports the second-order N-electron valence state perturbation theory (NEVPT2) :cite:`Angeli2001NEVPT`  
using the strongly contracted (SC) internal contraction scheme, :cite:`Angeli2001SCNEVPT,Angeli2002` which is an intruder-state-free MRPT.
SC-NEVPT2 can be applied to the CASCI/CASSCF wave function optimized by FCI and DMRG solvers. :cite:`Guo2016`
Root_ID needs to be specified for the state-specific NEVPT2 calculation as ``mrpt.NEVPT(mc,root=Root_ID)``.
(by default Root_ID=0 for the lowest root)

A simple example (see :source:`examples/mrpt/03-dmrg_nevpt2.py`) of running
SC-NEVPT2 calculations with FCI and DMRG solvers is 

.. literalinclude:: ../../examples/mrpt/03-dmrg_nevpt2.py
which outputs

.. code::

  CI NEVPT = -0.0655571220984536 -0.0805814569225376
  MPS NEVPT correlation E = -0.0655570680554582 -0.0913913618723217

namely, the second-order correlation energies for the ground and the first-excited states 
with the FCI solver, and those with the DMRG solver.

Compressed Perturber Functions
==============================

The bottleneck in SC-NEVPT2 is the evaluation of the energies of the perturber functions,
where up to the 4-particle reduced density matrix (4-RDM) appears.
In DMRG-SC-NEVPT2, this evaluation is done with the compressed bond dimension (M') which is smaller than 
the bond dimension in the DMRG energy optimization as
``mrpt.NEVPT(mc,root=Root_ID).compress_approx(maxM=M')``

More information can be found in Reference :cite:`Guo2016` 

References
==========

.. bibliography:: ref_mrpt.bib
   :style: unsrt
