.. _user_mrpt:

Multi-reference perturbation theory (MRPT)
******************************************

*Modules*: :mod:`mrpt`

Introduction
============

A second-order perturbative energy correction can be added on top of a multireference wave function.
The MRPT module in PySCF supports the second-order N-electron valence state perturbation theory (NEVPT2) :cite:`Angeli2001NEVPT`  
using the strongly contracted (SC) internal contraction scheme, :cite:`Angeli2001SCNEVPT,Angeli2002` which is an intruder-state-free MRPT.
SC-NEVPT2 can be applied to CASCI/CASSCF wave functions produced by the FCI or DMRG solvers. :cite:`Guo2016`
The number of the CI root needs to be specified for state-specific NEVPT2 calculations with ``mrpt.NEVPT(mc,root=Root_ID)``.
By default, the NEVPT2 calculation is performed for the lowest root, `Root_ID=0`.

A simple example of SC-NEVPT2 calculations with the FCI and DMRG solvers is given in :source:`examples/mrpt/03-dmrg_nevpt2.py` 

.. literalinclude:: ../../examples/mrpt/03-dmrg_nevpt2.py

which outputs

.. code::

  FCI NEVPT correlation E = -0.0655576579894365 -0.0913916717329482 
  MPS NEVPT correlation E = -0.0655570680554582 -0.0913913618723217

namely, the second-order correlation energies for the ground and the first-excited states 
with the FCI solver, and those with the DMRG solver.

Compressed Perturber Functions
==============================

The bottleneck in SC-NEVPT2 is the evaluation of the energies of the perturber functions,
where up to the 4-particle reduced density matrix (4-RDM) appears.
In DMRG-SC-NEVPT2, this evaluation is done with the compressed bond dimension (M') which is smaller than 
the bond dimension in the DMRG energy optimization to avoid the bottleneck.
It can be specified by ``mrpt.NEVPT(mc,root=Root_ID).compress_approx(maxM=M')``.

More information can be found in Reference :cite:`Guo2016` 

References
==========

.. bibliography:: ref_mrpt.bib
   :style: unsrt
