.. _pbc:

*******************************************
pbc --- Periodic boundary conditions
*******************************************
The :mod:`pbc` module provides electronic structure implementations with periodic boundary
conditions based on periodic Gaussian basis functions. The PBC implementation supports
both all-electron and pseudopotential descriptions.

In PySCF, the PBC implementation has closely related to the molecular implementation.
The module names, function names, and layout of the PBC code are the same as (or as close
as possible to) those of the molecular code.  The PBC code supports the use (and mixing)
of basis sets, pseudopotentials, and effective core potentials developed across the
materials science and quantum chemistry communities, offering great flexibility.  Moreover,
many post-mean-field methods defined in the molecular code can be seamlessly mixed with
PBC calculations performed at the Gamma point.  For example, one can perform a Gamma-point
Hartree-Fock calculation in a supercell, followed by a CCSD(T) calculation, which is
implemented in the molecular code.

In the PBC k-point module,
we implement minor changes to the Gamma-point data structures. The associated classes and methods
are prefixed by "K"; for example, the mean-field k-point restricted Hartree-Fock and Kohn-Sham modules are KRHF and KRKS. 
On top of these KSCF methods, one can find many correlation methods, such as k-point CCSD and k-point EOM-CCSD methods.
Other post-mean-field methods can be also be developed to explicitly enforce translational symmetry through k-point sampling.

When using results of this code for publication, please cite the following papers:

1) "Gaussian-Based Coupled-Cluster Theory for the Ground-State and Band Structure of Solids" J. McClain, Q. Sun, G. K.-L. Chan, and T. C. Berkelbach, J. Chem. Theory Comput. 13, 1209 (2017).

2) "Gaussian and plane-wave mixed density fitting for periodic systems" Q. Sun, T. C. Berkelbach, J. McClain, G. K.-L. Chan, J. Chem. Phys. 147, 164119 (2017).

The list of modules described in this chapter is:

.. toctree::

   pbc/gto.rst
   pbc/ao2mo.rst
   pbc/scf.rst
   pbc/dft.rst
   pbc/df.rst
   pbc/tdscf.rst
   pbc/mp.rst
   pbc/ci.rst
   pbc/cc.rst
   pbc/mpicc.rst
   pbc/tools.rst
   pbc/mix_mol.rst
   pbc/prop.rst
   pbc/gw.rst
