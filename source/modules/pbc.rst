.. _pbc:

:mod:`pbc` --- Periodic boundary conditions
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
   :maxdepth: 1

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
   pbc/x2c.rst

Examples
========

* :source:`examples/pbc/00-input_cell.py`
* :source:`examples/pbc/04-input_basis.py`
* :source:`examples/pbc/05-input_pp.py`
* :source:`examples/pbc/06-input_ke_cutoff.py`
* :source:`examples/pbc/09-band_ase.py`
* :source:`examples/pbc/09-init_from_ase.py`
* :source:`examples/pbc/09-talk_to_ase.py`
* :source:`examples/pbc/10-gamma_point_scf.py`
* :source:`examples/pbc/11-gamma_point_all_electron_scf.py`
* :source:`examples/pbc/12-gamma_point_post_hf.py`
* :source:`examples/pbc/20-k_points_scf.py`
* :source:`examples/pbc/21-k_points_all_electron_scf.py`
* :source:`examples/pbc/22-k_points_ccsd.py`
* :source:`examples/pbc/22-k_points_mp2.py`
* :source:`examples/pbc/22-k_points_tddft.py`
* :source:`examples/pbc/23-smearing.py`
* :source:`examples/pbc/24-k_points_vs_gamma.py`
* :source:`examples/pbc/25-k_points_mpi_ccsd.py`
* :source:`examples/pbc/26-linear_dep.py`
* :source:`examples/pbc/27-multigrid.py`
* :source:`examples/pbc/28-charged_system.py`
* :source:`examples/pbc/29-eom_ccsd_Ta.py`
* :source:`examples/pbc/30-ao_integrals.py`
* :source:`examples/pbc/30-ao_value_on_grid.py`
* :source:`examples/pbc/30-mo_integrals.py`
* :source:`examples/pbc/30-overlap_periodic_cell.py`
* :source:`examples/pbc/31-low_dimensional_pbc.py`
* :source:`examples/pbc/31-pbc_0D_as_mol.py`
* :source:`examples/pbc/32-graphene.py`
* :source:`examples/pbc/33-soc_integrals.py`
* :source:`examples/pbc/34-pm_localization.py`
* :source:`examples/pbc/35-gaussian_density_fit.py`
* :source:`examples/pbc/36-ccsd_level_shift.py`
* :source:`examples/pbc/40-customize_kpts_hamiltonian.py`
* :source:`examples/pbc/41-pbc1d_real_space_sum.py`
* :source:`examples/pbc/42-k2gamma.py`
* :source:`examples/pbc/44-k_point_eom_ccsd_ghf.py`
* :source:`examples/pbc/46-k_point_eom_ccsd_rhf.py`

