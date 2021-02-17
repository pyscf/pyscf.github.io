.. _developer_scf:

***********************************
Self-consistent field (SCF) methods
***********************************

*Modules*: :mod:`scf`, :mod:`pbc.scf`, :mod:`soscf`

Overview
========
The SCF iterative loop is the (kernel) function in hf.py
which takes a mean-field object. The current mean-field classes in PySCF
are ...

The kernel function carries out the following steps
- generating an initial guess (code) get_init_guess
- building the Fock matrix, which consists of
  - building the 1-e part get_hcore
  - building the 2-e part get_veff (depends on the DM, and internally
    calls get_jk, which calls get_j and get_k)
  - obtaining the 1-e overlap
  - assembling the Fock matrix from hcore and veff and extrapolating
    by DIIS (get_fock)
  - testing for convergence (get_grad)
  - computing the new eigenvectors and eigenvalues, and filling the
    new DM (eig, get_occ) and assembling the new DM (make_rdm1)
  - updating the energy (energy_tot, which calls energy_elec and energy_nuc)
  - writing results to the checkfile (dump_chk)

Internally, different methods reuse this kernel by overwriting the
mean-field methods.

- DFT SCF is implemented by
specializing mf.get_veff and mf.energy_elec for the various KS objects (see e.g. RKS.py).

- Density Fitting is implemented by
  
- Gamma point PBC SCF re-implements (list of functions)

- k-point PBC SCF re-implements


Incore implementation
=====================
Define simplest one to modify: here 2e integrals are incore.

Custom Hamiltonians
===================

The HF approximation SCF procedure can be used for arbitrary Hamiltonians.
This is useful when considering model Hamiltonians. The key functions to supply
are the initial guess (which can also be supplied by supplying an initial DM),
get_ovlp, get_hcore, and the two electron integrals.
Here show examples of using user defined Hamiltonians






    



