.. _developer_scf:

***********************************
Self-consistent field (SCF) methods
***********************************

*Modules*: :mod:`scf`, :mod:`pbc.scf`, :mod:`soscf`

Overview
========
The SCF iterative loop is the :func:`pyscf.scf.hf.kernel` function
which takes a mean-field object. The current mean-field classes in PySCF
are ...

The :func:`kernel` function carries out the following steps:

- generating an initial guess -- :func:`get_init_guess`

- building the Fock matrix, which consists of

  - building the 1-e part -- :func:`get_hcore`

  - building the 2-e part -- :func:`get_veff` (depends on the DM, and internally
    calls :func:`get_jk`, which calls :func:`get_j` and :func:`get_k`)

  - obtaining the 1-e overlap -- :func:`get_ovlp`

  - assembling the Fock matrix from ``hcore`` and ``veff`` and extrapolating
    by DIIS -- :func:`get_fock`

  - testing for convergence -- :func:`get_grad`

  - computing the new eigenvectors and eigenvalues, and filling the
    new DM -- :func:`eig` and :func:`get_occ`, and assembling the new DM -- :func:`make_rdm1`

  - updating the energy -- :func:`energy_tot`, which calls :func:`energy_elec` and :func:`energy_nuc`

  - writing results to the checkfile -- :func:`dump_chk`

Internally, different methods reuse this kernel by overwriting the
mean-field methods.

- DFT SCF is implemented by
  specializing :meth:`get_veff` and :meth:`energy_elec` for the various KS objects 
  (see e.g. :func:`pyscf.dft.rks.get_veff` and :func:`pyscf.dft.rks.energy_elec`).

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
:meth:`get_ovlp`, :meth:`get_hcore`, and the two electron integrals 
(attribute :attr:`._eri` of the mean-field object).
Following shows an example of HF with a Hubbard model Hamiltonian::

    import numpy
    from pyscf import gto, scf, ao2mo, mcscf

    mol = gto.M()

    # incore_anyway=True ensures the customized Hamiltonian (the _eri attribute)
    # to be used.  Without this parameter, the MO integral transformation may
    # ignore the customized Hamiltonian if memory is not enough.
    mol.incore_anyway = True

    n = 12 # No. of sites
    mol.nelectron = n // 2 # half filling

    h1 = numpy.zeros((n,n))
    for i in range(n-1):
        h1[i,i+1] = h1[i+1,i] = -1.0
    h1[n-1,0] = h1[0,n-1] = -1.0  # PBC
    eri = numpy.zeros((n,n,n,n))
    for i in range(n):
        eri[i,i,i,i] = 4.0

    mf = scf.RHF(mol)
    mf.get_hcore = lambda *args: h1
    mf.get_ovlp = lambda *args: numpy.eye(n)
    # ao2mo.restore(8, eri, n) to get 8-fold permutation symmetry of the integrals
    # ._eri only supports the two-electron integrals in 4-fold or 8-fold symmetry.
    mf._eri = ao2mo.restore(8, eri, n)
    mf.init_guess = '1e'
    mf.kernel()

    # post-SCF calculation
    mycas = mcscf.CASSCF(mf, 4, 4)
    mycas.kernel()
