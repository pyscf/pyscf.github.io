.. _developer_scf:

***********************************
Self-consistent field (SCF) methods
***********************************

*Modules*: :mod:`scf`, :mod:`pbc.scf`, :mod:`soscf`

Overview
========
The class :class:`pyscf.scf.hf.SCF` is the base class for all mean-field methods,
which are implemented as derived classes. These include

=================================== ======
:class:`pyscf.scf.hf.RHF`           RHF
:class:`pyscf.scf.uhf.UHF`          UHF
:class:`pyscf.scf.rohf.ROHF`        ROHF
:class:`pyscf.scf.ghf.GHF`          GHF
:class:`pyscf.scf.dhf.DHF`          Dirac-HF
:class:`pyscf.scf.dhf.RDHF`         restricted Dirac-HF
:class:`pyscf.dft.rks.RKS`          RKS-DFT
:class:`pyscf.dft.uks.UKS`          UKS-DFT
:class:`pyscf.dft.roks.ROKS`        ROKS-DFT
:class:`pyscf.dft.gks.GKS`          GKS-DFT
:class:`pyscf.dft.dks.DKS`          Dirac-KS-DFT
:class:`pyscf.dft.dks.RDKS`         restricted Dirac-KS-DFT
:class:`pyscf.pbc.scf.hf.SCF`       base class for SCF methods with PBC at :math:`\Gamma` point
:class:`pyscf.pbc.scf.hf.RHF`       :math:`\Gamma`-point RHF
:class:`pyscf.pbc.scf.uhf.UHF`      :math:`\Gamma`-point UHF
:class:`pyscf.pbc.scf.rohf.ROHF`    :math:`\Gamma`-point ROHF
:class:`pyscf.pbc.scf.ghf.GHF`      :math:`\Gamma`-point GHF
:class:`pyscf.pbc.dft.rks.RKS`      :math:`\Gamma`-point RKS-DFT
:class:`pyscf.pbc.dft.uks.UKS`      :math:`\Gamma`-point UKS-DFT
:class:`pyscf.pbc.dft.roks.ROKS`    :math:`\Gamma`-point ROKS-DFT
:class:`pyscf.pbc.scf.khf.KSCF`     base class for SCF methods with PBC and with k-point sampling 
:class:`pyscf.pbc.scf.khf.KRHF`     KRHF
:class:`pyscf.pbc.scf.kuhf.KUHF`    KUHF
:class:`pyscf.pbc.scf.krohf.KROHF`  KROHF
:class:`pyscf.pbc.scf.kghf.KGHF`    KGHF
:class:`pyscf.pbc.dft.krks.KRKS`    KRKS-DFT
:class:`pyscf.pbc.dft.kuks.KUKS`    KUKS-DFT
:class:`pyscf.pbc.dft.kroks.KROKS`  KROKS-DFT
=================================== ======

The key attributes of the SCF class include

==================  ====================
:attr:`init_guess`  initial guess method
:attr:`DIIS`        DIIS method, can be :class:`pyscf.scf.diis.DIIS`, :class:`pyscf.scf.diis.EDIIS`, or :class:`pyscf.scf.diis.ADIIS`
:attr:`mo_coeff`    saved MO coefficients
:attr:`mo_energy`   saved MO energies
:attr:`mo_occ`      saved MO occupations
==================  ====================

The SCF iterative loop is the :func:`pyscf.scf.hf.kernel` function
which takes a mean-field object.
The :func:`kernel` function carries out the following steps:

- generating an initial guess -- :func:`get_init_guess`

- building the Fock matrix, which consists of

  - building the 1-e part -- :func:`get_hcore`

  - building the 2-e part -- :func:`get_veff` (takes the DM as input, and internally
    calls :func:`get_jk`, which calls :func:`get_j` and :func:`get_k`)

  - assembling the Fock matrix from :math:`h_{\rm core}` and :math:`v_{\rm eff}` and extrapolating
    by DIIS -- :func:`get_fock`

- obtaining the 1-e overlap -- :func:`get_ovlp`

- computing the new eigenvectors and eigenvalues -- :func:`eig`, 
  and filling and assembling the new DM -- :func:`get_occ` and :func:`make_rdm1`

- testing for convergence -- :func:`get_grad`

- updating the energy -- :func:`energy_tot`, which calls :func:`energy_elec` and :func:`energy_nuc`

- optionally writing results to the checkfile -- :func:`dump_chk`

Internally, different methods reuse this kernel by overwriting the
mean-field methods.

- DFT SCF is implemented by
  specializing :meth:`get_veff` and :meth:`energy_elec` for the various KS objects 
  (see e.g., :func:`pyscf.dft.rks.get_veff` and :func:`pyscf.dft.rks.energy_elec`).

- Density Fitting (DF) is applied by overwriting the :meth:`get_jk` method 
  with the DF three center 2-e integrals. (see e.g., :func:`pyscf.df.df_jk.density_fit`)
  
- :math:`\Gamma`-point PBC SCF re-implements

  - :meth:`get_ovlp` and :meth:`get_hcore` with crystal orbital 1-e integrals and optionally PBC pseudopotentials

  - :meth:`get_jk` with PBC DF integrals (e.g., :func:`pyscf.pbc.df.df_jk.get_jk`)

  - :meth:`energy_nuc` with Ewald summation (:func:`pyscf.pbc.gto.ewald`)

- Based on the :math:`\Gamma`-point PBC SCF, 
  k-point PBC SCF computes :meth:`make_rdm1`, :meth:`get_ovlp`, :meth:`get_hcore`, :meth:`get_jk`, :meth:`eig` at each k-point,
  and re-implements :meth:`energy_elec` with a k-point summation.


Integrals
=========
The molecular HF operator relies on the following underlying AO integral functions:

- :meth:`get_hcore` calls 1-e integral functions 
  :meth:`mol.intor('int1e_kin')` and :meth:`mol.intor('int1e_nuc')` (or :meth:`mol.intor('ECPscalar')` if using ECP)
  to compute the kinetic energy and nuclear attraction matrices.

- :meth:`get_jk` computes Coulomb and exchange matrices for the input DM.
  If enough memory is provided, the two-electron repulsion integral is computed (:meth:`mol.intor('int2e')`)
  and saved as the attribute :attr:`_eri`, 
  and is then contracted with the DM (:func:`pyscf.scf.hf.dot_eri_dm`) (incore algorithm).
  Otherwise, a "direct" algorithm, where the AO integrals are computed on the fly,
  is used (:func:`pyscf.scf._vhf.direct`). The incore path can be forced by 
  setting :attr:`.incore_anyway` of the :class:`Mole` object to ``True``. 

Custom Hamiltonians
===================
The HF approximation SCF procedure can be used for arbitrary Hamiltonians.
This is useful when considering model Hamiltonians. The key functions to supply
are the initial guess (which can also be generated by supplying an initial DM),
:meth:`get_ovlp`, :meth:`get_hcore`, and the two electron integrals 
(attribute :attr:`._eri` of the mean-field object).
The following shows an example of HF with a Hubbard model Hamiltonian::

    import numpy
    from pyscf import gto, scf, ao2mo, mcscf

    mol = gto.M()

    # incore_anyway=True ensures the customized Hamiltonian (the _eri attribute)
    # is used.  Without this parameter,  MO integral transformation used in
    # subsequent post-HF calculations may 
    # ignore the customized Hamiltonian if there is not enough memory.
    mol.incore_anyway = True

    n = 10
    mol.nelectron = n

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
    # ._eri only supports two-electron integrals with 4-fold or 8-fold symmetry.
    mf._eri = ao2mo.restore(8, eri, n)
    mf.init_guess = '1e'
    mf.kernel()

    # post-SCF calculation
    mycas = mcscf.CASSCF(mf, 4, 4)
    mycas.kernel()
