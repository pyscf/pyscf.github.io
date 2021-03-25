
Quickstart
**********

The present tutorial is meant to provide a brief introduction to the use of PySCF to run a multitude of quantum chemical calculations. Starting with input parsing and uncorrelated Hartree-Fock theory, we'll incrementally touch upon how to use the majority of methods and features offered by PySCF through a number of simple examples, which all make reference to specific use cases within the dedicated `examples <https://github.com/pyscf/pyscf/tree/master/examples>`_ directory. Please also note that the cells below often share objects in-between one another.

.. _INPUT:

Input Parsing
=============

Molecules (or cells, cf. the :ref:`periodic <PBC>` section) are typically loaded in any of two ways (`gto/00-input_mole.py <https://github.com/pyscf/pyscf/blob/master/examples/gto/00-input_mole.py>`_):

  >>> from pyscf import gto
  >>> mol_h2o = gto.Mole()
  >>> mol_h2o.atom = '''O 0 0 0; H 0 1 0; H 0 0 1'''
  >>> mol_h2o.basis = 'ccpvdz'
  >>> mol_h2o.build()

or - using the convenient shortcut function - as  

  >>> mol_h2o = gto.M(atom = 'O 0 0 0; H 0 1 0; H 0 0 1', basis = 'ccpvdz')

Calling ``build()`` initializes a bunch of internal control parameters. Whenever you change the value of the attributes of :class:`Mole`, you need to call this function again to refresh the internal data of the object.

Symmetry may be specified in the ``Mole.symmetry`` attribute as either ``True`` or ``False`` (default is ``False``, i.e., off). Alternatively, a particular subgroup can be specified by a string argument (`gto/13-symmetry.py <https://github.com/pyscf/pyscf/blob/master/examples/gto/13-symmetry.py>`_):

  >>> mol_c2 = gto.M(atom = 'C 0 0 .625; O 0 0 -.625', symmetry = 'd2h')
  
Numerous other ways of inputting a molecular or crystalline geometry also exist (e.g., by means of Z-matrices or reading in .xyz files), cf. the complete suite of `gto <https://github.com/pyscf/pyscf/blob/master/examples/gto>`_ examples.

.. _MF:

Mean-Field Theory
=================

.. _HF:

Hartree-Fock
------------

A mean-field calculation is now trivially executed upon having initialized a :class:`Mole` object. For instance, a simple Hartree-Fock calculation on the H\ :sub:`2`\ O geometry of the :ref:`above <INPUT>` section will read (cf. `scf/00-simple_hf.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/00-simple_hf.py>`_):

  >>> from pyscf import scf
  >>> mf_h2o_rhf = scf.RHF(mol_h2o)
  >>> e_h2o = mf_h2o_rhf.kernel()

Besides the final converged ground-state energy, the mean-field object will further store the accompanying MO coefficients, occupations, etc. To illustrate how open-shell, possibly spin-polarized calculations are performed, different Hartree-Fock simulations of the O\ :sub:`2` dimer - with its triplet ground state - are given as (cf. `scf/02-rohf_uhf.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/02-rohf_uhf.py>`_):

  >>> mol_o2 = gto.M(atom='O 0 0 0; O 0 0 1.2', spin=2) # (n+2 alpha, n beta) electrons
  >>> mf_o2_uhf = scf.UHF(mol_o2)
  >>> mf_o2_uhf.kernel()
  >>> mf_o2_rohf = scf.ROHF(mol_o2)
  >>> mf_o2_rohf.kernel()

Finally, a second-order, Newton-Raphson SCF method has been implemented, which is also applicable with most XC functionals, cf. :ref:`below <KSDFT>`. This algorithm needs orthonormal orbitals and their corresponding occupancies as an initial guess:

  >>> mf_h2o_rhf = mf_h2o_rhf.newton()
  >>> e_h2o = mf_h2o_rhf.kernel()

.. _KSDFT:

Kohn-Sham Density Functional Theory
-----------------------------------

Running a KS-DFT calculation is every bit as straightforward as what we discussed :ref:`above <HF>` for an HF counterpart (cf. `dft/00-simple_dft.py <https://github.com/pyscf/pyscf/blob/master/examples/dft/00-simple_dft.py>`_):

  >>> from pyscf import dft
  >>> mf_h2o_rks = dft.RKS(mol) # likewise for UKS and ROKS
  >>> mf_h2o_rks.xc = 'b3lyp'
  
Besides the use of predefined XC functionals (cf. `pyscf/dft/libxc.py <https://github.com/pyscf/pyscf/blob/master/pyscf/dft/libxc.py>`_ and `pyscf/dft/xcfun.py <https://github.com/pyscf/pyscf/blob/master/pyscf/dft/xcfun.py>`_ for the complete list of
available functionals), these can also be fully defined by the user (`dft/24-custom_xc_functional.py <https://github.com/pyscf/pyscf/blob/master/examples/dft/24-custom_xc_functional.py>`_), as can the angular and radial grids (`dft/11-grid_scheme.py <https://github.com/pyscf/pyscf/blob/master/examples/dft/11-grid_scheme.py>`_):

  >>> mf_h2o_rks.xc = '.2 * HF + .08 * LDA + .72 * B88, .81 * LYP + .19 * VWN' # B3LYP
  >>> mf_h2o_rks.grids.atom_grid = (100, 770)
  >>> mf_h2o_rks.grids.prune = None
  
The use of a combination of dense and sparse grids are particularly important whenever XC functionals with non-local correlation calculation are employed (cf. `dft/33-nlc_functionals.py <https://github.com/pyscf/pyscf/blob/master/examples/dft/33-nlc_functionals.py>`_):

  >>> mf_c2_rks = dft.RKS(mol_c2)
  >>> mf_c2_rks.xc = 'wb97m_v'
  >>> mf_c2_rks.nlc = 'vv10'
  >>> mf_c2_rks.grids.atom_grid = (99,590)
  >>> mf_c2_rks.grids.prune = None
  >>> mf_c2_rks.nlcgrids.atom_grid = (50,194)
  >>> mf_c2_rks.nlcgrids.prune = dft.gen_grid.sg1_prune

.. _TDMF:

Time-Dependent Mean-Field Theory
--------------------------------

Linear response theory has been implemented for both HF and KS-DFT (cf. `tddft/00-simple_tddft.py <https://github.com/pyscf/pyscf/blob/master/examples/tddft/00-simple_tddft.py>`_):

  >>> from pyscf import tdscf
  >>> tdmf_h2o = tdscf.TDHF(mf_h2o_rhf) # or tdscf.TDDFT(mf_h2o_rks)
  >>> tdmf_h2o.nstates = 6
  >>> tdmf_h2o.kernel()

From a converged time-dependent mean-field calculation, the corresponding natural transition orbitals for a particular excited state may be recovered as (cf. `tddft/01-nto_analysis.py <https://github.com/pyscf/pyscf/blob/master/examples/tddft/01-nto_analysis.py>`_):

  >>> weights, nto = mytd.get_nto(state=2)
  
As an alternative to response theory, :math:`\Delta`-SCF with Gill's maximium occupation method has been implemented for calculating specific excited states, cf. `scf/50-mom-deltaSCF.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/50-mom-deltaSCF.py>`_.

.. _LOC:

Spatially Localized Molecular Orbitals
--------------------------------------

PySCF offers a number of different standard schemes for localizing MOs, e.g., Pipek-Mezey, Foster-Boys, and Edmiston-Ruedenberg (cf. `local_orb/03-split_localization.py <https://github.com/pyscf/pyscf/blob/master/examples/local_orb/03-split_localization.py>`_):

  >>> from pyscf import lo
  >>> occ_orbs = mf_h2o_rhf.mo_coeff[:, mf_h2o_rhf.mo_occ > 0.]
  >>> fb_h2o = lo.Boys(mol_h2o, occ_orbs, mf_h2o_rhf) # Foster-Boys
  >>> loc_occ_orbs = fb.kernel()
  >>> virt_orbs = mf_h2o_rhf.mo_coeff[:, mf_h2o_rhf.mo_occ == 0.]
  >>> pm_h2o = lo.Boys(mol_h2o, virt_orbs, mf_h2o_rhf) # Pipek-Mezey
  >>> loc_virt_orbs = pm.kernel()
  
In addition, Knizia's intrinsic bond orbitals have been implemented (cf. `local_orb/04-ibo_benzene_cubegen.py <https://github.com/pyscf/pyscf/blob/master/examples/local_orb/04-ibo_benzene_cubegen.py>`_):

  >>> iao = lo.iao.iao(mol, occ_orbs)
  >>> iao = lo.vec_lowdin(iao, mf_h2o_rhf.get_ovlp())
  >>> ibo = lo.ibo.ibo(mol, occ_orbs, iaos=iao)

.. _REL:

Relativistic Effects
--------------------

PySCF implements a Dirac-Hartree-Fock solver for including relativistic effects, in possible combination with Breit Gaunt interactions (cf. `scf/05-breit_gaunt.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/05-breit_gaunt.py>`_):

  >>> mf_c2_dhf = scf.DHF(mol_c2)
  >>> mf_c2_dhf.with_gaunt = True
  >>> mf_c2_dhf.with_breit = True
  >>> mf_c2_dhf.kernel()

As a popular alternative, scalar relativistic effects may be applied to a mean-field treatment by decorating the a :class:`SCF` object (either HF or KS-DFT) with the ``.x2c`` method (cf. `scf/21-x2c.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/21-x2c.py>`_), on top of which a correlated calculation may follow:

  >>> mf_o2_x2c = scf.UKS(mol_o2).x2c()
  >>> mf_o2_x2c.kernel()

.. _SYM:

Symmetry Handling
-----------------

Wave function symmetry may be explicitly controlled in an SCF calculation on the C\ :sub:`2` geometry of the :ref:`above <INPUT>` section by specifying frozen occupancy through the ``irrep_nelec`` attribute (`scf/13-symmetry.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/13-symmetry.py>`_):

  >>> mf_c2_rhf = scf.RHF(mol_c2)
  >>> mf_c2_rhf.irrep_nelec = {'Ag': 4, 'B1u': 4, 'B2u': 2, 'B3u': 2}
  >>> e_c2 = mf_c2_rhf.kernel()
  
Likewise, the final orbital symmetries may be probed from the MO coefficients (`symm/32-symmetrize_natural_orbital <https://github.com/pyscf/pyscf/blob/master/examples/symm/32-symmetrize_natural_orbital.py>`_):

  >>> from pyscf import symm
  >>> orbsym = symm.label_orb_symm(mol_c2, mol_c2.irrep_id, mol_c2.symm_orb, mf_c2_rhf.mo_coeff)

Integrals & Density Fitting
===========================

.. _INT:

1- and 2-Electron Integrals
---------------------------

A typical use case of for the integral code in PySCF is the integral transformation for a given set of orbitals to arrive at 1- and 2-electron integrals in a chosen MO basis, with the latter stored as (ij|kl) with 4-fold symmetry (cf. also `ao2mo/00-mo_integrals.py <https://github.com/pyscf/pyscf/blob/master/examples/ao2mo/00-mo_integrals.py>`_):

  >>> import numpy as np
  >>> from pyscf import ao2mo
  >>> hcore_ao = mol_h2o.intor_symmetric('int1e_kin') + mol_h2o.intor_symmetric('int1e_nuc')
  >>> hcore_mo = np.einsum('pi,pq,qj->ij', mf_h2o_rhf.mo_coeff, hcore_ao, mf_h2o_rhf.mo_coeff)
  >>> eri_4fold_ao = mol_h2o.intor('int2e_sph', aosym=4)
  >>> eri_4fold_mo = ao2mo.incore.full(eri_4fold_ao, mf_h2o_rhf.mo_coeff)
  
If desired, the transformed 2-electron integrals may also be saved to and read from a file in HDF5 format (`ao2mo/01-outcore.py <https://github.com/pyscf/pyscf/blob/master/examples/ao2mo/01-outcore.py>`_):

  >>> import tempfile
  >>> import h5py
  >>> ftmp = tempfile.NamedTemporaryFile()
  >>> ao2mo.kernel(mol_h2o, mf_h2o_rhf.mo_coeff, ftmp.name)
  >>> with h5py.File(ftmp.name) as f:
  >>>     eri_4fold = f['eri_mo']
  
User-defined Hamiltonians can also be used in PySCF, e.g., as input to a mean-field calculation and subsequent correlated treatment (`mcscf/40-customizing_hamiltonian.py <https://github.com/pyscf/pyscf/blob/master/examples/mcscf/40-customizing_hamiltonian.py>`_):

  >>> # 1D anti-PBC Hubbard model at half filling
  >>> n, u = 12, 2.
  >>> mol_hub = gto.M()
  >>> mol_hub.nelectron = n // 2
  >>> mol_hub.incore_anyway = True
  >>> h1 = np.zeros([n] * 2, dtype=np.float64)
  >>> for i in range(n-1):
  >>>     h1[i, i+1] = h1[i+1, i] = -1.
  >>> h1[n-1, 0] = h1[0, n-1] = -1.
  >>> eri = np.zeros([n] * 4, dtype=np.float64)
  >>> for i in range(n):
  >>>     eri[i, i, i, i] = u
  >>> mf_hub = scf.RHF(mol_hub)
  >>> mf_hub.get_hcore = lambda *args: h1
  >>> mf_hub.get_ovlp = lambda *args: np.eye(n)
  >>> mf_hub._eri = ao2mo.restore(8, eri, n) # 8-fold symmetry
  >>> mf_hub.init_guess = '1e'
  >>> mf_hub.kernel()

.. _DF:

Density Fitting Techniques
--------------------------

Density fitting of 2-electron integrals is most conveniently invoked by means of two main channels (cf. `df/00-with_df.py <https://github.com/pyscf/pyscf/blob/master/examples/df/00-with_df.py>`_):

  >>> mf_c2_rhf_df = mf_c2_rhf.density_fit(auxbasis='def2-universal-jfit') # option 1
  >>> from pyscf import df
  >>> mf_c2_rhf_df = df.density_fit(scf.RHF(mol_c2), auxbasis='def2-universal-jfit') # option 2
  
In the former of these two option, decoration by the ``scf.density_fit`` function generates a new object that works in exactly the
same way as the regular :class:`SCF` object, but which is entirely independent of the original ``mf_c2_rhf`` object.

For a discussion on how to use density fitting alongside the :ref:`Newton-Raphson SCF algorithm <HF>` and :ref:`scalar relativistic effects <REL>`, please see `scf/23-decorate_scf.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/23-decorate_scf.py>`_.

Correlated Wave Function Theory
===============================

.. _MPCCADC:

Perturbation Theory, Coupled Cluster, and Algebraic Diagrammatic Constructions
------------------------------------------------------------------------------

PySCF offers both second-order Møller-Plesset, coupled cluster, and algebraic diagrammatic construction functionalities. The former of these are are implemented both with and without :ref:`density fitting <DF>`, again depending on the ``with_df`` attribute of the underlying mean-field object (cf. `mp/00-simple_mp2.py <https://github.com/pyscf/pyscf/blob/master/examples/mp/00-simple_mp2.py>`_):

  >>> from pyscf import mp
  >>> mp2_c2 = mp.MP2(mf_c2_rhf)
  >>> e_c2 = mp2_c2.kernel()[0]
  >>> mp2_c2_df = mp.MP2(mf_c2_rhf_df)
  >>> e_c2_df = mp2_c2_df.kernel()[0]
  
At the coupled cluster level of theory, CCD, CCSD, and CCSD(T) calculation can be performed for both closed- and open-shell systems (cf. `cc/00-simple_ccsd_t.py <https://github.com/pyscf/pyscf/blob/master/examples/cc/00-simple_ccsd_t.py>`_):

  >>> from pyscf import cc
  >>> ccsd_h2o = cc.CCSD(mf_h2o_rhf)
  >>> ccsd_h2o.direct = True # AO-direct algorithm to reduce I/O overhead
  >>> ccsd_h2o.frozen = 1 # frozen core
  >>> e_ccsd = ccsd_h2o.kernel()[1]
  >>> e_ccsd_t = e_ccsd + ccsd_h2o.ccsd_t()

As for MP2, this CCSD calculation will employ density fitting depending on the respective settings of ``mf_h2o_rhf``. This is also for subsequent EOM-CCSD calculations (cf. `cc/20-ip_ea_eom_ccsd.py <https://github.com/pyscf/pyscf/blob/master/examples/cc/20-ip_ea_eom_ccsd.py>`_):

  >>> e_ip_ccsd = ccsd_h2o.ipccsd(nroots=1)[0]
  >>> e_ea_ccsd = ccsd_h2o.eaccsd(nroots=1)[0]
  >>> e_ee_ccsd = ccsd_h2o.eeccsd(nroots=1)[0]
  
Finally, the ADC(2), ADC(2)-X, and ADC(3) schemes have all been implemented using a similar API (cf. `adc/01-closed_shell.py <https://github.com/pyscf/pyscf/blob/master/examples/adc/01-closed_shell.py>`_):

  >>> from pyscf import adc
  >>> adc_h2o = adc.ADC(mf_h2o_rhf)
  >>> e_ip_adc2 = adc_h2o.kernel()[0] # IP-ADC(2) for 1 root
  >>> adc_h2o.method = "adc(2)-x"
  >>> adc_h2o.method_type = "ea"
  >>> e_ea_adc2x = adc_h2o.kernel()[0] # EA-ADC(2)-x for 1 root
  >>> adc_h2o.method = "adc(3)"
  >>> adc_h2o.method_type = "ea"
  >>> e_ea_adc3 = adc_h2o.kernel(nroots = 3)[0] # EA-ADC(3) for 3 roots

Please note that all of these codes are written in pure Python (with calls to BLAS) and neither of them make use of point group symmetry.

.. _FCI:

Full Configuration Interaction
------------------------------

In contrast to the correlation methods discussed :ref:`above <MPCCADC>`, PySCF offer a number of powerful kernels (written in optimized C) for performing exact diagonalization of all kinds of Hamiltonians and systems of arbitrary spin. For standard cases, in which all electrons of a given systems are correlated among all MOs, the syntax follows that of other correlation methods for closed- and open-shell systems (cf. `fci/00-simple_fci.py <https://github.com/pyscf/pyscf/blob/master/examples/fci/00-simple_fci.py>`_):

  >>> from pyscf import fci
  >>> fci_h2o = fci.FCI(mf_h2o_rhf)
  >>> e_fci = fci_h2o.kernel()[0]
  
However, the various FCI solvers (tabulated in `pyscf/fci/__init__.py <https://github.com/pyscf/pyscf/blob/master/pyscf/fci/__init__.py>`_) further allow for user-defined 1- and 2-electron Hamiltonians (cf. `fci/01-given_h1e_h2e.py <https://github.com/pyscf/pyscf/blob/master/examples/fci/01-given_h1e_h2e.py>`_):

  >>> fs = fci.direct_spin1.FCI() # direct_spin0 instead for singlet system ground states
  >>> e, fcivec = fs.kernel(h1, h2, N, 8) # 8 electrons in N orbitals
  >>> e, fcivec = fs.kernel(h1, h2, N, (5,4))  # (5 alpha, 4 beta) electrons
  >>> e, fcivec = fs.kernel(h1, h2, N, (3,1))  # (3 alpha, 1 beta) electrons
  
The individual solvers can yield more than a single (ground) states by setting ``fs.nroots > 1``, and 1- to 4-electron density matrices, alongside 1- and 2-electron transition density matrices, can be computed at differing cost (cf. `fci/14-density_matrix.py <https://github.com/pyscf/pyscf/blob/master/examples/fci/14-density_matrix.py>`_):

  >>> rdm1 = fs.make_rdm1(fcivec, N, (5, 4)) # spin-traced 1-electron density matrix
  >>> rdm1a, rdm1b = fs.make_rdm1s(fcivec, norb, (5, 4)) # alpha and beta 1-electron density matrices
  >>> t_rdm1 = fs.trans_rdm1(fcivec0, fcivec1, N, (5, 4)) # spin-traced 1-electron transition density matrix
  
In addition, the FCI code is accompanied by a wealth of library tools for inspecting individual wave function expansions, assigning spin states and :ref:`symmetry <SYM>`, explicitly constructing the full Hamiltonian, etc.

Multiconfigurational Methods
============================

.. _CAS:

Complete Active Space Configuration Interaction & Self-Consistent Field
-----------------------------------------------------------------------

The powerful FCI solvers discussed :ref:`above <FCI>` further act as engines for the various complete active space methods in PySCF, which all share a similar API in common (cf. `mcscf/00-simple_casci.py <https://github.com/pyscf/pyscf/blob/master/examples/mcscf/00-simple_casci.py>`_ & `mcscf/00-simple_casscf.py <https://github.com/pyscf/pyscf/blob/master/examples/mcscf/00-simple_casscf.py>`_):

  >>> from pyscf import mcscf
  >>> casci_h2o = mcscf.CASCI(mf_h2o_rhf, 6, 8)
  >>> e_casci = casci_h2o.kernel()[0]
  >>> casscf_h2o = mcscf.CASSCF(mf_h2o_rhf, 6, 8)
  >>> e_casscf = casscf_h2o.kernel()[0]

While CASCI calculations may also be performed using the FCI codes, the API around the CASCI and CASSCF classes allows for easier control over aspects, such as, :ref:`density fitting of 2-electron integrals <DF>` (cf. `mcscf/16-density_fitting.py <https://github.com/pyscf/pyscf/blob/master/examples/mcscf/16-density_fitting.py>`_) and the standard frozen-core approximation (cf. `mcscf/19-frozen_core.py <https://github.com/pyscf/pyscf/blob/master/examples/mcscf/19-frozen_core.py>`_):

  >>> casscf_h2o_df = mcscf.DFCASSCF(mf_h2o_rhf, 6, 8, auxbasis='ccpvtzfit')
  >>> casscf_h2o_df.frozen = 1 # frozen core
  >>> e_casscf_df = casscf_h2o_df.kernel()[0]

In the case of CASSCF calculations, these may be performed in a state-specific or state-averaged manner (cf. `mcscf/41-state_average.py <https://github.com/pyscf/pyscf/blob/master/examples/mcscf/41-state_average.py>`_):

  >>> casscf_c2 = mcscf.CASSCF(mf_c2_rhf, 8, 8)
  >>> solver_t = fci.direct_spin1_symm.FCI(mol_c2)
  >>> solver_t.spin = 2
  >>> solver_t.nroots = 1
  >>> solver_t = fci.addons.fix_spin(solver_t, shift=.2, ss=2) # 1 triplet
  >>> solver_s = fci.direct_spin0_symm.FCI(mol_c2) # 2 singlets
  >>> solver_s.spin = 0
  >>> solver_s.nroots = 2
  >>> mcscf.state_average_mix_(casscf_c2, [solver_t, solver_s], np.ones(3) / 3.)
  >>> casscf_c2.kernel()
  
Finally, additional dynamic correlation may be added by means of second-order perturbation theory in the form of NEVPT2 (cf. `mrpt/02-cr2_nevpt2/cr2-scan.py <https://github.com/pyscf/pyscf/blob/master/examples/mrpt/02-cr2_nevpt2/cr2-scan.py>`_):

  >>> from pyscf import mrpt
  >>> e_nevpt2 = mrpt.NEVPT(casscf_h2o).kernel()

.. _EXTFCI:

External Approximate Full Configuration Interaction Solvers
-----------------------------------------------------------

Besides the exact solvers discussed :ref:`earlier <FCI>`, a number of highly efficient approximate solvers for use in CASCI and CASSCF calculations may also be employed via their interfaces in PySCF. For instance, the `StackBlock <https://github.com/sanshar/StackBlock>`_ code can be used as an optimized DMRG solver to perform parallel DMRGSCF calculations across several processes (cf. `dmrg/01-dmrg_casscf_with_stackblock.py <https://github.com/pyscf/pyscf/blob/master/examples/dmrg/01-dmrg_casscf_with_stackblock.py>`_):

  >>> from pyscf import dmrgscf
  >>> import os
  >>> from pyscf.dmrgscf import settings
  >>> if 'SLURMD_NODENAME' in os.environ: # slurm system
  >>>     settings.MPIPREFIX = 'srun'
  >>> elif 'PBS_NODEFILE' in os.environ: # PBS system
  >>>     settings.MPIPREFIX = 'mpirun'
  >>> else: # MPI on single node
  >>>     settings.MPIPREFIX = 'mpirun -np 4'
  >>> mc = dmrgscf.DMRGSCF(mf_c2_rhf, 8, 8)
  >>> mc.state_average_([.5] * 2)
  >>> mc.fcisolver.memory = 4 # in GB
  >>> mc.fcisolver.num_thrds = 8 # number of threads to spawn on each MPI process
  >>> e_dmrgscf = mc.kernel()
  
Likewise, similar interfaces furthermore exist to enable the execution of i-FCIQMC (using `NECI <https://github.com/ghb24/NECI_STABLE>`_) and SHCI (using either `Dice <https://github.com/sanshar/Dice>`_ or `Arrow <https://github.com/QMC-Cornell/shci>`_) calculations from within PySCF.

Geometry Optimization Techniques
================================

Solvent Effects
===============

Polarizable Continuum Methods
-----------------------------

Quantum Mechanics/Molecular Mechanics Methods
---------------------------------------------

Semi-Empirical Methods
======================

.. _PBC:

Periodic Boundary Conditions
============================

df/00-with_df.py, pbc/11-gamma_point_all_electron_scf.py

Miscellaneous Library Tools
===========================


