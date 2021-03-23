
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

  >>> mol_c2 = gto.M(atom = 'C 0 0 0.625; O 0 0 -0.625', symmetry = 'd2h')
  
Numerous other ways of inputting a molecular or crystalline geometry also exist (e.g., by means of Z-matrices or from files), cf. the complete suite of `gto <https://github.com/pyscf/pyscf/blob/master/examples/gto>`_ examples.

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

Finally, a second-order SCF method has been implemented, which is also applicable with most XC functionals, cf. :ref:`below <KSDFT>`. This algorithm needs orthonormal orbitals and their corresponding occupancies as an initial guess:

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
  >>> mol_hub = gto.M()
  >>> mol_hub.nelectron = 6
  >>> mol_hub.incore_anyway = True
  >>> n = 12
  >>> h1 = np.zeros([n] * 2, dtype=np.float64)
  >>> for i in range(n-1):
  >>>     h1[i, i+1] = h1[i+1, i] = -1.
  >>> h1[n-1, 0] = h1[0, n-1] = -1.
  >>> eri = np.zeros([n] * 4, dtype=np.float64)
  >>> for i in range(n):
  >>>     eri[i, i, i, i] = 2.
  >>> mf_hub = scf.RHF(mol_hub)
  >>> mf_hub.get_hcore = lambda *args: h1
  >>> mf_hub.get_ovlp = lambda *args: np.eye(n)
  >>> mf_hub._eri = ao2mo.restore(8, eri, n) # 8-fold symmetry
  >>> mf_hub.init_guess = '1e'
  >>> mf_hub.kernel()

.. _DF:

Density Fitting Techniques
--------------------------

Correlated Wave Function Theory
===============================

Møller-Plesset Perturbation Theory
----------------------------------

We can compute the correlation energy at the second-order
Møller-Plesset level of theory with :mod:`mp.mp2`::

  >>> from pyscf import mp
  >>> mp2 = mp.MP2(m)
  >>> print('E(MP2) = %.9g' % mp2.kernel()[0])
  E(MP2) = -0.379359288

Coupled Cluster
---------------

Algebraic Diagrammatic Construction
-----------------------------------

Full Configuration Interaction
------------------------------

Multiconfigurational Methods
============================

Complete Active Space Configuration Interaction
-----------------------------------------------

CASCI and CASSCF calculations can be run with similar inputs::

  >>> from pyscf import mcscf
  >>> mc = mcscf.CASCI(m, 4, 6)
  >>> print('E(CASCI) = %.9g' % mc.casci()[0])
  E(CASCI) = -149.601051
  >>> mc = mcscf.CASSCF(m, 4, 6)
  >>> print('E(CASSCF) = %.9g' % mc.kernel()[0])
  E(CASSCF) = -149.613191

In this example, the CAS space is (6e, 4o), that is, six electrons in
four orbitals.

Complete Active Space Self-Consistent Field
-------------------------------------------

Density Matrix Renormalization Group
------------------------------------

Full Configuration Interaction Quantum Monte Carlo
--------------------------------------------------

Multireference Perturbation Theory
----------------------------------

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

Miscellaneous Library Tools
===========================


