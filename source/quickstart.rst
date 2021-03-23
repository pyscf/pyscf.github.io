
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

Calling ``build()`` initializes a bunch of internal control parameters. Whenever you change the value of the attributes of :class:`Mole <https://github.com/pyscf/pyscf/blob/master/pyscf/gto/mole.py#L1941>`, you need to call this function again to refresh the internal data of the object.

Symmetry may be specified in the ``Mole.symmetry`` attribute as either ``True`` or ``False`` (default is ``False``, i.e., off). Alternatively, a particular subgroup can be specified by a string argument (`gto/13-symmetry.py <https://github.com/pyscf/pyscf/blob/master/examples/gto/13-symmetry.py>`_):

  >>> mol_c2 = gto.M(atom = 'C 0 0 0.625; O 0 0 -0.625', symmetry = 'd2h')
  
Numerous other ways of inputting a molecular or crystalline geometry also exist (e.g., by means of Z-matrices or from files), cf. the complete suite of `gto <https://github.com/pyscf/pyscf/blob/master/examples/gto>`_ examples.

.. _MF:

Mean-Field Theory
=================

.. _HF:

Hartree-Fock
------------

A mean-field calculation is now trivially executed upon having initialized a :class:`Mole <https://github.com/pyscf/pyscf/blob/master/pyscf/gto/mole.py#L1941>` object. For instance, a simple Hartree-Fock calculation on the H\ :sub:`2`\ O geometry of the :ref:`above <INPUT>` section will read (cf. `scf/00-simple_hf.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/00-simple_hf.py>`_):

  >>> from pyscf import scf
  >>> mf_h2o_rhf = scf.RHF(mol_h2o)
  >>> e_h2o = mf_h2o_rhf.kernel()

Besides the final converged ground-state energy, the :class:`SCF <https://github.com/pyscf/pyscf/blob/master/pyscf/scf/hf.py#L1308>` object will further store the accompanying MO coefficients, occupations, etc. To illustrate how open-shell, possibly spin-polarized calculations are performed, different Hartree-Fock simulations of the O\ :sub:`2` dimer - with its triplet ground state - are given as (cf. `scf/02-rohf_uhf.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/02-rohf_uhf.py>`_):

  >>> from pyscf import gto, scf
  >>> mol_o2 = gto.M(atom='O 0 0 0; O 0 0 1.2', spin=2) # (n+2,n) electrons of (alpha,beta) spin
  >>> mf_o2_uhf = scf.UHF(mol_o2)
  >>> mf_o2_uhf.kernel()
  >>> mf_o2_rohf = scf.ROHF(mol_o2)
  >>> mf_o2_rohf.kernel()

.. _KSDFT:

Kohn-Sham Density Functional Theory
-----------------------------------

Running a KS-DFT calculation is every bit as straightforward as what we discussed :ref:`above <HF>` for an HF counterpart (cf. `dft/00-simple_dft.py <https://github.com/pyscf/pyscf/blob/master/examples/dft/00-simple_dft.py>`_):

  >>> from pyscf import dft
  >>> mf_h2o_rks = dft.RKS(mol) # likewise for UKS and ROKS
  >>> mf_h2o_rks.xc = 'b3lyp'
  
Besides the use of predefined XC functionals (cf. `pyscf/dft/libxc.py <https://github.com/pyscf/pyscf/blob/master/pyscf/pyscf/dft/libxc.py>`_ and `pyscf/dft/xcfun.py <https://github.com/pyscf/pyscf/blob/master/pyscf/pyscf/dft/xcfun.py>`_ for the complete list of
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

.. _LOC:

Spatially Localized Molecular Orbitals
--------------------------------------

.. _REL:

Relativistic Effects
--------------------

.. _SYM:

Symmetry Handling
-----------------

Wave function symmetry may be explicitly controlled in an SCF calculation on the C\ :sub:`2` geometry of the :ref:`above <INPUT>` section by specifying frozen occupancy through the ``irrep_nelec`` attribute (`scf/13-symmetry.py <https://github.com/pyscf/pyscf/blob/master/examples/scf/13-symmetry.py>`_):

  >>> mf_c2 = scf.RHF(mol_c2)
  >>> mf_c2.irrep_nelec = {'Ag': 4, 'B1u': 4, 'B2u': 2, 'B3u': 2}
  >>> e_c2 = mf_c2.kernel()
  
Likewise, the final orbital symmetries may be probed from the MO coefficients (`symm/32-symmetrize_natural_orbital <https://github.com/pyscf/pyscf/blob/master/examples/symm/32-symmetrize_natural_orbital>`_):

  >>> from pyscf import symm
  >>> orbsym = symm.label_orb_symm(mol_c2, mol_c2.irrep_id, mol_c2.symm_orb, mf_c2.mo_coeff)

Integrals & Density Fitting
===========================

1- and 2-Electron Integrals
---------------------------

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


