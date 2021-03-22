
Quickstart
**********

The present tutorial is meant to provide a brief introduction to the use of PySCF to run a multitude of quantum chemical calculations. Starting with input parsing and uncorrelated Hartree-Fock theory, we'll incrementally touch upon how to use the majority of methods and features offered by PySCF through a number of simple examples, which all make reference to specific use cases within the dedicated `examples <https://github.com/pyscf/pyscf/tree/master/examples>`_ directory. Please also note that the cells below often share objects in-between one another.

Input Parsing
=============

Molecules (or cells, cf. the :ref:`periodic<PBC>` section) are typically loaded in any of two ways:

  >>> from pyscf import gto, scf
  >>> mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='ccpvdz')

Mean-Field Theory
=================

Hartree-Fock
------------

Here is an example to run a HF calculation on the hydrogen molecule::

  >>> from pyscf import gto, scf
  >>> mol = gto.M(atom='H 0 0 0; H 0 0 0.74', basis='ccpvdz')
  >>> mf = scf.RHF(mol)
  >>> mf.kernel()
  converged SCF energy = -1.12870009355644
  -1.1287000935564409

The coordinates are given by default in Ångström. One can also define
the geometry in atomic units with::

  >>> from pyscf import gto, scf
  >>> mol = gto.M(atom='H 0 0 0; H 0 0 1.4', unit='bohr', basis='ccpvdz')
  >>> mf = scf.RHF(mol)
  >>> mf.kernel()
  converged SCF energy = -1.12870944897989
  -1.1287094489798912

As an example of spin-polarized calculations, let's take the oxygen
molecule, O2. O2 has a triplet ground state, so we need to set the
"spin" to ``2`` (a triplet has two more alpha electrons than beta
electrons)::

  >>> from pyscf import gto, scf
  >>> mol = gto.M(atom='O 0 0 0; O 0 0 1.2', basis='ccpvdz', spin=2)
  >>> mf = scf.UHF(mol)
  >>> mf.kernel()
  converged SCF energy = -149.62899231417  <S^2> = 2.0326472  2S+1 = 3.0216864
  -149.62899231416952

Molecular point-group symmetry is not used by default; however, you
can turn it on with the ``symmetry`` variable::

  >>> from pyscf import gto, scf
  >>> mol = gto.M(atom='O 0 0 0; O 0 0 1.2', basis='ccpvdz', spin=2, symmetry=True)
  >>> mf = scf.UHF(mol)
  >>> mf.kernel()
  converged SCF energy = -149.628992314169  <S^2> = 2.0326472  2S+1 = 3.0216864
  -149.6289923141695

Kohn-Sham Density Functional Theory
-----------------------------------

Time-Dependent Mean-Field Theory
--------------------------------

Spatially Localized Molecular Orbitals
--------------------------------------

Relativistic Effects
--------------------

Symmetry Handling
-----------------

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


