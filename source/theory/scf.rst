.. _theory_scf:

***********************************
Self-consistent field (SCF) methods
***********************************

Introduction
============
In self-consistent field (SCF) methods, the electron interactions are treated in a mean-field way.
Here, the SCF methods include both Hartree-Fock (HF) theory and Kohn-Sham (KS) density functional theory (DFT).
This Chapter summarizes the general SCF capabilities of PySCF. 
For more details specific to DFT, see :ref:`theory_dft`.

A minimal example of using the SCF module is as follows::

    from pyscf import gto, scf
    mol = gto.M(
        atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
        basis = 'ccpvdz',
        symmetry = True,
    )
    myhf = scf.HF(mol)
    myhf.kernel()

This will run a HF calculation for the hydrogen fluoride molecule using the default SCF settings.


Theory
======




Initial guess
=============




Converging SCF calculations
===========================
