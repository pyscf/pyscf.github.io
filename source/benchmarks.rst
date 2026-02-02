.. _benchmarks:

Benchmarks of PySCF
*******************

Below are the timings of a suite of benchmark applications for PySCF. These
results are provided as a rough guide for the timings that you can expect when
running PySCF. All calculations were run on [insert hardware description].

.. warning::

  This page is a work in progress.

Hartree-Fock theory
===================

.. list-table:: (H\ :sub:`2`\ O)\ :sub:`5`\  / cc-pvdz
   :widths: 25 75
   :header-rows: 0

   * - Total wall time (s)
     - 0.24
   * - Number of threads 
     - 12
   * - Number of electrons
     - 50
   * - Number of AOs
     - 120
   * - Input file
     - path/to/input.py
   * - Geometry file
     - molecules/water_clusters/003.xyz
   * - PySCF version 
     - 2.9.0
   * - BLAS for NumPy
     - scipy-openblas
   * - BLAS for PySCF
     - libmkl_intel_lp64.so, libmkl_sequential.so, libmkl_core.so 

Time-dependent Hartree-Fock theory
==================================


Density functional theory
=========================

