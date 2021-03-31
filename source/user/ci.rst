.. _theory_ci:

Configuration interaction (CISD and FCI)
****************************************
*Modules*: :mod:`ci`, :mod:`fci`, :mod:`pbc.ci`

PySCF has separate molecular implementations for configuration
interaction singles and doubles (CISD) (:mod:`ci`) and full
configuration interaction (FCI) (:mod:`fci`). The functionalities of the
CISD implementation are similar to the functionalities of MP2
(:ref:`user_mp2`) and CCSD (:ref:`user_cc`) in PySCF.


Introduction
============

Configuration interaction is a post-Hartree--Fock method, diagonalising
the many-electron Hamiltonian matrix. FCI includes all Slater determinants of
appropriate symmetry in the eigenvector basis. CISD includes only those
that differ from the Hartree--Fock determinant by at most two occupied
spinorbitals.


Configuration interaction singles and doubles (CISD)
====================================================

A simple example (see :source:`examples/ci/00-simple_cisd.py`) of running
a restricted and an unrestricted CISD calculation is

.. literalinclude:: ../../examples/ci/00-simple_cisd.py

which outputs

.. code::

  converged SCF energy = -99.9873974403488
  E(RCISD) = -100.196182975762  E_corr = -0.2087855354132202
  RCISD correlation energy -0.2087855354132202

namely, the (restricted) Hartree--Fock energy, the RCISD energy and their
difference, the RCISD correlation energy, as well as

.. code::

  converged SCF energy = -99.9873974403486  <S^2> = 5.3734794e-13  2S+1 = 1
  E(UCISD) = -100.196182975774  E_corr = -0.2087855354254317
  UCISD correlation energy -0.20878553542543174

namely, the (unrestricted) Hartree--Fock energy, the UCISD energy and
their difference, the UCISD correlation energy.

.. note::

  .. code::

    mycc = mf.CISD().run()
  
  could have been replaced by

  .. code::

    myci = pyscf.ci.CISD(mf)
    myci.kernel()

  or

  .. code::

    myci = pyscf.ci.CISD(mf).run()
  
  for the same result.

.. note::

  Density fitting is not implemented for CISD.


Spin symmetry
-------------

The CISD module in PySCF supports a number of reference wavefunctions with
broken spin symmetry.  In particular, CISD can be performed with a
spin-restricted, spin-unrestricted, and general (spin-mixed) Hartree-Fock
solution, leading to the RCISD, UCISD, and GCISD methods. The small
example above shows this for RCISD and UCISD.

The module-level ``ci.CISD(mf)`` constructor can infer the correct method based
on the level of symmetry-breaking in the mean-field argument.  For more explicit
control or inspection, the respective classes and functions can be found in
``cisd.py`` (restricted), ``ucisd.py`` (unrestricted), and ``gcisd.py``
(general).

For example, a spin-unrestricted calculation on triplet oxygen can be performed
as follows::

    from pyscf import gto, scf, ci
    mol = gto.M(
        atom = 'O 0 0 0; O 0 0 1.2',  # in Angstrom
        basis = 'ccpvdz',
        spin = 2
    )
    mf = scf.HF(mol).run() # this is UHF
    myci = ci.CISD(mf).run() # this is UCISD
    print('UCISD total energy = ', myci.e_tot)


Properties
----------

A number of properties are available at the CISD level.

Unrelaxed 1- and 2-electron reduced density matrices can be calculated. 
They are returned in the MO basis::

    dm1 = myci.make_rdm1()
    dm2 = myci.make_rdm2()

Analytical nuclear gradients can be calculated
(see :cite:`Osamura1987,Yamaguchi2011` and references therein)::

    mygrad = myci.nuc_grad_method().run()


Frozen orbitals
---------------

By default, CISD calculations in PySCF correlate all electrons in all available
orbitals. To freeze the lowest-energy core orbitals,
use the ``frozen`` keyword argument::

    myci = ci.CISD(mf, frozen=2).run()

To freeze occupied and/or unoccupied orbitals with finer control, a list of
0-based orbital indices can be provided as the ``frozen`` keyword argument::
    
    # freeze 2 core orbitals
    myci = ci.CISD(mf, frozen=[0,1]).run()
    # freeze 2 core orbitals and 3 unoccupied orbitals
    myci = ci.CISD(mf, frozen=[0,1,16,17,18]).run()


Wavefunction overlap
--------------------

The following example shows how to evaluate the overlap of two different
CISD wavefunctions.

.. literalinclude:: ../../examples/ci/32-wfn_overlap.py


References
==========

.. bibliography:: ref_ci.bib
   :style: unsrt
