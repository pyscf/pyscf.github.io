-.. _getting_started:


Tutorial
********

This tutorial shows simple examples of using PySCF to run quantum
chemical calculations. An Ipython notebook of user-guide can be found
at https://github.com/nmardirossian/PySCF_Tutorial.



Simple examples
===============


Hartree-Fock calculations
-------------------------

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


An MP2 calculation
------------------

We can compute the correlation energy at the second-order
Møller-Plesset level of theory with :mod:`mp.mp2`::

  >>> from pyscf import mp
  >>> mp2 = mp.MP2(m)
  >>> print('E(MP2) = %.9g' % mp2.kernel()[0])
  E(MP2) = -0.379359288


Coupled-cluster calculations
----------------------------


Density functional calculations
-------------------------------

CASCI and CASSCF
----------------

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


