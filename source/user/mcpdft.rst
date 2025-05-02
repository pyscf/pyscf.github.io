.. _user_mcpdft:

Multi-configuration pair-density functional theory (MC-PDFT)
************************************************************

*Modules*: :py:mod:`pyscf.mcpdft`


Introduction
------------

Multiconfiguration pair-density functional theory (MC-PDFT) refers to methods in which the total electronic energy is obtained or derived from a functional of the on-top pair density as well as the total electron density, and these densities are in turn obtained from a multiconfigurational self-consistent field (MCSCF) wave function or wave functions of some sort.
The on-top pair density is the probability of two electrons existing simultaneously at one point in space.
The non-classical part of the MC-PDFT energy is called the "on-top energy" and is directly analogous to the exchange-correlation energy in standard density functional theory.
In practice, the on-top energy functional is always a generalization or a "translation" of an exchange-correlation functional from Kohn-Sham DFT.
MC-PDFT furnishes a practical and effective generalization of KS-DFT functionals from single Slater determinants to general multiconfigurational wave functions because the pair density, unlike the spin density, can simultaneously encode both strong electron correlation effects and spin symmetry constraints.
PySCF implements the MC-PDFT protocol of (CITE), in which the molecular orbital and/or configuration interaction coefficients are chosen to minimize the expectation value of the molecular Hamiltonian in the underlying MCSCF wave function, not the MC-PDFT energy itself.

The relationship of :mod:`mcpdft` to :mod:`mcscf` is roughly similar to the relationship between :mod:`dft` and :mod:`scf`, in that instances of MC-PDFT methods are also instances of MCSCF methods.
:mod:`mcpdft` thus inherits :mod:`mcscf`'s compatibility with density fitting (:ref:`user_df`) and x2c (:ref:`user_x2c`) as well as CASSCF's ``state_average_`` and ``state_average_mix_`` extensions.
As in :mod:`dft`, the MC-PDFT object has a ``grids`` attribute which stores information about the quadrature grid.
The on-top functional is identified with a string stored in the ``otxc`` attribute of the MC-PDFT object, which is passed as the required second argument to the method constructor; see below.

CASCI-PDFT
""""""""""
.. code-block:: python

  from pyscf import gto, mp, mcpdft
  mol = gto.M(
      atom = 'O 0 0 0; O 0 0 1.2',
      basis = 'ccpvdz',
      spin = 2)
  myhf = mol.RHF().run()
  # Use MP2 natural orbitals to define the active space for the single-point CAS-CI calculation
  mymp = mp.UMP2(myhf).run()

  noons, natorbs = mcscf.addons.make_natural_orbitals(mymp)
  ncas, nelecas = (6,8)
  otfnal = 'tPBE'
  mycas = mcpdft.CASCI(myhf, otfnal, ncas, nelecas)
  mycas.kernel(natorbs)


CASSCF-PDFT
"""""""""""
.. code-block:: python

  import pyscf
  from pyscf import mcpdft
  mol = pyscf.M(
      atom = 'O 0 0 0; O 0 0 1.2',
      basis = 'ccpvdz',
      spin = 2)
  myhf = mol.RHF().run()
  ncas, nelecas = (6,(5,3))
  otfnal = 'tPBE'
  # For convenience, mycas.grids.level can be set at construction using an optional kwarg
  mycas = mcpdft.CASSCF(myhf, otfnal, ncas, nelecas, grids_level=6).run()

When state-averaged CASSCF wave functions are used, the MC-PDFT energy is computed separately for each state.

Available on-top functionals
----------------------------

t, ft

Multi-state extension methods
-----------------------------

CMS-PDFT, L-PDFT




