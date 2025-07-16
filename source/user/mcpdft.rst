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
PySCF implements the MC-PDFT protocol of Reference :cite:`Manni2014`, in which the molecular orbital and/or configuration interaction coefficients are chosen to minimize the expectation value of the molecular Hamiltonian in the underlying MCSCF wave function, not the MC-PDFT energy itself.

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
When either state-averaged CASSCF or multiple roots with CASCI are calculated, the ordinal label of each state is as assigned by the underlying MCSCF calculation, regardless of the relative MC-PDFT energies.
In other words, state 2 may sometimes have a lower MC-PDFT energy than state 1.

Available on-top functionals
----------------------------

"Translated" on-top functionals corresponding to any pure LDA or GGA exchange-correlation functional :cite:`Manni2014`, as well as meta-GGA functionals :cite:`Bao2025` which depend on the kinetic energy, are available.
These are specified by prepending the letter ``t`` to the name of the base DFT functional, as in ``tPBE`` or ``tBLYP``.
"Fully-translated" on-top functionals :cite:`Carlson2015`  corresponding to pure LDA or GGA functionals are also available and are specified by prepending the letters ``ft`` to the name of the base functional.

With the exception of PBE0, for which the translated and fully-translated generalizations are implemented precisely as stipulated in Reference :cite:`Pandharkar2020`, prepending ``t`` or ``ft`` to the names of hybrid DFT functionals in this way will not work, because there are multiple forms of hybrid on-top functional in the literature (for instance, compare References :cite:`Pandharkar2020` and :cite:`Mostafanejad2020`).
However, it is still possible to specify a variety of global hybrid on-top functionals:

.. code-block:: python

  from pyscf import mcpdft

  # Construct a custom hybrid functional with 20% CASSCF energy (JPCL 11, 10158, 2020) 
  
  my_otxc = 't' + mcpdft.hyb ('BLYP', .2)
  mc = mcpdft.CASSCF (mf, my_otxc, 6, 8).run ()
  
  # "lambda-MC-PDFT" of JCTC 16, 2274, 2020.
  
  my_otxc = 't' + mcpdft.hyb('PBE', .5, hyb_type='lambda')
  mc.run (otxc=my_otxc)


Multi-state extension methods
-----------------------------

MC-PDFT evaluates the total energy of a single electronic state through a nonlinear functional expression involving density variables.
For approximate functionals evaluated at approximate densities, it cannot be guaranteed that potential energy surfaces computed in this way for states of the same symmetry only cross along `3N-8`-dimensional conical intersection seams, which is a requirement for qualitatively accurate non-adiabatic molecular dynamics simulations.
Multi-state methods are extensions of MC-PDFT which address this difficulty by defining the electronic energies of a set of states as eigenvalues of a small Hamiltonian matrix, of which the elements are related in one way or another to the MC-PDFT energy functional expression.
In PySCF, multi-state MC-PDFT methods are implemented as child classes of state-averaged CASSCF methods.

The preferred multi-state MC-PDFT method is L-PDFT :cite:`Hennefarth2023`, which is performed in the following way:

.. code-block:: python

  from pyscf import mcpdft

  mc = mcpdft.CASSCF(mf, 'tPBE', 2, 2)
  mc.fix_spin_(ss=0) # often necessary!
  mc_ms = mc.multi_state ([.5, .5], method="lin")
  mc_ms.kernel () # Performs SA-CASSCF orbital optimization followed by L-PDFT

Note that the first argument of ``multi_state`` is a sequence of weights, exactly as in ``state_average_``.
The second argument stipulates the type of multi-state MC-PDFT.
XMS-PDFT :cite:`Bao2020_XMS` and CMS-PDFT :cite:`Bao2020_CMS` are available by specifying ``method="XMS"`` and ``method="CMS"``, respectively, in place of ``method="lin"``.
For L-PDFT only, it is also possible to carry out calculations involving multiple different point group symmetries or spins using ``multi_state_mix``, similar to ``state_average_mix_``.



