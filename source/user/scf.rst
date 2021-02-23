.. _theory_scf:

***********************************
Self-consistent field (SCF) methods
***********************************

*Modules*: :mod:`scf`, :mod:`pbc.scf`, :mod:`soscf`

Introduction
============
The SCF methods include both Hartree-Fock (HF) theory and Kohn-Sham (KS) density functional theory (DFT). This chapter summarizes
the general SCF capabilities in PySCF.
Details specific to DFT can be found in :numref:`theory_dft`.

A minimal example showing how to use the :mod:`scf` module is::

    from pyscf import gto, scf
    mol = gto.M(
        atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
        basis = 'ccpvdz',
        symmetry = True,
    )
    mf = scf.HF(mol)
    mf.kernel()

This will run a HF calculation with the default SCF settings.


Background
==========
In HF and KS-DFT, the ground-state wavefunction is expressed as a single Slater determinant :math:`\Phi_0` of molecular orbitals (MOs) :math:`\psi`,
:math:`\Phi_0 = \mathcal{A}|\psi_1(1)\psi_2(2) \ldots \psi_N(N)|`:

.. math::

   |\Phi_0\rangle = \frac{1}{\sqrt{N!}}
   \begin{vmatrix} 
   \psi_1(1) &\psi_2(1) &\dots  &\psi_N(1)\\
   \psi_1(2) &\psi_2(2) &\dots  &\psi_N(2)\\
   \vdots               &\vdots               &\ddots &\vdots\\
   \psi_1(N) &\psi_2(N) &\dots  &\psi_N(N)
   \end{vmatrix} \;.

The total electronic energy :math:`E=\langle\Psi_0|\hat{H}|\Psi_0\rangle` 
is then minimized, subject to orbital orthogonality. This is equivalent
to treating the electron interaction via a mean-field approximation.

The minimum is obtained by solving the Fock equation

.. math::

    \hat{f} \psi_i = \varepsilon_i \psi_i

with the Fock operator :math:`\hat{f}` defined as

.. math::

    \hat{f} = -\frac{1}{2}\nabla^2 + \hat{v}_{\rm ext} + \hat{J} + \hat{K} \;,

where :math:`-\frac{1}{2}\nabla^2` is the noninteracting kinetic energy operator,
:math:`\hat{v}_{\rm ext}` is the external potential,
:math:`\hat{J}` is the Coulomb operator, and
:math:`\hat{K}` is the exact exchange operator. :cite:`SzaOst2012`

Variants
========
The general spin-orbital can be written as
    .. math::
        \psi_i(1) = \phi_{i\alpha}(r)\alpha + \phi_{i\beta}(r)\beta \;, 

PySCF supplies four variants of the HF and KS-DFT methods  implementing different
constraints on :math:`\psi(1)`. 

* Restricted (RHF/RKS)

  Orbitals either have alpha or beta spin :math:`\psi_i =\phi_i(r)\alpha` or `\psi_i = \phi_i(r)\beta`, and
  for every orbital of alpha spin, there is one of beta spin with the same spatial component. The
  closed shell determinant is thus :math:`\Phi=\mathcal{A}|\phi_1(r_1)\alpha \phi_1(r_2)\beta \ldots \phi_{N/2}(r_{N-1})\alpha \phi_{N/2}(r_N)\beta|`
  and :math:`S=0`.

* Unrestricted (UHF/UKS)
  
  Orbitals either have alpha or beta spin, and alpha and beta spin orbitals can have different spatial components. The determinant is
  thus :math:`\Phi=\mathcal{A}|\phi_1(r_1)\sigma_1 \phi_2(r_2)\sigma_2 \ldots \phi_{N}(r_N)\sigma_N|` where :math:`\sigma \in \{\alpha,\beta\}`.
  This introduces spin contamination for states not of maximal :math:`S_z`.

* Restricted open-shell (ROHF/ROKS)
  For :math:`N_\alpha > N_\beta`. The first :math:`N_\beta` orbitals have the same spatial components
  for :math:`\alpha` and :math:`\beta` spin, while the remaining orbitals are of :math:`\alpha` spin.
  i.e. :math:`\Phi=\mathcal{A}|\phi_1 \alpha \phi_1\beta \ldots \phi_{N_\beta} \alpha \phi_{N_\beta}\beta \phi_{N_\beta+1}\alpha \ldots \phi_{N}\alpha|`
  The final wavefunction is an eigenfunction of the :math:`\hat{S}^2` operator with :math:`S_z=S`.

* Generalized (GHF/GKS)

  The general form of the spin-orbital :math:`\psi` is used.
  This is useful when the previous methods do not provide stable solutions 
  (see :source:`examples/scf/17-stability.py`)
  or when the Hamiltonian does not commute with :math:`\hat{S}_z` (e.g. with spin-orbit coupling)
  (see :source:`examples/scf/44-soc_ecp.py`).
  
Calculations with these methods can be invoked by creating an instance of the corresponding class::

    mf = scf.RHF(mol).run()
    mf = scf.UHF(mol).run()
    mf = scf.ROHF(mol).run()
    mf = scf.GHF(mol).run()
    mf = scf.RKS(mol).run()
    mf = scf.UKS(mol).run()
    mf = scf.ROKS(mol).run()
    mf = scf.GKS(mol).run()

More examples can be found in
:source:`examples/scf/00-simple_hf.py`,
:source:`examples/scf/01-h2o.py`,
:source:`examples/scf/02-rohf_uhf.py`, and
:source:`examples/scf/02-ghf.py`.

Controllable parameters
=======================

Initial guess
-------------
PySCF provides several options for the initial guess to solve the
SCF problem; see :cite:`Leh2019` for a review and assessment of
initial guesses. These can be specified by setting the attribute
:attr:`.init_guess` to the following values:

* ``'minao'`` (default)

    Superposition of atomic densities projected from the atomic natural orbital (ANO) basis.

* ``'1e'``

    The core Hamiltonian is diagonalized to get the initial MOs. For quantum chemistry problems,
    the use of the 1e guess is not recommended, because the guess is very bad.

* ``'atom'``

    Superposition of atomic HF density matrices. The atomic HF calculations are spin-restricted and employ spherically averaged occupations
    with ground states determined in :cite:`Leh2020`.

* ``'huckel'``

    A Hückel guess based on on-the-fly atomic HF calculations like in ``'atom'``. :cite:`Leh2019`

* ``'vsap'``

    Superposition of atomic potentials. Note this is only available for DFT calculations. :cite:`Leh2019`
    
* ``'chk'``

    Read the existing SCF results from the checkpoint file as the initial guess.

Alternatively, the user can manually set the initial guess density matrix for an SCF calculation 
through the ``dm0`` argument. 
For example, the following script first computes the HF density matrix for the :math:`\rm Cr^{6+}` cation,  
which is then used as an initial guess for a HF calculation of the :math:`\rm Cr` atom. ::

    #
    # use cation to produce initial guess
    #
    mol = gto.Mole()
    mol.build(
        symmetry = 'D2h',
        atom = [['Cr',(0, 0, 0)], ],
        basis = 'cc-pvdz',
        charge = 6,
        spin = 0,
    )

    mf = scf.RHF(mol)
    mf.kernel()
    dm1 = mf.make_rdm1()

    mol.charge = 0
    mol.spin = 6
    mol.build(False,False)

    mf = scf.RHF(mol)
    mf.kernel(dm0=dm1)

More examples can be found in 
:source:`examples/scf/15-initial_guess.py`, and
:source:`examples/scf/31-cr_atom_rohf_tune_init_guess.py`.

Converging SCF iterations
-------------------------
PySCF implements two algorithms to converge the SCF iteration, namely,
direct inversion in the iterative subspace (DIIS) and second-order SCF (SOSCF).

* DIIS (default)

    With DIIS, the Fock matrix at each iteration is extrapolated using  Fock matrices from the previous iterations,
    by minimizing the norm of the commutator :math:`[\mathbf{F},\mathbf{PS}]`. :cite:`Pul1980,Pul1982`
    Two variants of DIIS are  implemented in PySCF, namely, EDIIS :cite:`KudScuCan2002` 
    and ADIIS :cite:`HuYan2010`.
    Examples of selecting different DIIS schemes can be found in
    :source:`examples/scf/24-tune_diis.py`.

* SOSCF

    To achieve quadratic convergence for orbital optimizations, 
    PySCF implements a general second-order solver called the
    co-iterative augmented hessian (CIAH) method. :cite:`Sun2016,Sun2017`
    This can be invoked by decorating the SCF objects with the :func:`.newton` method::

        mf = scf.RHF(mol).newton()

    More examples can be found in 
    :source:`examples/scf/22-newton.py`.

* Damping

    Damping of the Fock matrix can be applied before DIIS starts.
    This is invoked by setting attributes :attr:`.damp` and :attr:`.diis_start_cycle`.
    For example, ::

        mf.damp = 0.5
        mf.diis_start_cycle = 2

    implies that DIIS will start at the second cycle, 
    and that the Fock matrix is damped at the first cycle.

* Level shifting

    A level shift forces a gap between the occupied and virtual Fock eigenvalues.
    Applying a level shift can help to converge SCF for small gap systems.
    This is invoked by setting the attribute :attr:`.level_shift`.
    See examples in 
    :source:`examples/scf/03-level_shift.py`, and
    :source:`examples/scf/52-dynamically_control_level_shift.py`.

* Fractional occupation

    Fractional occupation can be invoked to converge SCF for small gap systems.
    See the example in
    :source:`examples/scf/54-fractional_occupancy.py`.

* Smearing

    Smearing sets fractional occupancies according to a temperature function. See the example
    :source:`examples/pbc/23-smearing.py`.
	    
.. _stability_analysis:

Stability analysis
==================
PySCF allows detection of both internal and external instabilities 
for a given SCF calculation. See examples in 
:source:`examples/scf/17-stability.py`.

Property calculation
====================
Various properties can be computed by calling the corresponding functions,
for example, 

* dipole moment::
 
    mf.dip_moment()

* Mulliken population:: 

    mf.mulliken_pop()

* nuclear gradients::

    g = mf.Gradients()
    g.kernel()

    
References
==========
.. bibliography:: ref_scf.bib
   :style: unsrt
