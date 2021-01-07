.. _theory_scf:

***********************************
Self-consistent field (SCF) methods
***********************************

*Modules*: :mod:`scf`, :mod:`pbc.scf`, :mod:`soscf`

Introduction
============
The SCF methods include both Hartree-Fock (HF) theory and Kohn-Sham (KS) density functional theory (DFT).
This Chapter summarizes the general SCF capabilities in PySCF.
More details specific to DFT are introduced in :numref:`theory_dft`.

A minimal example of using the :mod:`scf` module is as follows::

    from pyscf import gto, scf
    mol = gto.M(
        atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
        basis = 'ccpvdz',
        symmetry = True,
    )
    mf = scf.HF(mol)
    mf.kernel()

This will run a HF calculation with the default SCF settings.


Theory
======
In HF and KS-DFT, the ground-state wavefunction is expressed as a single Slater determinant 
of molecular orbitals (MOs) :math:`\psi`,

.. math::

   |\Phi_0\rangle = \frac{1}{\sqrt{N!}}
   \begin{vmatrix} 
   \psi_1(\mathbf{r}_1) &\psi_2(\mathbf{r}_1) &\dots  &\psi_N(\mathbf{r}_1)\\
   \psi_1(\mathbf{r}_2) &\psi_2(\mathbf{r}_2) &\dots  &\psi_N(\mathbf{r}_2)\\
   \vdots               &\vdots               &\ddots &\vdots\\
   \psi_1(\mathbf{r}_N) &\psi_2(\mathbf{r}_N) &\dots  &\psi_N(\mathbf{r}_N)
   \end{vmatrix} \;.

In this way, the electron correlation is treated equivalently as a mean-field Coulomb repulsion.
The total electronic energy :math:`E=\langle\Psi_0|\hat{H}|\Psi_0\rangle` 
is then minimized by solving the Hartree-Fock equations

.. math::

    \mathbf{FC} = \mathbf{SC\varepsilon} \;

with :math:`\mathbf{F}`, :math:`\mathbf{C}`, and :math:`\mathbf{S}` 
denoting the Fock matrix, 
the expansion coefficients of MOs,
and the overlaps among basis functions, respectively. :cite:`Pop1954`

Methods
=======
Based on different treatments to spin symmetries, PySCF implements 
four variants of HF and KS-DFT methods.

* Restricted (RHF/RKS)

    The spatial parts of alpha and beta spin orbitals are constrained to be identical.
    This is appropriate for systems with closed-shell ground states.

* Unrestricted (UHF/UKS)

    The alpha and beta spin orbitals are computed independently by
    solving two sets of Fock equations, with the possibility of introducing spin contaminations.
    This is appropriate for most open-shell systems.

* Restriced open-shell (ROHF/ROKS)

    Like RHF (RKS), only one set of orbitals is computed, where the open-shell orbitals are effectively 
    fractionally occupied. :cite:`Roo1960` 
    The final wavefunction is an eigenfunction of the :math:`\hat{S}^2` operatoer,
    thus no spin contamination is introduced.
    This is appropriate for open-shell systems when unrestricted calculations give large spin contaminations,
    *e.g.*, in the case of bond stretching.

* Generalized (GHF/GKS)

    In the methods above, the spin orbitals are treated as real orbitals and are eigenfunctions of the 
    :math:`\hat{S}_z` operator. More generally, the spin orbitals can be expressed as :cite:`SeePop1977`

    .. math::
        \psi(r,\xi) = \psi_\alpha(r)\alpha(\xi) + \psi_\beta(r)\beta(\xi) \;, 

    where :math:`\psi_\alpha(r)` and :math:`\psi_\beta(r)` are complex spatial orbitals, and 
    :math:`\alpha(\xi)` and :math:`\beta(\xi)` are eigenfunctions of the :math:`\hat{S}_z` operator.
    This is useful when the previous methods do not provide stable solutions 
    (see :download:`examples/scf/17-stability.py </../examples/scf/17-stability.py>`)
    or when spin-orbit coupling is considered 
    (see :download:`examples/scf/44-soc_ecp.py </../examples/scf/44-soc_ecp.py>`).

Calculations with these methods can be invoked by creating an instance of the corresponding class::

    mf = scf.RHF(mol)
    mf = scf.UHF(mol)
    mf = scf.ROHF(mol)
    mf = scf.GHF(mol)
    mf.kernel()

More examples can be found in
:download:`examples/scf/00-simple_hf.py </../examples/scf/00-simple_hf.py>`,
:download:`examples/scf/01-h2o.py </../examples/scf/01-h2o.py>`,
:download:`examples/scf/02-rohf_uhf.py </../examples/scf/02-rohf_uhf.py>`,
:download:`examples/scf/02-ghf.py </../examples/scf/02-ghf.py>`.

Controllable parameters
=======================

Initial guess
-------------
PySCF provides several options as the initial guesses for solving the 
SCF problem. These can be specified by setting the attribute :attr:`.init_guess`
to the following values:

* ``'minao'`` (default)

    Superpostion of atomic density projected from the atomic natural orbital (ANO) basis.

* ``'1e'``

    The core Hamiltonian is diagonalized to get the initial MOs. 

* ``'atom'``

    Superposition of atomic HF density matrix.

* ``'huckel'``

    A Hückel guess based on on-the-fly atomic HF calculations. :cite:`Leh2019`

* ``'vsap'``

    Superposition of atomic potentials. Note this is only available for DFT calculations. 
    
* ``'chk'``

    Read the existing SCF results from the checkpoint file as the initial guess.

Alternatively, the user could manually set the initial guess density matrix for an SCF calculation 
by setting the ``dm0`` argument. 
For example, the followings script first computes the HF density matrix for :math:`\rm Cr^{6+}` cation,  
which is then used as the initial guess for the HF calculation of :math:`\rm Cr` atom. ::

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
:download:`examples/scf/15-initial_guess.py </../examples/scf/15-initial_guess.py>`,
:download:`examples/scf/31-cr_atom_rohf_tune_init_guess.py </../examples/scf/31-cr_atom_rohf_tune_init_guess.py>`.

Converging SCF iterations
-------------------------
PySCF implements two types of algorithms to converge the SCF iterations, namely,
direct Inversion in the iterative subspace (DIIS) and second-order SCF (SOSCF).

* DIIS (default)

    With DIIS, the Fock matrix at each iteration is extrapolated using the Fock matrices from the previous iterations,
    by minimizing the norm of the commutator :math:`[\mathbf{F},\mathbf{PS}]`. :cite:`Pul1980,Pul1982`
    Two variants of DIIS are also implemented in PySCF, namely, EDIIS :cite:`KudScuCan2002` 
    and ADIIS :cite:`HuYan2010`, where the objective functions to be minimized  
    are expressed as energy funtions. 
    Examples of selecting different DIIS schemes can be found in
    :download:`examples/scf/24-tune_diis.py </../examples/scf/24-tune_diis.py>`.

* SOSCF
    To achieve quadratic convergence for orbital optimizations, 
    PySCF implements a general second-order solver called the
    Co-iterative augmented hessian (CIAH) method. :cite:`Sun2016,Sun2017`
    This can be invoked by decorating the SCF objects with the :func:`.newton` method::

        mf = scf.RHF(mol).newton()

    More examples can be found in 
    :download:`examples/scf/22-newton.py </../examples/scf/22-newton.py>`.

* Damping

    Damping of the Fock matrix can be applied before the DIIS starts.
    This is invoked by setting the attributes :attr:`.damp` and :attr:`.diis_start_cycle`.
    For example, ::

        mf.damp = 0.5
        mf.diis_start_cycle = 2

    implies that the DIIS will start at the second cycle, 
    and that the Fock matrix is dampped at the first cycle.

* Level shifting

    Applying level shift can help converge SCF for small gap systems.
    This is invoked by setting the attribute :attr:`.level_shift`.
    See examples in 
    :download:`examples/scf/03-level_shift.py </../examples/scf/03-level_shift.py>`,
    :download:`examples/scf/52-dynamically_control_level_shift.py </../examples/scf/52-dynamically_control_level_shift.py>`.

* Fractional occupation

    Fractional occupation can be invoked to converge SCF for small gap systems.
    See examples in
    :download:`examples/scf/54-fractional_occupancy.py </../examples/scf/54-fractional_occupancy.py>`.

Stability analysis
==================
PySCF allows detection of both internal and external instabilities 
for a given SCF calculation. See examples in 
:download:`examples/scf/17-stability.py </../examples/scf/17-stability.py>`.

Properties calculation
======================
Various properties can be computed by calling the corresponding functions,
*e.g.*, dipole moment (:func:`.dip_moment`), nuclear gradients (:func:`.nuc_grad_method`),
and Mülliken population (:func:`.mulliken_pop`), *etc.*

References
==========
.. bibliography:: ref_scf.bib
   :style: unsrt
