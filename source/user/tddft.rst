.. _theory_tdscf:

Time-dependent Hartree-Fock and Density Functional Theory
*********************************************************

*Modules*: :mod:`tdscf`, :mod:`pbc.tdscf`

Introduction
============
PySCF implements the time-dependent Hartree-Fock (TDHF) and 
time-dependent density functional theory (TDDFT) (frequency domain)
linear response theories to compute excited-state energies 
and transition properties in the :mod:`tdscf` module. 
A minimal example that runs a TDDFT calculation is as follows ::

    from pyscf import gto, scf, dft, tddft
    mol = gto.Mole()
    mol.build(
        atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
        basis = '631g',
        symmetry = True,
        verbose = 4,
    )

    mf = dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.kernel()

    mytd = tddft.TDDFT(mf)
    mytd.nstates = 10
    mytd.kernel()
    mytd.analyze()

The example above computes the excitation energies, oscillator strengths
and transition dipole moments of the ten lowest singlet exicted states.

Theory
======
Using first-order time-dependent perturbation theory within HF or KS theory,
one obtains the non-Hermitian TDHF or TDDFT equations for the excitation energies
:cite:`DreHea2005`:

.. math::

    \left(\!\!\begin{array}{ll} 
        \mathbf{A} & \mathbf{B} \\
        \mathbf{B}^\ast & \mathbf{A}^\ast
    \end{array}\!\!\right)
    \left(\!\!\begin{array}{c}
        \mathbf{X} \\ \mathbf{Y}
    \end{array}\!\!\right) = \omega 
    \left(\!\!\begin{array}{lr}
        \mathbf{1} & \mathbf{0} \\
        \mathbf{0} & -\mathbf{1}
    \end{array}\!\!\right)
    \left(\!\!\begin{array}{c}
        \mathbf{X} \\ \mathbf{Y}
    \end{array}\!\!\right) \;,
    
where :math:`\mathbf{A}` and :math:`\mathbf{B}` are the orbital
hessians which also appear in the stability analysis for reference states (see :numref:`stability_analysis`),
:math:`\omega` is the excitation energy, and 
:math:`\mathbf{X}` and :math:`\mathbf{Y}` represent the response of the density matrix.
In cases where the system possesses
a degenerate ground state or has triplet instabilities, the algorithms used to solve the 
above equations may be unstable. This can be solved by applying the Tamm-Dancoff approximation (TDA) :cite:`HirHea1999`,
which simply neglects the :math:`\mathbf{B}` and :math:`\mathbf{Y}` matrices and leads to a 
Hermitian eigenvalue problem

.. math::

    \mathbf{AX} = \omega\mathbf{X} \;.

.. Note that when TDA is applied, TDHF reduces to the same formalism as configuration interaction singles (CIS),
.. because the Coulomb operator projected onto the space of singly-excited Slater determinants
.. yields terms that are equivalent to the linear response of the HF exchange and Coulomb potentials.
.. However, in the case of TDDFT, linear response of the exchange-correlation (XC) potential leads to 2nd order derivative of 
.. the XC functional, which does not appear in the ground-state DFT.


Methods
=======
For TDHF or TDDFT calculations,
the reference state can be either restricted or unrestricted::

    mytd = mol.RKS().run().TDDFT().run()
    mytd = mol.UKS().run().TDDFT().run()

By default, only singlet excited states are computed. 
In order to compute triplet excited states, one needs to set the 
attribute :attr:`.singlet` to ``False``::

    mytd.singlet = False
    mytd.kernel()

One can also perform symmetry analysis by calling the :func:`.analyze()` method,
which also computes the oscillator strengths and dipole moments::

    mytd.analyze(verbose=4)


Property calculation
====================

Oscillator strengths
--------------------
Oscillator strengths for each excited state can be computed in 
both length and velocity gauges::

    mytd.oscillator_strength(gauge='length')
    mytd.oscillator_strength(gauge='velocity')

Higher order corrections :cite:`LesEgiLi2015` 
to the oscillator strength can also be included::

    #include corrections due to magnetic dipole and electric quadruple
    mytd.oscillator_strength(gauge='velocity', order=1)
    #also include corrections due to magnetic quadruple and electric octupole
    mytd.oscillator_strength(gauge='velocity', order=2)

Transition moments
------------------
PySCF implements various types of transition moments between the reference SCF state and 
the TDHF or TDDFT excited states. These include:

* electric dipole, quadrupole and octupole transition moments in both length and velocity gauges::

    mytd.transition_dipole()
    mytd.transition_velocity_dipole()
    mytd.transition_quadrupole()
    mytd.transition_velocity_quadrupole()
    mytd.transition_octupole()
    mytd.transition_velocity_octupole()

* magnetic dipole and quadrupole transition moments::

    mytd.transition_magnetic_dipole()
    mytd.transition_magnetic_quadrupole()

Nuclear gradients
-----------------
Analytic nuclear gradients are available for TDHF and TDDFT, 
and they can be computed as follows::

    tdg = mytd.Gradients()
    g1 = tdg.kernel() #default will compute the gradients of first excited state
    g1 = tdg.kernel(state=1) #first excited state
    g2 = tdg.kernel(state=2) #second excited state

Natural transition orbital analysis
-----------------------------------
Natural transition orbitals (NTOs) can be computed by 
singular value decomposition of the transition density matrix.
In PySCF, these orbitals can be obtained as follows::

    weights, nto_coeff = mytd.get_nto(state=1)

where ``nto_coeff`` are the coefficients for NTOs represented in AO basis,
and they are ordered as occupied orbitals followed by virtual orbitals.

References
==========
.. bibliography:: ref_tddft.bib
   :style: unsrt
