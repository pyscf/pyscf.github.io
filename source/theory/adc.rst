.. _theory_adc:

************************************************
Algebraic diagrammatic construction (ADC) scheme
************************************************

*Modules*: :mod:`adc`

Background: Propagator Theory
=============================

The algebraic diagrammatic construction (ADC) theory has a close connection with the propagator theory, where the central object of interest is a propagator of the many-body system. 
The general expression for a retarded propagator that describes response of a many-electron system to an external perturbation with frequency :math:`\omega` is:

.. math::

        G_{\mu\nu}(\omega) &= G_{\mu\nu}^{+}(\omega) \pm G_{\mu\nu}^{-}(\omega)\\
               &=\langle\Psi_{0}^{N}|q_{\mu} (\omega - H + E_{0})^{-1}q_{\nu}^{\dagger}|\Psi_{0}^{N}\rangle\\
               &\pm \langle\Psi_{0}^{N}|q_{\nu}^{\dagger} (\omega + H - E_{0})^{-1}q_{\mu}|\Psi_{0}^{N}\rangle,

where :math:`G_{\mu\nu}^{+}(\omega)` and :math:`G_{\mu\nu}^{-}(\omega)` are the forward and backward components of :math:`G_{\mu\nu}(\omega)`, 
:math:`{H}` is the electronic Hamiltonian, :math:`|\Psi_{0}^{N}\rangle` and :math:`{E}_{0}^{N}` are the ground-state :math:`{N}`-electron wavefunction and energy, respectively. 
The frequency can be expressed as :math:`\omega \equiv \omega'` + :math:`i\eta`, where :math:`\omega'` is the real part of :math:`\omega` and :math:`i\eta` is an infinitesimal imaginary number. 
Depending on the form of the operators :math:`q_{\nu}^{\dagger}` and :math:`q_{\mu}`, the propagator :math:`G_{\mu\nu}(\omega)` can describe various spectroscopic processes. 
Choosing :math:`q_{\nu}^{\dagger} = a_{p}^{\dagger}a_{q} - \langle\Psi_{0}|a_{p}^{\dagger}a_{q}|\Psi_{0}\rangle`, where the operators :math:`a_{p}` and :math:`a_{p}^{\dagger}` 
are the usual one-electron annihilation and creation operators, corresponds to a polarization propagator that provides information about electronic excitations in optical (e.g., UV/vis)
spectroscopy. Alternatively, a propagator called the one-electron Green's Function with :math:`q_{\nu}^{\dagger} = a_{p}^{\dagger}` describes electron attachment (EA) and ionization (IP) processes.

An alternative way of expressing any propagtaor is its Lehmann (or spectral) representation. As an example, the Lehmann (or spectral) representation of the one-electron Green's Function is as follows:

.. math::
        G_{pq}(\omega) &= \sum_{n}\frac{\langle\Psi_{0}^{N}|a_{p}|\Psi_{n}^{N+1}\rangle \langle\Psi_{n}^{N+1}|a_{q}^{\dagger}|\Psi_{0}^{N}\rangle}{\omega - E_{n}^{N+1}+ E_{0}^{N}}\\
                       &+  \sum_{n}\frac{\langle\Psi_{0}^{N}|a_{q}^{\dagger}|\Psi_{n}^{N-1}\rangle \langle\Psi_{n}^{N-1}|a_{p}|\Psi_{0}^{N}\rangle}{\omega + E_{n}^{N-1} - E_{0}^{N}}

where the forward component corresponds to EA and the backward component to IP.
Here, the resolution of identity is carried out over all exact eigenstates :math:`|{\Psi_{n}^{N+1}}\rangle (|{\Psi_{n}^{N-1}}\rangle)` of the electron-attached (-ionized) system 
with energies :math:`E_{n}^{N+1} (E_{n}^{N-1})`. The poles of :math:`G_{pq}(\omega)` provide information about vertical electron attachment 
:math:`(\omega_{n} = E_{n}^{N+1} - E_{0}^{N})` and ionization :math:`(\omega_{n} = E_{0}^{N} - E_{n}^{N - 1})` energies, while the residues 
(e.g., :math:`\langle\Psi_{0}^{N}|a_{p}|\Psi_{n}^{N+1}\rangle \langle\Psi_{n}^{N+1}|a_{q}^{\dagger}|\Psi_{0}^{N}\rangle`) describe probabilities of the corresponding transitions.

For each component of the propagator, the Lehmann representation can be compactly written in a matrix form as:

.. math::
    \textbf{G}_{\pm}(\omega) = \mathbf{\tilde{X}}_{\pm}(\omega - \boldsymbol{\tilde{\Omega}}_{\pm})^{-1}\mathbf{\tilde{X}}_{\pm}^{\dagger},

where :math:`\boldsymbol{\tilde{\Omega}}_{\pm}` are diagonal matrices of the transition energies (:math:`\omega_{n}`) and :math:`\mathbf{\tilde{X}}_{\pm}` are matrices of the transition amplitudes. For example, for the one-electron Green's function, :math:`\boldsymbol{\tilde{\Omega}}_{\pm}` correspond to vertical attachment/ionization energies, while transition amplitudes are defined as :math:`\tilde{X}_{+pn} = \langle\Psi_{0}^{N}|a_{p}|\Psi_{n}^{N+1}\rangle` and :math:`\tilde{X}_{-qn} = \langle\Psi_{0}^{N}|a_{q}^{\dagger}|\Psi_{n}^{N-1}\rangle`.

Algebraic diagrammatic construction
===================================
The propagator :math:`\textbf{G}_{\pm}(\omega)` shown above is exact if expressed using the exact eigenstates of the system. 
However, computing the exact propagator is computationally very expensive and many approximations have been introduced to reduce the computational cost.
In the ADC theory, :math:`\textbf{G}_{\pm}(\omega)` is expressed in a basis of many-electron wavefunctions where the matrix representation of the Hamiltonian :math:`H` 
is no longer diagonal and each component of the propagator can be written as:

.. math::
    \textbf{G}_{\pm}(\omega) = \textbf{T}_{\pm}(\omega\textbf{S}_{\pm}-\textbf{M}_{\pm})^{-1}\textbf{T}_{\pm}^{\dagger}

Here, :math:`\textbf{M}_{\pm}` and :math:`\textbf{T}_{\pm}` are called the effective Hamiltonian and transition moments matrices, respectively. 
For generality, we assume that the new basis of electronic states is non-orthogonal with the overlap matrix :math:`\textbf{S}_{\pm}`. 

In the non-Dyson ADC methods, the propagators :math:`\textbf{G}_{\pm}(\omega)` are evaluated by expanding the matrices :math:`\textbf{M}`, :math:`\textbf{S}`, and :math:`\textbf{T}` 
in a perturbative series and truncating this expansion at the :math:`n`-th order. For example,

.. math::
    \textbf{M}_{\pm} &\approx \textbf{M}^{(0)}_{\pm} + \textbf{M}^{(1)}_{\pm} + ...+ \textbf{M}^{(n)}_{\pm}\\

Plugging the above equations into :math:`\textbf{G}_{\pm}(\omega)` leads to equations of the :math:`n`-th order ADC approximation (ADC(n)). 
The ADC(n) vertical transition energies can be computed by solving the generalized eigenvalue problem:

.. math::
    \textbf{M}_{\pm}\textbf{Y}_{\pm}=\textbf{S}_{\pm}\textbf{Y}_{\pm}\boldsymbol{\Omega}_{\pm}

where :math:`\boldsymbol\Omega_{\pm}` is a diagonal matrix of eigenvalues and :math:`\textbf{Y}_{\pm}` is the matrix of eigenvectors that can be used to compute spectroscopic amplitudes 

.. math::
    \textbf{X}_{\pm}=\textbf{T}_{\pm}\textbf{S}_{\pm}^{-\frac{1}{2}}\textbf{Y}_{\pm}

which provide information about transition intensities. 
Using :math:`\boldsymbol{\Omega}_{\pm}` and :math:`\mathbf{X}_{\pm}` allows to compute the ADC(n) propagator and density of states

.. math::
        \mathbf{G}_{\pm}(\omega) &= \mathbf{X}_{\pm} \left(\omega - \boldsymbol{\Omega}_{\pm}\right)^{-1}  \mathbf{X}_{\pm}^\dagger \\
        A(\omega) &= -\frac{1}{\pi} \mathrm{Im} \left[ \mathrm{Tr} \, \mathbf{G}_{\pm}(\omega) \right]

Methods
========
At the present moment, the ADC implementation in Pyscf can be used to perform calculations of electron affinities (EA) and ionization potentials (IP) for closed- and open-shell molecules starting with the RHF, ROHF, or UHF reference SCF wavefunctions. Three ADC approximations are currently available: ADC(2), ADC(3), as well as the extended second-order approximation ADC(2)-X. The capabilities of the ADC module at present are summarized in the following table:

========= ==============  ================  ===========
Method    Reference       Spin-adaptation   Properties
--------- --------------  ----------------  -----------
ADC(2)    RHF, ROHF, UHF     No                IP, EA
ADC(2)-X  RHF, ROHF, UHF     No                IP, EA
ADC(3)    RHF, ROHF, UHF     No                IP, EA
========= ==============  ================  ===========

When using results of this code for publications, please cite the following paper: "Third-order algebraic diagrammatic construction theory for electron attachment and ionization energies: Conventional and Green’s function implementation", S. Banerjee and A.Y. Sokolov, J. Chem. Phys. 151, 224112 (2019).
