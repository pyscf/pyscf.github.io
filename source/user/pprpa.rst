.. _user_pprpa:

particle-particle Random Phase Approximation (ppRPA)
************************************************************

*Modules*: :py:mod:`pyscf.pprpa`


Introduction
============

particle-particle Random Phase Approximation (ppRPA):cite:`vanaggelenExchangecorrelationEnergyPairing2013`,
which was originally used to describe the nuclear many-body correlation,
has been developed to predict ground-state and excited-state properties of molecular and bulk systems.
The ppRPA correlation energy is exact up to the second order electron-electron interaction and is equivalent to ladder coupled cluster doubles.
For calculations of excitation energies in ppRPA,
the excitation energies of the N-electron system can be calculated as the differences between the two-electron addition
energies of the (N-2)-electron system from the particle-particle channel.
Similarly,
the excitation energies can also be obtained from the differences between the two-electron removal energies of the (N+2)-electron system from the hole-hole channel.
The choice of the pp or the hh channels enhances the flexibility of the ppRPA method.

The ppRPA classes take molecule and supercell objects on the equal footing,
which can take ``RHF``, ``UHF``, ``RKS``, ``UKS`` from either a molecule or a supercell calculation.

Correlation Energy
==================
The evaluation of the ppRPA requires the full diagonalization of the ppRPA matrix,
which scales as :math:`N^6`.

Excitation Energy
=================

The formal scaling of ppRPA for computing excitation energies is :math:`N^4` with the Davidson algorithm:cite:`yangExcitationEnergiesParticleparticle2014`.
The computational cost can be further significantly reduced by using active-space approach:cite:`liLinearScalingCalculations2023`,
which directly truncates ppRPA matrix in the molecular orbital space without loss of accuracy.

In the current implementation,
the active-space approach can be used in both restricted and unrestricted calculation,
the Davidson algorithm is only support in the restricted calculation.

particle-particle channel
-------------------------
To get two-electron addition energies in the particle-particle channel,
set ``charge=2`` in the ``mol`` or ``cell`` object to create an (N-2)-electron system.

hole-hole channel
-----------------
To get two-electron removal energies in the particle-particle channel,
set ``charge=-2`` in the ``mol`` or ``cell`` object to create an (N+2)-electron system.
