.. _user_adc:

*****************************************
Algebraic diagrammatic construction (ADC)
*****************************************

*Modules*: :py:mod:`pyscf.adc`

Introduction
=============================
Algebraic diagrammatic construction theory (ADC) is a post-Hartree-Fock
approach for computing excited electronic states.
:cite:`Schirmer1983,Schirmer1998,Dreuw:2014p82,Banerjee:2023p3037`
This is done by approximating a linear-response function
(so-called propagator) to the :math:`n`-th order in perturbation theory,
defining the hierarchy of ADC(n) methods.
Depending on the property being investigated, propagators can be of different
types. Some common examples include the polarization propagator for neutral
electronic excitations, the one-particle Green's function for charged excitations,
and the two-electron propagator for Auger electron spectroscopy.
The different propagators lead to different variants of the ADC methods.

At present, the `adc` module in PySCF can be used to calculate neutral electronic
excitations from the polarization propagator and charged excitations
(ionization potentials, electron affinities) from the one-particle
Green's function.
For each property (EE, IP, and EA) three ADC approximations are available:
strict second-order (ADC(2)), extended second-order (ADC(2)-X), and
strict third-order (ADC(3)) methods. :cite:`Dreuw:2014p82,Banerjee:2023p3037`
Each method provides access to the excitation energies and transition intensities
via a 'one-shot' calculation.

A conventional (single-reference) ADC(n) calculation involves computing the
ground-state energy and wavefunction that correspond to those of the
:math:`n`-th order Møller--Plesset perturbation theory (MPn), followed by evaluating
the ADC(n) excitation energies and the corresponding transition probabilities.

Some examples of ADC calculations for the HF molecule are shown below::

        from pyscf import gto, scf, adc
        mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')
        mf = scf.RHF(mol).run()

        myadc = adc.ADC(mf)
        myadc.verbose = 4

        #IP-ADC(3) for 1 root
        myadc.method = "adc(3)"
        myadc.method_type = "ip"
        eip,vip,pip,xip = myadc.kernel()

        #EA-ADC(2)-x for 1 root
        myadc.method = "adc(2)-x"
        myadc.method_type = "ea"
        eea,vea,pea,xea = myadc.kernel()

        #EE-ADC(2) for 3 roots
        myadc.compute_properties = False
        myadc.method = "adc(2)"
        myadc.method_type = "ee"
        eee,vee,pee,xee = myadc.kernel(nroots = 3)

Here, ``method`` specifies the level of ADC approximation,
``method_type`` defines the type of excitation (EE, IP, and EA), and
``nroots`` defines the number of excited states to be computed (default is 1).

More examples can be found in
:source:`examples/adc`.


Spin-restricted and spin-unrestricted calculations
==========================================================================
The `adc` module can be used to compute EE's, IP's, and EA's
for closed- and open-shell molecules starting with the restricted (RHF) or
unrestricted (UHF) Hartree-Fock reference wavefunctions using the RADC or UADC
modules.
Additionally, open-shell ADC calculations can be carried out using the
restricted open-shell (ROHF) reference wavefunction, which are supported by
the UADC module.
See :ref:`user_scf` to learn more about the different reference wavefunctions.

Shown below is an example of the IP- and EA-UADC(2) calculations for the
open-shell OH radical::

        from pyscf import gto, scf, adc
        mol = gto.Mole()
        mol.atom = [
        ['O', ( 0., 0.    , 0.   )],
        ['H', ( 0., 0.    , 0.969)],]
        mol.basis = 'aug-cc-pvdz'
        mol.spin  = 1
        mol.build()

        mf = scf.UHF(mol).run()

        #IP-UADC calculation for 4 roots using the UHF reference
        myadc = adc.ADC(mf)
        myadc.verbose = 4

        myadc.method = "adc(2)"
        myadc.method_type = "ip"
        eip,vip,pip,xip = myadc.kernel(nroots = 4)

        mf = scf.ROHF(mol).run()

        #EA-UADC calculation for 4 roots using the ROHF reference
        myadc = adc.ADC(mf)
        myadc.verbose = 4

        myadc.method = "adc(2)"
        myadc.method_type = "ea"
        eea,vea,pea,xea = myadc.kernel(nroots = 4)


More examples can be found in
:source:`examples/adc`.


Excited-state properties
=========================
The `adc` module supports computing transition and excited-state properties.
For EE, oscillator strengths are available (currently, only in the UADC
implementation). :cite:`Dreuw:2014p82`
The IP- and EA-ADC methods provide spectroscopic factors, which measure the
probabilities of charged excitations in photoelectron spectra. :cite:`Banerjee:2023p3037`
Computation of transition properties is performed by default and can be switched
off by setting ``compute_properties = False`` ::

        myadc.compute_properties = False
        myadc.method = "adc(3)"
        myadc.method_type = "ip"
        myadc.kernel(nroots = 3)

Open-shell calculations using EE-UADC further support evaluating the spin
square operator expectation values for the excited states (<S^2>). 
The <S^2> values are not computed by default, they can be requested using
the ``compute_spin_square = True`` flag. See a relevant example
for more details: :source:`examples/adc/08-open_shell_spin_square.py`.

For the IP- and EA-ADC methods, the `adc` module can be used to compute
the Dyson orbitals :cite:`Oana2007` visualizing the wavefunction of
ionized hole or attached electron for states with large
spectroscopic factors (> 0.5). ::

        dyson_orb = myadc.compute_dyson_mo()

Additionally, the ADC codes enable computing the one-particle reduced
density matrices (1-RDMs) for the ground and excited electronic states.
To obtain the ground-state 1-RDM, run the ``make_ref_rdm1()`` function
after a successful ADC calculation::

        myadc.kernel(nroots = 3)
        rdm1_ref = myadc.make_ref_rdm1()

The excited-state 1-RDMs can be calculated as follows::

        rdm1_exc = myadc.make_rdm1()

The 1-RDMs can be used to compute excited-state one-particle properties,
such as dipole moments, see the following example:
:source:`examples/adc/07-closed_shell_1RDMS.py`.
The 1-RDM functionality is currently limited, see below.

Transition 1-RDMs between the ground and excited states are also available
for all methods, except EE-RADC.
They are returned as spectroscopic amplitudes at the end of a successful ADC
calculation with ``compute_properties = True``, for example::

        from pyscf import gto, scf, adc
        mol = gto.Mole()
        mol.atom = [
        ['O', ( 0., 0.    , 0.   )],
        ['H', ( 0., 0.    , 0.969)],]
        mol.basis = 'aug-cc-pvdz'
        mol.spin  = 1
        mol.build()

        mf = scf.UHF(mol).run()

        # Transition 1-RDMs (alpha and beta spin) computed using EE-UADC(2)
        myadc = adc.ADC(mf)
        myadc.verbose = 4
        myadc.method = "adc(2)"
        myadc.method_type = "ee"
        trdm_a, trdm_b = myadc.kernel(nroots = 4)[3]


Analysis of excited-state properties
=====================================
The `adc` module allows to analyze the ADC(n) eigenvectors
for characterizing the nature of electronic transitions. When
``compute_properties = True`` is set, this analysis will also display the largest
contributions to the spectroscopic factors (IP, EA) or amplitudes (EE) for singly
excited states. The analysis of ADC(n) eigenvectors
and spectroscopic factors can be invoked using the ``analyze()`` function::

        myadc.kernel(nroots = 3)
        myadc.analyze()


Core excitation energies and spectra
=====================================
The IP-ADC code supports calculating core ionization energies and
X-ray photoelectron spectra (XPS) using the core-valence separation technique
(CVS). To invoke CVS, specify the ``ncvs`` parameter 
of IP-ADC object before running the ``kernel()`` function.
The ``ncvs`` parameter should be set to an integer defining the index of the core
orbital that is expected to be ionized first in the XPS spectrum.
For example, probing the 1s orbital of C in CO can be done by setting
``ncvs = 2`` since the C 1s orbitals are higher in energy than O 1s.


Algorithms and job control
===========================

The current capabilities of the `adc` module are summarized in the
following table:

============= ================ =====================================================
 Method        Reference        Available properties
------------- ---------------- -----------------------------------------------------
 EE-ADC(2)     RHF, UHF, ROHF   EEs, osc. strengths (UADC), spin square (UADC)
 EE-ADC(2)-X   RHF, UHF, ROHF   EEs, osc. strengths (UADC), spin square (UADC)
 EE-ADC(3)     RHF, UHF, ROHF   EEs, osc. strengths (UADC), spin square (UADC)
 IP-ADC(2)     RHF, UHF, ROHF   IPs, core IPs, spec. factors, Dyson orb
 IP-ADC(2)-X   RHF, UHF, ROHF   IPs, core IPs, spec. factors, Dyson orb
 IP-ADC(3)     RHF, UHF, ROHF   IPs, core IPs, spec. factors, Dyson orb
 EA-ADC(2)     RHF, UHF, ROHF   EAs, spec. factors, Dyson orb
 EA-ADC(2)-X   RHF, UHF, ROHF   EAs, spec. factors, Dyson orb
 EA-ADC(3)     RHF, UHF, ROHF   EAs, spec. factors, Dyson orb
============= ================ =====================================================

The ADC(n) calculations are performed using one of the three algorithms for
handling the two-electron integrals:

* In-core

  All tensors such as two-electron integrals and
  amplitudes are stored in memory. This is the default algorithm used when
  sufficient memory is available.


* Out-of-core

  Uses disk to store expensive tensors.
  This algorithm is invoked when storing two-electron integrals requires
  more memory than specified by the ``max_memory`` keyword (in MB).
  See :source:`examples/adc/05-outcore.py`


* Density fitting (DF)

  The memory and disk usage can be greatly reduced by approximating the
  two-electron integrals with density fitting. An example of the ADC(2)
  calculation with density fitting is provided below::

     from pyscf import gto, scf, adc, df
     mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')

     mf = scf.RHF(mol).density_fit('ccpvdz-jkfit').run()
     myadc = adc.ADC(mf).density_fit('ccpvdz-ri')
     eip,vip,pip,xip = myadc.kernel()

  Density fitting introduces small errors in excitation energies
  ($\sim 10^{-3}$ eV), provided that the appropriate auxiliary basis
  set is used. :cite:`Banerjee2021`
  For calculations with more than 300 orbitals, using density fitting
  is strongly recommended.

  More examples can be found in: :source:`examples/adc/06-dfadc.py`.


Current limitations
====================

Some limitations of current implementation are listed below:

* The EE-RADC code does not support calculating oscillator strengths.
  This property can be computed using the EE-UADC code (i.e., by using
  the UHF reference) for a closed- or open-shell molecule.

* The EE- and IP/EA-RADC codes compute only the states of lowest spin:
  S = 0 (singlet) and S = 1/2 (doublet), respectively. Using the
  corresponding UADC code allows to compute excitations with $\Delta$(S) = 0,
  $\pm$1 for EE and $\Delta$(S) = $\pm$1/2, $\pm$3/2 for IP and EA.

* Computing spin square expectation values is currently only available for
  EE-UADC.

* The reference and excited-state 1-RDMs are not implemented for EE-RADC.
  Also, the reference 1-RDMs are not available for any UADC method.

* The EE-UADC(3) calculations of excited-state one-particle reduced density
  matrices include correlation contributions up to EE-UADC(2)-X, i.e. the
  third-order terms are missing from the singles-singles and singles-doubles
  coupling sectors.

* The EE-UADC(3) oscillator strengths do not include the contributions from
  third-order amplitudes.

