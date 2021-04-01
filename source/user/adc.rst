.. _user_adc:

************************************************
Algebraic diagrammatic construction (ADC) scheme
************************************************

*Modules*: :mod:`adc`

Introduction
=============================
The algebraic diagrammatic construction theory (ADC(n)) is a post-Hartree-Fock
method used for computing correlated excited states of molecules.
:cite:`Schirmer1983,Schirmer1998`
The ADC methods involve a perturbative expansion of a propagator followed
by truncation at a certain order :math:`n` that defines the ADC(n) approximation.
Depending upon the property being investigated, propagators can be of different
types. Some common examples include the polarization propagator for neutral
electronic excitations, the one-particle Green's function for charged excitations,
and the two-electron propagator for Auger electron spectroscopy.
The different propagators lead to different variants of the ADC method.

At present, the `adc` module in PySCF can be used to calculate the one-particle
Green's function, that provides access to charged excitations,
such as ionization potentials (IP-ADC) and
electron affinities (EA-ADC) with different ADC(n) approximations
(ADC(2), ADC(2)-X and ADC(3)). :cite:`Trofimov2005,Banerjee2019`
The ADC(n) methods provide access to the
excitation energies and their corresponding transition intensities via a
'one-shot' calculation.

A simple ADC(n) computation involves the calculation of the ground-state energy
and wavefunction that correspond to those of the :math:`n`-th order
Møller--Plesset perturbation theory (MPn), followed by the evaluation of the
ADC(n) excitation energies (e.g., IP or EA) and the corresponding transition
probabilities.

An example of IP and EA computation using the default ADC method (ADC(2))
is shown below::

        from pyscf import gto, scf, adc
        mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')
        mf = scf.RHF(mol).run()
        myadc = adc.ADC(mf)
        myadc.kernel_gs()
        e_ip,v_ip,p_ip,x_ip,adc_es = myadc.ip_adc()
        e_ea,v_ea,p_ea,x_ea,adc_ea = myadc.ea_adc()

In the example shown above, the ground state calculation is performed by
running the function ``kernel_gs()``, whereas the excited state calculations
are carried out using the ``ip_adc()`` and ``ea_adc()`` functions.

Alternatively, both the ground and excited state calculations can be run
together using the ``kernel()`` function::

        from pyscf import gto, scf, adc
        mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')
        mf = scf.RHF(mol).run()
        myadc = adc.ADC(mf)
        myadc.kernel()

By default, the ``kernel()`` function performs a IP-ADC(2) calculation. One can specify the type of charged
excitation and order of the desired ADC computation::

        myadc.method = "adc(3)"
        myadc.method_type = "ea"
        myadc.kernel()

The ADC functions return the ``nroots`` lowest-energy eigenvalues. The
default value of ``nroots`` is set to 1. More roots can be requested using::

        myadc.kernel(nroots=3)

More examples can be found in
:source:`examples/adc/01-closed_shell.py`,
:source:`examples/adc/03-closed_shell_different_setup.py`.


Spin-restricted and spin-unrestricted calculations
==========================================================================
The `adc` module can be used to perform calculations of IP's and EA's of closed- and
open-shell molecules starting with the RHF and UHF reference
wavefunctions, leading to the RADC(n) and UADC(n) methods, respectively.
:cite:`Banerjee2021`
See :ref:`user_scf` to learn more about the different reference wavefunctions.

Shown below is an example of the IP- and EA-UADC(2) calculations for the
open-shell OH radical::

        from pyscf import gto, scf, adc
        mol = gto.M(atom='O 0 0 0; H 0 0 1', basis='aug-cc-pvdz')
        mol.spin  = 1
        mf = scf.UHF(mol).run()

        myadc = adc.ADC(mf)  # This is a UADC calculation
        myadc.kernel_gs()
        eip,vip,pip,xip,ip_es = myadc.ip_adc()
        eea,vea,pea,xea,ea_es = myadc.ea_adc()

More examples can be found in
:source:`examples/adc/ 02-open_shell.py`,
:source:`examples/adc/`04-open_shell_different_setup.py`.


Spectroscopic properties
=========================
The `adc` module supports calculation of the spectroscopic factors, which provide
information about probabilities of transitions in the photoelectron spectra. :cite:`Banerjee2021`
Computation of spectroscopic factors is performed by default and can be switched
off by setting ``compute_properties = False`` ::

        myadc.compute_properties = False
        myadc.method = "adc(3)"
        myadc.method_type = "ip"
        myadc.kernel(nroots = 3)

After the ADC calculation is performed, the `adc` module can be used to compute
the Dyson orbitals :cite:`Oana2007` corresponding to ionized and electron-attached states::

        dyson_orb = myadc.compute_dyson_mo()


Analysis of spectroscopic properties
=====================================
The `adc` module allows to perform the analysis of the ADC(n) eigenvectors, that
can be useful for characterizing the nature of electronic transitions. When
``compute_properties`` is set to True, this analysis will also display the largest
contributions to the spectroscopic factors. The analysis of the ADC(n) eigenvectors
and spectroscopic factors can be invoked using the ``analyze()`` function::

        myadc.kernel(nroots = 3)
        myadc.analyze()


Algorithms and job control
===========================

The capabilities of the `adc` module at present are summarized in in the
following table:

========== ========== ==================== ===============================
 Method     Reference  Spin-adaptation        Properties
---------- ---------- -------------------- -------------------------------
 ADC(2)     RHF, UHF    Yes                IP, EA, spectroscopic factors, Dyson orb
 ADC(2)-X   RHF, UHF    Yes                IP, EA, spectroscopic factors, Dyson orb
 ADC(3)     RHF, UHF    Yes                IP, EA, spectroscopic factors, Dyson orb
========== ========== ==================== ===============================

The ADC(n) calculations can be performed using different algorithms, depending on
the available memory controlled by the ``max_memory`` keyword:

* In-core

  All tensors such as two-electron integrals and
  amplitudes are stored in memory. This is the default algorithm used when
  sufficient memory is available.


* Out-of-core

  Use of disk to store the expensive tensors.
  This algorithm is invoked by setting ``max_memory`` to a small value.
  See :source:`examples/adc/05-outcore.py`


* Density-fitted (DF) algorithm

 The memory and disk usage can be greatly reduced by approximating the
 two-electron integrals with density-fitting. A simple example of a
 DF-ADC(2) calculation is::

    from pyscf import gto, scf, adc, df
    mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')

    mf = scf.RHF(mol).density_fit('ccpvdz-jkfit').run()
    myadc = adc.ADC(mf).density_fit('ccpvdz-ri')
    eip,vip,pip,xip = myadc.kernel()

More examples can be found in:
:source:`examples/adc/06-dfadc.py`.


References
==========
.. bibliography:: ref_adc.bib
   :style: unsrt
