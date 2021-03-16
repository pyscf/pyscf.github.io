.. _theory_adc:

************************************************
Algebraic diagrammatic construction (ADC) scheme
************************************************

*Modules*: :mod:`adc`

Introduction
=============================
The algebraic diagrammatic construction theory (ADC(n)) is a post-Hartree-Fock
method used for computing correlated excited states of molecules.
The ADC methods involve perturbative expansion of a propagator followed
by truncation at a certain order :math:`n`. The ADC(n) methods provide access to
excitation energies and their corresponding transition intensities via a
'one-shot' calculation. At present, the `adc` module in PySCF contains different
variants of the method (ADC(2), ADC(2)-X and ADC(3)) that can be used to calculate
charged excitations such as ionization potentials (IP-ADC) and
electron affinities (EA-ADC).

A simple ADC calculation involves the calculation of ground-state ADC(n) energy
and amplitudes that correspond to those of the :math:`n`-th order
Møller--Plesset perturbation theory (MPn), followed by the calculation of the
ADC(n) spectroscopic properties (e.g., ionization/electron-attachment energies
and their corresponding transition intensities).

An example of IP and EA computation using the default ADC method (ADC(2))
is shown below::

        from pyscf import gto, scf, adc
        mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')
        mf = scf.RHF(mol).run()
        myadc = adc.ADC(mf)
        myadc.kernel_gs()
        e_ip,v_ip,p_ip,x_ip,adc_es = myadc.ip_adc()
        e_ea,v_ea,p_ea,x_ea,adc_ea = myadc.ea_adc()

Alternatively, one can specify the type of charged excitation and order of
the desired ADC computation::

        myadc.method = "adc(3)"
        myadc.method_type = "ea"
        myadc.kernel()

The ADC functions return the ``nroots`` eigenvalues with the lowest energy. The
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
See :ref:`theory_scf` to know more about the different reference wavefunctions.

Shown below is an example of a IP-/EA-UADC(2) calculation::

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


Properties
==========
A number of properties such as spectroscopic amplitudes and spectroscopic factors
are available at the ADC(n) level via the flag ``compute_properties``, which is
is set to True by default::

        myadc.compute_properties = True
        myadc.method = "adc(3)"
        myadc.method_type = "ip"
        myadc.kernel(nroots = 3)

Additionally, the `adc` module can be used to compute Dyson orbitals corresponding
to ionized and electron-attached states::

        dyson_orb = myadc.compute_dyson_mo()


Analysis of spectroscopic properties
=====================================
The `adc` module provides analysis of the ADC(n) eigenvectors by default,
whereas an additional analysis of spectroscopic factors is invoked when the flag
``compute_properties`` is set to True. The analysis of the ADC(n) eigenvectors
and spectroscopic factors can be performed as::

        myadc.kernel(nroots = 3)
        myadc.analyze()


The capabilities of the `adc` module at present are summarized in in the
following table:

========== ========== ==================== ===============================
 Method     Reference  Spin-adaptation        Properties
---------- ---------- -------------------- -------------------------------
 ADC(2)     RHF, UHF    Yes                IP, EA, spec factors, Dyson orb
 ADC(2)-X   RHF, UHF    Yes                IP, EA, spec factors, Dyson orb
 ADC(3)     RHF, UHF    Yes                IP, EA, spec factors, Dyson orb
========== ========== ==================== ===============================


Algorithms and job control
===========================

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

 Memory and disk usage greatly reduced by approximating
 two-electron integrals. A simple example of a
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
