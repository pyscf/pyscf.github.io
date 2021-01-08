.. _adc:

:mod:`adc` --- Algebraic diagrammatic construction (ADC) scheme
***************************************************************

.. module:: adc
   :synopsis: Algebraic diagrammatic construction energies and propties

The :mod:`adc` module is an implementation of algebraic diagrammatic construction theory for computation of spectroscopic properties. At the present moment the :mod:`adc` module supports computation of ionization potentials (IP) and electron affinties (EA) using three ADC approximations: ADC(2), ADC(2)-X, and ADC(3).

Computations using the :mod:`adc` module should be performed in three steps:

1. Carry out a mean-field computation using the :mod:`scf` module to obtain ground-state (reference) molecular orbitals.

2. Compute the ADC(n) ground-state energy and amplitudes. For the ADC(n) method, the ground-state energy and amplitudes
   correspond to those of the :math:`n`-th order Møller--Plesset perturbation theory (MPn).

3. Compute the ADC(n) spectroscopic properties (e.g., electron attachment energies
   and the corrresponding transition intensities)

An example of IP and EA computation using the default ADC method (ADC(2)) is shown below::

        from pyscf import gto, scf, adc
        mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')
        mf = scf.RHF(mol).run()
        myadc = adc.ADC(mf)
        myadc.kernel()
        e_ip,v_ip,p_ip = myadc.ip_adc(nroots=3)
        e_ea,v_ea,p_ea = myadc.ea_adc(nroots=3)

Here, the kernel() function executes computation of the ADC(2) ground-state energy and amplitudes, 
which are used to perform calculation of the three lowest-energy IP and EA. 
Other methods of the ADC family can be requested by setting the method variable to ``ADC(2)-X`` or ``ADC(3)``, e.g.::

        from pyscf import gto, scf, adc
        mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')
        mf = scf.RHF(mol).run()
        myadc = adc.ADC(mf)
        myadc.method = "adc(3)"
        myadc.kernel()
        e_ip,v_ip,p_ip = myadc.ip_adc(nroots=1)
        e_ea,v_ea,p_ea = myadc.ea_adc(nroots=1)

An alternative way to run ADC computations is to directly call the kernel function of the corresponding ADC classes::

        myadc = adc.ADC(mf)
        myadc.kernel()
        myadcip = adc.uadc.UADCIP(myadc)
        myadcea = adc.uadc.UADCEA(myadc)
        e_ip,v_ip,p_ip = myadcip.kernel(nroots=3)
        e_ea,v_ea,p_ea = myadcea.kernel(nroots=3)

The ADC implementation in Pyscf can be used to perform calculations of electron affinities (EA) and ionization potentials (IP) for closed- and open-shell molecules starting with the RHF, ROHF, or UHF reference SCF wavefunctions, e.g.::

        from pyscf import gto, scf, adc
        mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz', charge = 1, spin = 1)
        mf = scf.UHF(mol).run()
        myadc = adc.ADC(mf)
        myadc.method = "adc(3)"
        myadc.kernel()
        e_ip,v_ip,p_ip = myadc.ip_adc(nroots=1)
        e_ea,v_ea,p_ea = myadc.ea_adc(nroots=1)

When using results of this code for publications, please cite the following paper: "Third-order algebraic diagrammatic construction theory for electron attachment and ionization energies: Conventional and Green’s function implementation", S. Banerjee and A.Y. Sokolov, J. Chem. Phys. 151, 224112 (2019).

Examples
========

* :source:`examples/adc/01-closed_shell.py`
* :source:`examples/adc/02-open_shell.py`
* :source:`examples/adc/03-closed_shell_different_setup.py`

Program reference
=================

.. automodule:: pyscf.adc
   :members:

adc.uadc module and UADC class
--------------------------------

The :class:`pyscf.adc.uadc.UADC` class is the object to hold the unrestricted ADC environment
attributes and results. The environment attributes are the parameters for the ground state as well as the excited state calculations. 

.. automodule:: pyscf.adc.uadc
   :members:
