.. _theory_gw:

GW approximation
****************

*Modules*: :mod:`gw`, :mod:`pbc.gw`


Introduction
============

The GW approximation is a Green's function-based method that calculates charged
excitation energies, i.e. ionizations potentials (IPs) and electron affinities
(EAs).  PySCF implements the G\ :sub:`0`\ W\ :sub:`0` approximation, in which the
self-energy is built with mean-field orbitals and orbital energies.  Therefore,
the results depend on the mean-field starting point, which can be Hartree-Fock
or density functional theory.  As described below, PySCF has three
implementations of the GW approximation, all of which are "full-frequency".

An example GW calculation is shown below:

.. literalinclude:: ../../examples/gw/00-simple_gw.py

In this example, the ``orbs`` keyword argument is
provided to select which GW orbital energies are requested.  By default,
all orbital energies (occupied and unoccupied) will be calculated, which
will increase the cost and may have errors, depending on the method
of frequency integration.

Frequency integration
=====================

The frequency integration needed for the GW approximation can be done in three
ways, controlled by the ``freq_int`` keyword argument: by analytic continuation 
(AC, ``freq_int='ac'``), contour deformation (CD, ``freq_int='cd'``), or
exactly (Exact, ``freq_int='exact'``).  The first two are much more affordable 
and typically provide sufficient accuracy.  GW-AC supports spin-restricted and
spin-unrestricted calculations; GW-CD and GW-Exact only support spin-restricted
calculations.

Analytic continuation
---------------------

Integration via analytic continuation is implemented in the ``GWAC`` module
that is accessed with ``freq_int='ac'``, which is also
the default ``GW`` module.  GW-AC has :math:`N^4` scaling and is recommended for
valence states only.
The analytic continuation can be done using a Pade
approximation (default, more reliable) or a two-pole model, controlled by the ``ac``
attribute.
GW-AC supports frozen core orbitals for reducing computational cost, 
controlled by the ``frozen`` attribute (number of frozen core MOs neglected in GW-AC
calculation). Frozen core orbitals are not supported by other GW methods currently.
There are two ways to compute GW orbital energies, controlled by the ``linearized`` attribute: 
``linearized=False`` (default) solves the quasiparticle equation through a Newton solver 
self-consistently, while ``linearized=True`` employs a linearization approximation::

    mygw = gw.GW(mf) # same as freq_int='ac' or GWAC module
    # mygw.ac = 'pade' # default
    mygw.ac = 'twopole'
    mygw.frozen = 1 # default is None
    mygw.linearized = False # default
    mygw.kernel(orbs=range(nocc-3,nocc+3))

Contour deformation
-------------------

Integration via contour deformation is implemented in the ``GWCD`` module
that is accessed with ``freq_int='cd'``.
GW-CD has :math:`N^4` scaling and is slower, but more robust, than GW-AC.
GW-CD is particularly recommended for accurate core and high-energy states::

    mygw = gw.GW(mf, freq_int='cd').run(orbs=[0,1])

Exact
-----

Exact frequency integration can be carried out analytically and is implemented
in the ``GWExact`` module that is accessed with ``freq_int='exact'``.  Exact integration 
requires complete diagonalization of the RPA matrix, which has :math:`N^6`
scaling.  However, all orbital energies can be readily obtained without error::

    mygw = gw.GW(mf, freq_int='exact').run()

By default, GW-Exact, like the other GW implementations, use the direct random-phase 
approximation (dRPA) to screen the Coulomb interaction.
Within GW-Exact, any alternative time-dependent mean-field theory (TDHF, TDDFT, etc.) can
be also used.  The instance of an executed ``tdscf`` method can be provided as a
keyword argument:

.. literalinclude:: ../../examples/gw/10-custom_screening.py
