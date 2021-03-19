.. _user_agf2:

***********************************************************
Auxiliary second-order Green's function perturbation theory
***********************************************************

*Modules*: :mod:`agf2`

Introduction
============

Auxiliary second-order Green's function perturbation theory (AGF2) 
:cite:`Backhouse2020a,Backhouse2020b`, is a post-Hartree--Fock method
primarily intended for the calculation of ionization potentials (IPs) 
and electron affinities (EAs). AGF2 claculations can be performed with
or without density fitting, and for restricted or unrestricted 
Hartree--Fock references.

Basic usage of AGF2 as given in :source:`examples/agf2/00-simple_agf2.py`
is

.. literalinclude::
    ../../examples/agf2/00-simple_agf2.py

Theory
======

In AGF2, one iteratively solves the Dyson equation

.. math::
    G(\omega) = G_0(\omega) + G_0(\omega) \Sigma(\omega) G(\omega),

using the MP2 self-energy, which possess the frequency-dependent 
spectral function

.. math::
    \Sigma_{xy}(\omega)
    &= \sum_{ija} \frac{ (xi|ja) [ 2 (yi|ja) - (yj|ia) ] }
                       { E_i + E_j - E_a } \\
    &+ \sum_{abi} \frac{ (xa|bi) [ 2 (ya|bi) - (yb|ai) ] }
                       { E_a + E_b - E_i }.

In AGF2, the first two moments of the MP2 self-energy 

.. math::
    T_{xy}^{(0)} &= \sum_{ija} (xi|ja) [ 2 (yi|ja) - (yj|ia) ] \\
                 &+ \sum_{abi} (xa|bi) [ 2 (ya|bi) - (yb|ai) ] \\
    T_{xy}^{(1)} &= \sum_{ija} (xi|ja) [ 2 (yi|ja) - (yj|ia) ]
                    (E_i + E_j - E_a) \\
                 &+ \sum_{abi} (xa|bi) [ 2 (ya|bi) - (yb|ai) ]
                    (E_a + E_b - E_i),

are used in order to express this frequency dependence as a set of 
compressed :math:`2p1h` and :math:`1p2h` states, forming an effective 
single-particle Hamiltonian. This allows Dyson equation to be solved
via the eigenvalue problem

.. math::
    \begin{pmatrix} F & v \\ v^\dagger & \epsilon \end{pmatrix}
    \phi_{n} = \lambda_{n} \phi_{n},

where :math:`F` is the Fock matrix, :math:`\epsilon` the compressed 
states coupling to :math:`F` with coupling strength :math:`v`, and
:math:`\phi_n` are termed the *quasi-molecular orbitals* (QMOs). The
QMOs can be reinserted into the self-energy in place of the 
:math:`i,j,a` and :math:`a,b,i` indices, permitting self-consistency.

Additonally, once projected into the physical space, the QMOs define a
correlated (non-idempotent) density matrix, permiting a further
self-consistency in the Fock matrix in a similar fashion to a
Hartree--Fock calculation. One must ensure that the correct number of
electrons remain in the physical space throughout this process, which
is handled by optimising a chemical potential :math:`\mu`.

Photemission spectra
====================

The compression performed in AGF2 permits the calculation of full
photoemission spectra, as now additional computational effort is
required to calculation large numbers of excitations (i.e. as it would
be in a configuration interaction calculation which employs an
iterative eigensolver). An example of a calculation which provides the
spectral function can be found in 
:source:`examples/agf2/03-photoemission_spectra.py`:

.. literalinclude::
    ../../examples/agf2/03-photoemission_spectra.py

Additionally, the Dyson orbitals corresponding to excitations are
accessibly directly from the method, which can be seen in
:source:`examples/agf2/10-dyson_orbitals.py`.

Restart a calculation
=====================

The contents of an AGF2 calculation can be dumped to the disk via the
familiar PySCF ``chkfile`` utilities. By default, an AGF2 method will
inherit the ``chkfile`` attribute of the underlying Hartree--Fock
reference object. An example of restoring an AGF2 calculation from a
checkpoint file can be found in 
:source:`examples/agf2/04-restart.py`:

.. literalinclude::
    ../../examples/agf2/04-restart.py

Parallelisation
===============

The AGF2 module supports both MPI an OpenMP parallelisation in an aim to
provide a scalable method applicable to interesting problems in quantum
chemistry. Distribution of computational load is handled by the optional
dependency ``mpi4py``, and will run without MPI if an installation cannot
be found.

References
==========

.. bibliography:: ref_agf2.bib
    :style: unsrt
