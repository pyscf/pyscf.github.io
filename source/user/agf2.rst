.. _user_agf2:

******************************************************************
Auxiliary second-order Green's function perturbation theory (AGF2)
******************************************************************

*Modules*: :mod:`agf2`

Introduction
============

Auxiliary second-order Green's function perturbation theory (AGF2)
:cite:`Backhouse2020a,Backhouse2020b` is an iterative, :math:`\mathcal{O}[N^5]`
scaling post-Hartree--Fock method primarily intended for the calculation of
charged excitation spectra, ionisation potentials (IPs) and electron affinities
(EAs), with energetics and single-particle static properties also available.
It can also be considered loosely as an iterative approximate ADC(2) method.
One advantage is that the entire spectrum of charged excitations is found
simultaneously, as the full Green's function is self-consistently optimised.
AGF2 calculations can be performed with or without density fitting, and for
restricted or unrestricted Hartree--Fock references, and is available with 
hybrid parallelism capabilities.

Basic usage of AGF2 as given in :source:`examples/agf2/00-simple_agf2.py`, as

.. literalinclude::
    ../../examples/agf2/00-simple_agf2.py

Theory
======

In AGF2, one iteratively solves the Dyson equation

.. math::
    G(\omega) = G_0(\omega) + G_0(\omega) \Sigma(\omega) G(\omega),

using the self-energy correct through second-order perturbation theory, which
can be written as 

.. math::
    \Sigma_{xy}(\omega)
    &= \sum_{ija} \frac{ (xi|ja) [ 2 (yi|ja) - (yj|ia) ] }
                       { E_i + E_j - E_a } \\
    &+ \sum_{abi} \frac{ (xa|bi) [ 2 (ya|bi) - (yb|ai) ] }
                       { E_a + E_b - E_i }.

In AGF2, the first two spectral moments of the MP2 self-energy,

.. math::
    T_{xy}^{(0)} &= \sum_{ija} (xi|ja) [ 2 (yi|ja) - (yj|ia) ] \\
                 &+ \sum_{abi} (xa|bi) [ 2 (ya|bi) - (yb|ai) ] \\
    T_{xy}^{(1)} &= \sum_{ija} (xi|ja) [ 2 (yi|ja) - (yj|ia) ]
                    (E_i + E_j - E_a) \\
                 &+ \sum_{abi} (xa|bi) [ 2 (ya|bi) - (yb|ai) ]
                    (E_a + E_b - E_i),

are used as a conserved quantity in order to coarse-grain this frequency
dependence as a set of compressed and renormalized :math:`2p1h` and :math:`1p2h`
states. 
These are used to form an effective single-particle Hamiltonian whose size now
only scales with :math:`\mathcal{O}[N]`, rather than :math:`\mathcal{O}[N^3]` if
the full frequency-dependence was maintained. 
This allows Dyson equation to be solved via the hermitian eigenvalue problem

.. math::
    \begin{pmatrix} F & v \\ v^\dagger & \epsilon \end{pmatrix}
    \phi_{n} = \lambda_{n} \phi_{n},

where :math:`F` is the Fock matrix, :math:`\epsilon` the compressed states
coupling to :math:`F` with coupling strength :math:`v`, and :math:`\phi_n` are
termed the *quasi-molecular orbitals* (QMOs), which represent the Dyson orbitals
at each iteration, and span both the original MOs and these auxiliary functions.
The QMOs can be reinserted into the self-energy in place of the :math:`i,j,a`
and :math:`a,b,i` indices, permitting self-consistency to convergence in the 
Green's function with this coarse-grained description of the 
frequency-dependence of the self-energy.

Additionally, once projected into the physical space, the QMOs define a
correlated (non-idempotent) density matrix, permitting a further
self-consistency in the Fock matrix in a similar fashion to a Hartree--Fock 
calculation. 
The correct number of electrons remain in the physical space throughout this
process via optimisation of a chemical potential :math:`\mu`. 
Finally, we note that self-consistency based on higher-order spectral moments of
either the self-energy or resulting Green's function at each iteration is
possible, with a limiting behaviour to the traditional second-order Green's
function method (GF2). 
However, numerical investigations have shown that this higher-order 
self-consistency leads to a deterioration of results, and therefore the 
consistency based on the  first two self-energy spectral moments is the default
behaviour.

Photoemission spectra
====================

The compression of the effective dynamics performed in AGF2 permits the 
calculation of the full spectrum of charged excitations, as no additional 
computational effort is required for this (in contrast to iterative eigensolvers
of effective hamiltonians), and the computational effort is therefore 
independent of the number of excitations requested. 
An example of a calculation which provides the full spectral function can be
found in :source:`examples/agf2/03-photoemission_spectra.py`:

.. literalinclude::
    ../../examples/agf2/03-photoemission_spectra.py

Additionally, the Dyson orbitals corresponding to individual excitations are 
accessible directly from the method, which can be seen in
:source:`examples/agf2/10-dyson_orbitals.py`.

Restart a calculation
=====================

The contents of an AGF2 calculation can be dumped to the disk via the familiar
PySCF ``chkfile`` utilities. 
By default, an AGF2 method will inherit the :attr:`chkfile` attribute of the
underlying Hartree--Fock reference object.
An example of restoring an AGF2 calculation from a checkpoint file can be found
in :source:`examples/agf2/04-restart.py`:

.. literalinclude::
    ../../examples/agf2/04-restart.py

Parallelisation
===============

The AGF2 module supports both MPI an OpenMP parallelisation in an aim to provide
a scalable method applicable to interesting problems in quantum chemistry.
Furthermore, the dominant scaling step is embarrassingly parallel. 
Distribution of computational load is handled by the optional dependency 
``mpi4py``, and will run without MPI using OpenMP threads if an installation
cannot be found.

References
==========

.. bibliography:: ref_agf2.bib
    :style: unsrt
