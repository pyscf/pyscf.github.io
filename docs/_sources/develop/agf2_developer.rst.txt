.. _developer_agf2:

******************************************************************
Auxiliary second-order Green's function perturbation theory (AGF2)
******************************************************************

*Modules*: :mod:`agf2`

Overview
========

The function :func:`pyscf.agf2.AGF2` will direct the method to one of the 
following classes, depending on the nature of the Hartree--Fock object passed as
an argument

===================================== ==================================
:class:`pyscf.agf2.ragf2.RAGF2`       RHF
:class:`pyscf.agf2.uagf2.UAGF2`       UHF
:class:`pyscf.agf2.dfragf2.DFRAGF2`   RHF with density fitting
:class:`pyscf.agf2.dfuagf2.DFUAGF2`   UHF with density fitting
:class:`pyscf.agf2.ragf2_slow.RAGF2`  RHF with arbitrary moment orders
:class:`pyscf.agf2.uagf2_slow.UAGF2`  UHF with arbitrary moment orders
===================================== ==================================

The key attributes of the AGF2 classes include

================== =====================
:attr:`mo_coeff`   saved MO coefficients
:attr:`mo_energy`  saved MO energies
:attr:`mo_occ`     saved MO occupations
:attr:`diis`       DIIS method
:attr:`se`         SelfEnergy object
:attr:`gf`         GreensFunction object
================== =====================

The AGF2 calculation is dispatched via the :func:`kernel` function, which
carries out the following steps:

 - generate the MO basis integrals -- :func:`ao2mo`

 - build the zeroth-order Green's function via the Hartree--Fock states -- :func:`init_gf`

 - build the MP2 self-energy -- :func:`build_se`

 - run the Fock loop (updates the self-energy and Green's function) -- :func:`fock_loop`

 - update the 1-body energy -- :func:`energy_1body`

 - update the self-energy -- :func:`build_se`

 - update the 2-body energy -- :func:`energy_2body`

Auxiliary space objects
=======================

The :attr:`se` and :attr:`gf` attributes are both derived from the
:class:`pyscf.agf2.aux.AuxiliarySpace` object, and are used as
containers for the poles and residues of the self-energy and Green's function. 
The classes contained in :source:`pyscf/agf2/aux.py` are

======================================= ================================
:class:`pyscf.agf2.aux.AuxiliarySpace`  base class for auxiliary spaces
:class:`pyscf.agf2.aux.SelfEnergy`      self-energy
:class:`pyscf.agf2.aux.GreensFunction`  Green's function
======================================= ================================

:class:`AuxiliarySpace` contains most of the attributes and functions
common to both :class:`SelfEnergy` and :class:`GreensFunction`. The key
attributes are

================= ======================================================
:attr:`energy`    array of energies (positions) of the poles
:attr:`coupling`  array of couplings of the poles to a physical space
:attr:`chempot`   position of the Fermi energy
================= ======================================================

Key methods are

===================== ==================================================
:attr:`get_occupied`  get a copy of the object with only occupied poles
:attr:`get_virtual`   get a copy of the object with only virtual poles
:attr:`get_array`     get a dense representation of the Hamiltonian coupling to a physical space
:attr:`dot`           perform a dot-product of the result of :attr:`get_array` with a vector via the sparse representation
:attr:`moment`        return the :math:`n`-th spectral moment of the auxiliaries
===================== ==================================================

Methods specific to the :class:`SelfEnergy` are

============================ ===========================================
:attr:`compress`             compress the self-energy via the AGF2 compression algorithms
:attr:`get_greens_function`  diagonalise the self-energy along with a physical space matrix to get a :attr:`GreensFunction`
============================ ===========================================

Methods specific to the :class:`GreensFunction` are

=========================== ============================================
:attr:`make_rdm1`           get the associated one-particle reduced density matrix
:attr:`real_freq_spectrum`  express the Green's function as a spectrum on a real-valued frequency grid, with a broadening factor
=========================== ============================================

Fock loop
=========

The :func:`fock_loop` function is used to solve the self-consistent
Hartree--Fock-like iterations on the correlated density matrix. 
This step also simultaneously ensures that there is a correct number of electrons in the
physical space.
The chemical potential is first optimised by minimising the metric (in the 
restricted case)

.. math::
    x(\mu) = \Big( N_\mathrm{elec} - 
                   2 \sum_{p}^{n_\mathrm{MO}} 
                     \sum_{i}^{n_\mathrm{QMO}^\mathrm{occ}}
                     \phi_{pi} \phi_{pi}^*
             \Big)^{2},

where the chemical potential (:math:`\mu`) adjusts the positions of the poles of
the self-energy before diagonalisation, thereby adjusting the weight of the
occupied QMOs on the physical space and in turn the trace of the one-particle
density matrix.
The gradient of :math:`x` with respect to the chemical potential
:math:`\frac{\partial x(\mu)}{\partial \mu}` is also implemented such that a
truncated Newton algorithm can be used via :func:`scipy.optimize.minimize`.
The functions required to handle this step are contained in
:source:`pyscf/agf2/chempot.py`.

Memory requirements
===================

Due to the renormalisation of the second-order diagrams, the MO basis integrals
must be transformed into the QMO basis at subsequent iterations. 
This requires the combinations :math:`(xo|ov)` and :math:`(xv|vo)`, where 
:math:`x` are general MOs, :math:`o` are occupied QMOs and :math:`v` are virtual
QMOs.
The occupied QMOs number :math:`n_\mathrm{MO} + n_\mathrm{MO}^\mathrm{occ}`, and
virtual :math:`n_\mathrm{MO} + n_\mathrm{MO}^\mathrm{vir}`, and so these arrays
become prohibitively expensive with increasing system size.
This necessitates the use of density fitting for large system sizes, which
provides lower scaling memory requirements.
The density fitting implementation still requires the four-centre integrals to
be built in order to compute the moments of the self-energy, but they may be
computed block-wise in order that the memory overhead remains low.

Arbitrary order moment calculations
===================================

The :class:`pyscf.agf2.ragf2_slow.RAGF2` and 
:class:`pyscf.agf2.uagf2_slow.UAGF2` classes support AGF2 calculations beyond
the standard efficient AGF2(1,0) implementation.
These classes are largely unoptimised and do not support density fitting, and
are intended for developmental use only.
For a discussion of the order of the moments required as parameters to these
methods, please refer to our paper :cite:`Backhouse2020a`.
An example of higher-order AGF2 calculations can be found in
:source:`pyscf/examples/agf2/05-agf2_moments.py`, with a connection to the 
algebraic diagrammatic construction method.

References
==========

.. bibliography:: ../user/ref_agf2.bib
    :style: unsrt
