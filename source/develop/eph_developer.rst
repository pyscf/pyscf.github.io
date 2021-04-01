.. _developer_eph:

************************
Electron Phonon Coupling
************************

*Modules*: :mod:`eph`, :mod:`pbc.eph`

Overview
========

All EPH classes in :mod:`pyscf.eph` are implemented as derived class of the base Hessian class, eg, :class:`pyscf.hessian.rhf.Hessian`. These include

=================================== =====
:class:`pyscf.eph.rhf.EPH`           RHF
:class:`pyscf.eph.uhf.EPH`           UHF
:class:`pyscf.eph.rks.EPH`           RKS
:class:`pyscf.eph.uks.EPH`           UKS
=================================== ======

The key attributes of the EPH class include

===========================     ====================================
:attr:`cutoff_frequency`        cutoff frequency in wavenumber unit
:attr:`keep_imag_frequency`     whether to store the imaginary modes
===========================     ====================================

The main driver is the :func:`pyscf.eph.rhf.kernel` function
which takes a EPH object and carries out the following steps:

- building the nuclear hessian, which consists of

  - building the e-e part -- :func:`hess_elec`

  - building the e-nuc part -- :func:`hess_nuc`

- computing the vibrational frequencies and polarization vectors -- :func:`get_mode`, which calls :func:`solve_hmat` for diagonalization

- building the electron phonon coupling matrix -- :func:`get_eph`, which consists of

  - computing the e-nuc part -- :func:`vnuc_generator`

  - computing the e-e part from MO derivative -- :func:`rhf_deriv_generator`

  - computing the e-e part from orbital derivative -- :func:`_get_jk` and additionally :func:`_get_vxc_deriv1` for dft

  - assembling the response matrix with eigenmodes/eigenvectors to the final matrix form

Internally, different methods reuse this kernel by overwriting the
Hessian methods.

- DFT SCF is implemented by
  specializing :meth:`hess_elec` and :meth:`hess_nuc` for the various KS objects
  (see e.g., :func:`pyscf.hessian.rks.hess_elec` and :func:`pyscf.hessian.rks.hess_nuc`).
