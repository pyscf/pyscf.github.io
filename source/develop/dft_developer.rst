.. _developer_scf:

*******************************
Density functional theory (DFT)
*******************************

*Modules*: :mod:`dft`, :mod:`pbc.dft`

Overview
========
The DFT methods in modules :mod:`pyscf.dft` and :mod:`pyscf.pbc.dft` 
are implemented as derived classes of the base SCF class :class:`pyscf.scf.hf.SCF`. 
The major modifications of DFT against HF include re-implementation of the 
:meth:`get_veff` and :meth:`energy_elec` methods.
In addition to :meth:`get_jk`, :meth:`get_veff` also calls 
:func:`numint.nr_rks` (for spin-restricted cases) or 
:func:`numint.nr_uks` (for spin-unrestricted cases) 
to compute the exchange-correlation (XC) energy and potential.
The XC energy is added to the total electronic energy in :meth:`energy_elec`.

Numerical integration
=====================
The XC functionals are evaluated on numerical quadrature grids.
These grids are defined in modules :mod:`pyscf.dft.gen_grid` and
:mod:`pyscf.pbc.dft.gen_grid`  for molecules and solids, respectively.
The actual methods to evaluate those XC functionals and their related integrals
are implemented in modules :mod:`pyscf.dft.numint` and :mod:`pyscf.pbc.dft.numint`.
For example, the XC energy and potential matrix for a given density matrix are computed by 
:meth:`nr_rks` (or :meth:`nr_uks`), which internally calls

- :func:`eval_ao` -- to compute the atomic orbitals (AOs) and their derivatives on the grid

- :func:`eval_rho` -- to compute the electron density and density derivatives on the grid

- :func:`eval_xc` -- to compute the XC energy and potential through an interface to the external library Libxc (:func:`pyscf.dft.libxc.eval_xc`)
  or XCFun (:func:`pyscf.dft.xcfun.eval_xc`). 

Other convenient functions implemented in :mod:`numint` include

- :func:`eval_mat` -- evaluating the XC potential matrix with AO, electron density, and XC potential values on the grid as the input 

- :func:`nr_vxc` -- evaluating the XC energy and potential matrix with the density matrix as the input

- :func:`nr_sap_vxc` -- evaluating the superposition of atomic potentials matrix, which is used as the initial guess for :math:`v_{\rm eff}`
  when setting :attr:`mf.init_guess` to ``'vsap'``.

- :func:`nr_rks_fxc`, :func:`nr_uks_fxc` -- evaluating the XC kernel matrix

