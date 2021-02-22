.. _developer_scf:

*******************************
Density functional theory (DFT)
*******************************

*Modules*: :mod:`dft`, :mod:`pbc.dft`

Overview
========
The DFT methods in modules :mod:`pyscf.dft` and :mod:`pyscf.pbc.dft` 
are implemented as derived classes of the base SCF class :class:`pyscf.scf.hf.SCF`. 
The major modifications of DFT with respect to HF include the re-implementation of the 
:meth:`get_veff` and :meth:`energy_elec` methods.
In addition to :meth:`get_jk`, :meth:`get_veff` also calls 
:func:`numint.nr_rks` (for spin-restricted cases) or 
:func:`numint.nr_uks` (for spin-unrestricted cases) 
to compute the exchange-correlation (XC) energy and potential.
The XC energy is added to the total electronic energy in :meth:`energy_elec`.

The key attributes of the KS classes include

==========================  ===== 
:attr:`xc`                  names of the XC functionals
:attr:`nlc`                 names of the non-local correlation functionals
:attr:`omega`               :math:`\omega` of the range-separated Coulomb operator :math:`e^{-\omega r_{12}^2} / r_{12}`
:attr:`grids`               :class:`dft.gen_grid.Grids` (or :class:`pbc.dft.gen_grid.Grids`) objects which defines the DFT grids
==========================  =====

Numerical grids
===============
The XC functionals are evaluated on numerical quadrature grids.
These grids are defined in modules :mod:`dft.gen_grid` and
:mod:`pbc.dft.gen_grid` for molecules and solids, respectively.
Both uniform grids and Becke (atomic) grids are supported.

The base class for atomic grids is :class:`pyscf.dft.gen_grid.Grids`, which include
the following key attributes

======================  ====
:attr:`level`           grid level to control the number of radial and angular grids
:attr:`atomic_radii`    atomic radii, can be :attr:`dft.radi.BRAGG_RADII` (Bragg-Slater radii), 
                        :attr:`dft.radi.COVALENT_RADII` (covalent radii) or ``None`` (no atomic size adjustment in the atomic partition scheme of Becke)
:attr:`radii_adjust`    functions for atomic size adjustment, can be :func:`dft.radi.becke_atomic_radii_adjust`,
                        :func:`dft.radi.treutler_atomic_radii_adjust` or ``None``
:attr:`radi_method`     functions for radial grid schemes, can be :func:`dft.radi.treutler`, :func:`dft.radi.delley`, 
                        :func:`dft.radi.mura_knowles` or :func:`dft.radi.gauss_chebyshev`
:attr:`becke_scheme`    functions for atomic partition schemes, can be :func:`dft.gen_grid.original_becke` or :func:`dft.gen_grid.stratmann` 
:attr:`prune`           functions for grid pruning, can be :func:`dft.gen_grid.nwchem_prune`, :func:`dft.gen_grid.sg1_prune`,
                        :func:`dft.gen_grid.treutler_prune` or ``None`` (no pruning)
:attr:`atom_grid`       user defined number of radial and angular grids for specific atom types
:attr:`coords`          saved Cartesian coordinates of grid points
:attr:`weights`         saved weights of grid points
======================  ====

To generate the Becke grids, one can simply initialize the :class:`Grids` object and then call the :meth:`build` method::

    g = Grids()
    g.build()

Internally, :meth:`build` calls the following functions:

- :func:`pyscf.dft.gen_grid.gen_atomic_grids` -- which generates the coordinates and weights of the grid points 
  with respect to the atom center for each atom type

- :func:`pyscf.dft.gen_grid.get_partition` -- which generates the coordinates and weights of the molecular grid points

See :source:`examples/dft/11-grid_scheme.py` for examples of specifying various grid schemes.

The PBC implementation of the Becke grids is defined in class :class:`pyscf.pbc.dft.gen_grids.BeckeGrids`,
which is derived from the base class :class:`pyscf.dft.gen_grid.Grids`.
The main modification of the PBC implementation is that the grid points
contained in the reference unit cell (including those belonging to the
periodic images) are all included. (see :func:`pyscf.pbc.dft.gen_grids.get_becke_grids`)

The uniform grid is implemented in :class:`pyscf.pbc.dft.gen_grids.UniformGrids`, whose
:meth:`build` method internally calls :func:`pyscf.pbc.gto.cell.get_uniform_grids`.

Numerical integration
=====================
The actual methods to evaluate XC functionals and their related integrals
are implemented in modules :mod:`pyscf.dft.numint` and :mod:`pyscf.pbc.dft.numint`.
For example, the XC energy and potential matrix for a given density matrix are computed by 
:meth:`nr_rks` (or :meth:`nr_uks`), which internally calls

- :func:`eval_ao` -- to compute the atomic orbitals (AOs) and their derivatives on the grid

- :func:`eval_rho` -- to compute the electron density and density derivatives on the grid

- :func:`eval_xc` -- to compute the XC energy and potential through an interface to the external library Libxc (:func:`pyscf.dft.libxc.eval_xc`)
  or XCFun (:func:`pyscf.dft.xcfun.eval_xc`). 

Other useful functions implemented in :mod:`numint` include

- :func:`eval_mat` -- evaluating the XC potential matrix with AO, electron density, and XC potential values on the grid as the input 

- :func:`nr_vxc` -- evaluating the XC energy and potential matrix with the density matrix as the input

- :func:`nr_sap_vxc` -- evaluating the superposition of atomic potentials matrix, which is used as the initial guess for :math:`v_{\rm eff}`
  when setting :attr:`mf.init_guess` to ``'vsap'``.

- :func:`nr_rks_fxc`, :func:`nr_uks_fxc` -- evaluating the XC kernel matrix

Using values on the grid
========================
A few examples of evaluating various quantities on numerical grids are as follows:

- :source:`examples/dft/30-ao_value_on_grid.py`

- :source:`examples/dft/31-xc_value_on_grid.py`

- :source:`examples/pbc/30-ao_value_on_grid.py`

- :source:`examples/pbc/30-overlap_periodic_cell.py`

- :source:`examples/gto/24-ao_value_on_grid.py`
