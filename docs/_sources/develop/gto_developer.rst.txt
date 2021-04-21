.. _developer_gto:

Molecules and Crystal Unit Cells
********************************

*Modules*: :mod:`gto`, :mod:`pbc.gto`

Overview
--------

The :class:`pyscf.gto.Mole` class holds the information of 
molecular structures, basis sets, effective core potentials
and other global options. It also implements
interface functions to obtain AO integrals computed by 
`libcint <https://github.com/sunqm/libcint>`_ library.

The :class:`pyscf.pbc.gto.Cell` class, which is derived from the
:class:`~pyscf.gto.Mole` class, holds additional information for calculations of
extended systems with periodic boundary conditions (PBCs), such as
crystal lattice vectors, pseudopotentials, periodic dimensions, etc.

Internal data structure
-----------------------

The :class:`Mole` class has a three-layer data structure: 
user input, internal format, libcint argument.

===============   ================  =====================
user input        internal format   libcint argument
===============   ================  =====================
:attr:`.atom`     :attr:`._atom`    :attr:`._atm`, :attr:`._env`
:attr:`.basis`    :attr:`._basis`   :attr:`._bas`, :attr:`._env`
:attr:`.ecp`      :attr:`._ecp`     :attr:`._ecpbas`, :attr:`._env`
:attr:`.pseudo`   :attr:`._pseudo`  see :func:`pyscf.pbc.gto.pseudo.pp.get_pp`
===============   ================  =====================

The attributes :attr:`~Mole.atom` and :attr:`~Mole.basis` store the user input for
molecular geometries and basis sets, respectively. These are then converted 
to the internal formats through the methods :meth:`~pyscf.gto.Mole.format_atom` and 
:meth:`~pyscf.gto.Mole.format_basis`, respectively::

  >>> from pyscf import gto
  >>> mol = gto.Mole()
  >>> mol.atom = 'H 0 0 0; H 0 0 0.74'
  >>> mol.basis = {"H": "6-31g"}
  >>> print(mol.atom)
  H 0 0 0; H 0 0 0.74
  >>> print(mol.basis)
  {'H': '6-31g'}
  >>> mol._atom = mol.format_atom(mol.atom)
  >>> print(mol._atom)
  [('H', [0.0, 0.0, 0.0]), ('H', [0.0, 0.0, 1.3983973321781458])]
  >>> mol._basis = mol.format_basis(mol.basis)
  >>> print(mol._basis)
  {'H': [[0, [18.731137, 0.0334946], [2.8253937, 0.23472695], [0.6401217, 0.81375733]], [0, [0.1612778, 1.0]]]}

Finally, the data stored in the internal format layer is converted,
via the :meth:`~pyscf.gto.Mole.make_env` method, to the format suitable for libcint input::

  >>> mol._atm, mol._bas, mol._env = mol.make_env(mol._atom, mol._basis, mol._env[:gto.mole.PTR_ENV_START])
  >>> print(mol._atm)
  [[ 1 20  1 23  0  0]
   [ 1 24  1 27  0  0]]
  >>> print(mol._bas)
  [[ 0  0  3  1  0 28 31  0]
   [ 0  0  1  1  0 34 35  0]
   [ 1  0  3  1  0 28 31  0]
   [ 1  0  1  1  0 34 35  0]]
  >>> print(mol._env)
  [ 0.          0.          0.          0.          0.          0.
    0.          0.          0.          0.          0.          0.
    0.          0.          0.          0.          0.          0.
    0.          0.          0.          0.          0.          0.
    0.          0.          1.39839733  0.         18.731137    2.8253937
    0.6401217   0.76192622  1.292371    1.471319    0.1612778   0.64297778]

Each row of :attr:`~Mole._atm` represents an atom, and 
each column has the meaning as follows:

============   =======
Column index   Content
============   =======
0              Nuclear charge
1              Index of :attr:`_env` where Cartesian coordinates of atoms are stored 
               (see :func:`~pyscf.gto.Mole.atom_coords` and :func:`~pyscf.gto.Mole.bas_coord`)
2              Label of nuclear charge models (``1``: point charge; ``2``: Gaussian distribution)
3              Index of :attr:`_env` where the exponent (zeta) of the Gaussian charge model is stored 
               (see :func:`~pyscf.gto.Mole.set_nuc_mod`)
============   =======

Each row of :attr:`~Mole._bas` represents a basis shell, and
each column has the following meaning:

============   =======
Column index   Content
============   =======
0              Index of atoms
1              Angular momentum (:func:`~pyscf.gto.Mole.bas_angular`)
2              Number of primative functions (:func:`~pyscf.gto.Mole.bas_nprim`)
3              Number of contractions (:func:`~pyscf.gto.Mole.bas_nctr`)
4              Kappa (:func:`~pyscf.gto.Mole.bas_kappa`)
5              Index of :attr:`_env` where the exponents of basis functions are stored 
               (the exponents can be obtained by calling :func:`~pyscf.gto.Mole.bas_exp`)
6              Index of :attr:`_env` where the contraction coefficients of basis functions are stored 
               (the contraction coefficients can be obtained by calling :func:`~pyscf.gto.Mole.bas_ctr_coeff`)
============   =======

The :attr:`_env` attribute holds the parameters of basis functions, ECPs and 
other useful information:

======  =======
Offset  Content
======  =======
0       Cutoff of basis funtion exponents (see :func:`~pyscf.gto.Mole.with_integral_screen`)
1       Common origin for integrals of dipole, :math:`r \times p`, etc. (see :func:`~pyscf.gto.Mole.with_common_origin`)
4       Origin of :math:`1 / r` operator (see :func:`~pyscf.gto.Mole.with_rinv_origin`)
7       Exponent of Gaussian charge distribution (see :func:`~pyscf.gto.Mole.with_rinv_zeta`)
8       :math:`\omega` for range-separated Coulomb operator (see :func:`~pyscf.gto.Mole.with_range_coulomb`)
9       see :func:`~pyscf.gto.Mole.set_f12_zeta`
17      see :func:`~pyscf.gto.Mole.with_rinv_at_nucleus`
18      Offset for ECP parameters (see :func:`~pyscf.gto.Mole.intor`)
19      Length of :attr:`~Mole._ecpbas` (see :func:`~pyscf.gto.Mole.intor`)
20      Parameters of basis functions, ECPs, etc.
======  =======

Access AO integrals
-------------------

The :class:`Mole` class provides a few functions to access the integrals computed by 
lincint. These include:

========================================  =======
Function                                  Comment
:func:`~pyscf.gto.Mole.intor`             General integral generator
:func:`~pyscf.gto.Mole.intor_symmetric`   Generator for Hermitian 1-e integrals
:func:`~pyscf.gto.Mole.intor_asymmetric`  Generator for anti-Hermitian 1-e integrals
:func:`~pyscf.gto.Mole.intor_by_shell`    Generator for integrals of specific shells
:func:`pyscf.gto.mole.intor_cross`        Generator for 1-e integrals between two :class:`Mole` objects
========================================  =======

These functions internally call :func:`pyscf.gto.moleintor.getints` and 
:func:`pyscf.gto.moleintor.getints_by_shell` to interface with libcint.

AO evaluation on grids
----------------------

The AO functions can be evaluated on grids by the :func:`~pyscf.gto.Mole.eval_ao` function.

Serialization
-------------

The :class:`Mole` object can be serialized by the :func:`~pyscf.gto.Mole.dumps`
function into a string with JSON format, which is used by the :mod:`pyscf.lib.chkfile`
module. And it can be deserialized by the :func:`~pyscf.gto.Mole.loads` function.

Crystal unit cell
-----------------

The :class:`~pyscf.pbc.gto.Cell` class is defined as an extension of the
:class:`~pyscf.gto.Mole` class. They share most of the data structures and methods.
This gives the freedom to mix finite-size and PBC calculations (see :ref:`mix_mol`).

.. note::

  The serialization methods of the two classes are not completely compatible.  
  If the :func:`~pyscf.gto.Mole.dumps` and :func:`~pyscf.gto.Mole.loads` functions 
  in the :class:`Mole` class are called for a :class:`Cell` object,
  some information of the :class:`Cell` object may be lost.

Besides the methods and parameters provided by the :class:`Mole` class (see :ref:`gto`),
there are some others frequently used in the PBC code to access the
information of the crystal.

* :attr:`Gv` -- (N x 3) array

  - The plane wave bases for 3D-FFT transformation. 
    Given ``cell.mesh = [nx,ny,nz]``, the number of plane waves is ``N=nx*ny*nz``.
    ``Gv`` is obtained by the method :meth:`~pyscf.pbc.gto.Cell.get_Gv` 
    (or simply :meth:`Cell.Gv`).

* :attr:`vol` -- float

  - :attr:`Cell.vol` gives the volume of the unit cell (in atomic unit).

* :func:`~pyscf.pbc.gto.Cell.reciprocal_vectors` -- returns a (3 x 3) array

  - Each row is a reciprocal space primitive vector.

* :func:`~pyscf.pbc.gto.Cell.energy_nuc` (or :func:`~pyscf.pbc.gto.cell.ewald`)

  - Similar to :func:`~pyscf.gto.Mole.energy_nuc`. 
    The nuclear repulsion energy is computed with Ewald summation.
    It depends on three parameters: the truncation radius for real-space lattice summation
    :attr:`rcut`, the Gaussian model charge :attr:`ew_eta`, and the energy cutoff
    :attr:`ew_cut`. And these are determined by function :func:`~pyscf.pbc.gto.cell.get_ewald_params`.

* :func:`~pyscf.pbc.gto.Cell.pbc_intor`

  - PBC analytic integral driver. Note that the :meth:`~pyscf.gto.Mole.intor` method 
    is not overloaded in the :class:`Cell` class. This allows one to compute  
    both the periodic and open-boundary integrals within the :class:`Cell` object.
