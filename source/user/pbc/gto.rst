.. _user_pbc_gto:

Crystal structure
*****************

*Modules*: :mod:`pbc.gto`

Initializing a crystal
======================

Initializing a crystal unit cell is very similar to the initialization
of a molecule. Instead of the :class:`gto.Mole` class, 
one uses the :class:`pbc.gto.Cell` class to define a cell::

  >>> from pyscf.pbc import gto
  >>> cell = gto.Cell()
  >>> cell.atom = '''H  0 0 0; H 1 1 1'''
  >>> cell.basis = 'gth-dzvp'
  >>> cell.pseudo = 'gth-pade'
  >>> cell.a = numpy.eye(3) * 2
  >>> cell.build()

The other two equivalent ways to initialize a molecule 
introduced in :numref:`user_gto` also apply here::

  >>> from pyscf.pbc import gto
  >>> cell = gto.Cell()
  >>> cell.build(
  ...   atom = '''H  0 0 0; H 1 1 1''',
  ...   basis = 'gth-dzvp',
  ...   pseudo = 'gth-pade',
  ...   a = numpy.eye(3) * 2)

  >>> import pyscf
  >>> cell = pyscf.M(
  ...   atom = '''H  0 0 0; H 1 1 1''',
  ...   basis = 'gth-dzvp',
  ...   pseudo = 'gth-pade',
  ...   a = numpy.eye(3) * 2)

  >>> from pyscf.pbc import gto
  >>> cell = gto.M(
  ...   atom = '''H  0 0 0; H 1 1 1''',
  ...   basis = 'gth-dzvp',
  ...   pseudo = 'gth-pade',
  ...   a = numpy.eye(3) * 2)


Lattice vectors
---------------

The crystal initialization requires an additional parameter :attr:`~Cell.a`, which
represents the lattice vectors. The format of :attr:`~Cell.a` is array-like::
    
    cell.a = numpy.eye(3) * 2
    cell.a = [[2,0,0],[0,2,0],[0,0,2]]

Each row of the 3-by-3 matrix of :attr:`~Cell.a` represents a lattice vector
in Cartesian coordinate, with the same unit as the input :attr:`atom` parameter.

.. note:: The input lattice vectors (row by row) should form a right-handed coordinate system, as otherwise some integrals may be computed incorrectly in PySCF.

Basis set and pseudopotential
-----------------------------

PySCF uses the crystalline Gaussian-type orbitals as the basis functions
for solid calculations. The predefined such basis sets include 
the valence basis sets that are optimized for the GTH pseudopotentials 
(a whole list can be found in :source:`pyscf/pbc/gto/basis` and :source:`pyscf/pbc/gto/pseudo`).
The input format of `basis sets`_ for the :class:`Cell` object is the same 
as that for the :class:`Mole` object.
In addition, the predefined basis sets and ECPs for molecular calculations 
can be used in solid calculations as well

.. literalinclude:: /../examples/pbc/05-input_pp.py

Finally, custom basis sets can be defined just like that in molecular calculations

.. literalinclude:: /../examples/pbc/04-input_basis.py

Low-dimentional systems
-----------------------

The PySCF  :mod:`pbc` module also supports low-dimensional periodic systems. You can initialize
the attribute :attr:`Cell.dimension` to specify the dimension of the system::

  >>> cell.dimension = 2
  >>> cell.a = numpy.eye(3) * 2
  >>> cell.build()

When :attr:`~Cell.dimension` is smaller than 3, a vacuum of infinite size will be
applied in certain direction(s).  For example, when :attr:`~Cell.dimension`
is 2, the z-direction will be treated as infinitely large and the xy-plane
constitutes the periodic surface. When :attr:`~Cell.dimension` is 1, y and z axes
are treated as vacuum thus wire is placed in the x direction. 
When :attr:`~Cell.dimension` is 0, all three directions are treated as vacuum, and this is
equivalent to a molecular system.

K-points
--------

The k-points used in solid calculations can be obtained through the 
:meth:`Cell.make_kpts` method. The minial input is to specify the number of k-points
in each lattice vector direction::
    
    >>> kpts = cell.make_kpts([1,2,2])
    >>> print(kpts.shape)
    (4,3)

By default, this will return the shifted Monkhorst-pack mesh which is centered at 
the :math:`Gamma`-point. To get non-shifted Monkhorst-pack mesh, one can call::

    >>> kpts = cell.make_kpts([1,2,2], with_gamma_point=False)

To get shifted Monkhorst-pack mesh centered at a specific point, one can call::

    >>> kpts = cell.make_kpts([1,2,2], scaled_center=[0.,0.25,0.25])

where ``scaled_center`` is defined in the units of lattce vectors.

The obtained k-points are used as input for crystalline electronic structure calculations.
For example, k-point sampled RHF can be invoked as follows::

    >>> from pyscf.pbc import scf
    >>> kpts = cell.make_kpts([2,2,2])
    >>> kmf = scf.KRHF(cell, kpts = kpts)
    >>> kmf.kernel()

More details about k-point sampling for each method can be found in the corresponding chapter.

Other parameters
----------------
The other attributes of the :class:`Mole` class such as :attr:`verbose`,
:attr:`max_memory`, :attr:`spin`, etc., have the same meanings in the :class:`Cell` class.

.. note:: Currently, point group symmetries for crystals are not supported.

Accessing AO integrals
======================

The :func:`Mole.intor` method can only be used to evaluate the integrals with open boundary
conditions. When the periodic boundary conditions of crystalline systems are
studied, one needs to use the :func:`pbc.Cell.pbc_intor` function to evaluate the
integrals for short-range operators, such as the overlap and kinetic matrix::

  overlap = cell.pbc_intor('int1e_ovlp')
  kin = cell.pbc_intor('int1e_kin')

By default, the :func:`pbc.Cell.pbc_intor` function only returns the :math:`\Gamma`-point
integrals. If k-points are specified, it will evaluate the integrals at each k-point::

  kpts = cell.make_kpts([2,2,2])  # 2x2x2 Monkhorst-pack mesh
  overlap = cell.pbc_intor('int1e_ovlp', kpts=kpts)

.. note:: :func:`pbc.Cell.pbc_intor` can only be used to evaluate the short-range
  integrals. PBC density fitting methods have to be used to compute the integrals for
  long-range operators such as nuclear attraction and Coulomb repulsion.

The two-electron Coulomb integrals can be evaluated with the PBC density fitting
methods::

    from pyscf.pbc import df
    eri = df.DF(cell).get_eri()

See :numref:`user_pbc_df` for more details of the PBC density fitting methods.
