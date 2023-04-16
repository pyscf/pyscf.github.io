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
introduced in :ref:`user_gto` also apply here::

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

PySCF uses crystalline Gaussian-type orbitals as basis functions
for solid calculations. 
The predefined basis sets and ECPs for molecular calculations 
can be used in solid calculations as well.
In addition, the predefined basis sets include 
valence basis sets that are optimized for the GTH pseudopotentials 
(a whole list can be found in :source:`pyscf/pbc/gto/basis` and :source:`pyscf/pbc/gto/pseudo`).
The input format of basis sets for the :class:`Cell` object is the same
as that for the :class:`Mole` object.

.. literalinclude:: /../examples/pbc/05-input_pp.py

Finally, custom basis sets can be defined just like in molecular calculations

.. literalinclude:: /../examples/pbc/04-input_basis.py

Low-dimensional systems
-----------------------

The PySCF  :mod:`pbc` module also supports low-dimensional periodic systems. You can initialize
the attribute :attr:`Cell.dimension` to specify the dimension of the system::

  >>> cell.dimension = 2
  >>> cell.a = numpy.eye(3) * 2
  >>> cell.build()

When :attr:`~Cell.dimension` is smaller than 3, a vacuum of infinite size will be
applied in certain direction(s).  For example, when :attr:`~Cell.dimension`
is 2, the z-direction will be treated as infinitely large and the xy-plane
constitutes the periodic surface. When :attr:`~Cell.dimension` is 1, the y and z axes
are treated as vacuum and thus the system is a wire in the x direction. 
When :attr:`~Cell.dimension` is 0, all three directions are treated as vacuum, and this is
equivalent to a molecular system.

K-points
--------

The k-points used in solid calculations can be obtained through the 
:meth:`Cell.make_kpts` method. The minimal input is to specify the number of k-points
in each lattice vector direction::
    
    >>> kpts = cell.make_kpts([1,2,2])
    >>> print(kpts.shape)
    (4,3)

By default, this will return the shifted Monkhorst-Pack mesh which is centered at 
the :math:`Gamma`-point. To get the non-shifted Monkhorst-Pack mesh, one can call::

    >>> kpts = cell.make_kpts([1,2,2], with_gamma_point=False)

To get a shifted Monkhorst-pack mesh centered at a specific point, one can call::

    >>> kpts = cell.make_kpts([1,2,2], scaled_center=[0.,0.25,0.25])

where ``scaled_center`` is defined in the units of lattice vectors.

The obtained k-points are used as input for crystalline electronic structure calculations.
For example, k-point sampled RHF can be invoked as follows::

    >>> from pyscf.pbc import scf
    >>> kpts = cell.make_kpts([2,2,2])
    >>> kmf = scf.KRHF(cell, kpts = kpts)
    >>> kmf.kernel()

More details about k-point sampling for each method can be found in the corresponding chapters.

Spin
----

The attribute :attr:`spin` of :class:`Cell` class has a different meaning to the
:attr:`spin` of :class:`Mole` class. :attr:`Mole.spin` indicates the number of
unpaired electrons of a molecule. Most attributes of class :class:`Cell`
represents the information of a unit cell, :attr:`Cell.spin` is used to indicate
the number of unpaired electrons of a super cell. In a :math:`\Gamma`-point
calculation, one can set `cell.spin` to control the unpaired electrons of the
entire system. For calculations with k-point samplings, for example a 2x2x2
k-point calculation, `cell.spin=1` leads to one unpaired electron rather than 8
unpaired electrons.

If wrong `cell.spin` is set, calculations will still run and a warning
message for inconsistency between the electron number and spin may be raised.

Currently, the program does not support to assign the unpaired electrons per
unit cell. A temporary workaround is to set `cell.spin` to the product of the
number of unpaired electrons and k-points.
However, this setting only guarantees that the total numbers of alpha electrons
and beta electrons agree to the `cell.spin` cross all k-points. The occupancies
may be different at different k-points.


Other parameters
----------------
The attribute :attr:`precision` of :class:`Cell` object determines the integral accuracy.
The default value is ``1e-8`` hartree. When calling the :func:`cell.build()` method,
some parameters are set automatically based on the value of :attr:`precision`, including:

  * :attr:`mesh` -- length-3 list or 1x3 array of int

    - The numbers of grid points in the FFT-mesh in each direction.

  * :attr:`ke_cutoff` -- float

    - Kinetic energy cutoff of the plane waves in FFT-DF

  * :attr:`rcut` -- float

    - Cutoff radius (in Bohr) of the lattice summation in the integral evaluation

Other attributes of the :class:`Mole` class such as :attr:`verbose`,
:attr:`max_memory`, etc., have the same meanings in the :class:`Cell` class.

.. note:: Currently, point group symmetries for crystals are not supported.

Accessing AO integrals
======================

The :func:`Mole.intor` method can only be used to evaluate integrals with open boundary
conditions. When the periodic boundary conditions of crystalline systems are
used, one needs to use the :func:`pbc.Cell.pbc_intor` function to evaluate the
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

See :ref:`user_pbc_df` for more details of the PBC density fitting methods.
