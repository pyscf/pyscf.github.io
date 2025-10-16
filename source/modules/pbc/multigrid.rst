.. _pbc_dft_multigrid:

pbc.dft.multigrid --- Multigrid Integration
*******************************************

Multigrid is a numerical integration algorithm optimized for the computation of
the Coulomb matrix and exchange-correlation (XC) functional. It can take
advantage of the locality of density and orbitals when computing electron
density and matrix elements.

The Multigrid algorithm employs different integration grids for different types
of orbitals. Diffuse orbitals use sparse integration grids, while compact
orbitals utilize denser grids. The calculations for compact orbitals are only
performed within smaller spatial regions. The numerical values from different
grids are then aggregated in reciprocal space using the fast Fourier transform.

Compared to the `get_j` function provided by the FFTDF module and the XC matrix
evaluation functions offered by the `numint` module, the Multigrid algorithm
can achieve an order of magnitude improvement in the computation of the Coulomb
matrix and DFT XC matrix.

Supported Features
==================

The Multigrid algorithm is designed to accelerate DFT calculations.
It can also be utilized to speed up derived properties based on DFT, such as
TDDFT excited state calculations, stability analysis, analytical nuclear
gradients, and the computation of the orbital Hessian in second-order SCF
(SOSCF) methods. However, it does not support the computation of the HF exchange
matrix or the calculation of MO two-electron integrals in post-HF methods.

In terms of supported systems, the Multigrid algorithm can be applied to both
periodic systems and molecular calculations. For periodic systems, it supports
single k-point calculations as well as k-point sampling in DFT calculations. For
molecular systems, additional configurations are required to enable the multigrid
algorithm, see :ref:`mole_multigrid`.

Generally, the Multigrid algorithm should be used with pseudopotentials and the
corresponding basis sets to reduce computational costs and accuracy.
Using the Multigrid algorithm for all-electron calculations typically results in significant overhead.


How to enable MultiGrid algorithm
=================================

In periodic system calculations, enabling the Multigrid algorithm is
straightforward. The PBC SCF class provides a method `multigrid_numint()`.
This method creates a new SCF instance which is configured with the Multigrid
functionality. This new instance can then be used in the same way as a standard
SCF instance. For example::

    cell = pyscf.M(
        a = np.eye(3)*3.5668,
        atom = '''C     0.      0.      0.    
                  C     0.8917  0.8917  0.8917
                  C     1.7834  1.7834  0.    
                  C     2.6751  2.6751  0.8917
                  C     1.7834  0.      1.7834
                  C     2.6751  0.8917  2.6751
                  C     0.      1.7834  1.7834
                  C     0.8917  2.6751  2.6751''',
        basis = 'gth-dzv',
        pseudo = 'gth-pade',
    )
    mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts([3,3,3]))
    mf = mf.multigrid_numint()
    mf.run()

Once the Multigrid algorithm is enabled, it is automatically applied to other
methods based on that SCF instance, such as TDDFT and SOSCF. For example::

    mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts([3,3,3]))
    mf = mf.multigrid_numint()

    mf = mf.soscf().run()

    mf.TDA().run()

    mf.Gradients().run()

The Multigrid algorithm is enabled through the `._numint` attribute of the SCF
instance You can directly modify the `._numint` attribute of an SCF instance to
apply the Multigrid algorithm.::

    from pyscf.pbc.dft.numint import MultiGridNumInt
    mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts([3, 3, 3]))
    mf._numint = MultiGridNumInt(cell)
    mf.run()

Configurable Options
====================

In most scenarios, the Multigrid module can automatically configure parameters such
as energy cutoff, radius cutoff, sub tasks to balance computational load and
accuracy. If you need more precise control over the Multigrid algorithm,
you can modify the attributes of the `mf._numint`. Here are some configurable options.

* `mesh`: It's a 3-element tuple, which controls the upper limit of the mesh
  points utilized by the algorithm. By default, this value is estimated based on
  the basis set, the shape, and the size of the lattice. By default, the error in
  the electron density at each mesh point is controlled to be below `cell.precision`.

* `ntasks`: The maximum number of subgroups. The Multigrid algorithm divides
  orbital pairs into subgroups and processes the integration for each subgroup
  separately. Increasing the number of tasks can better utilize the locality of
  the orbital functions, reducing the number of floating-point operations
  (FLOPS). However, it can also increase the overhead of other operations, such
  as FFT. The recommended value for `ntasks` is typically set between 3 and 6.
  The default value is 4.

* `xc_with_j`: Determines whether the calculation of the XC matrix should also
  include the calculation of the Coulomb energy and the Coulomb matrix.
  Combining the calculation of the XC and Coulomb can reduce the computation cost.
  For GGA functionals, it can reduce the computational load by 30% - 40%.
  The default setting is `True`, and it is generally recommended to enable this option.

* `ke_ratio`: Controls the ratio of the energy cutoffs among the different
  subtasks when partitioning the orbital pairs. `ke_ratio` is a number greater
  than 1. When approaching 1, the task partitioning is fine, and the orbital
  locality is better utilized. However, this also increases the number of tasks.
  Regardless of the configuration of this parameter, the maximum number of tasks
  generated is still controlled by the `ntasks` attribute.


Compatibility with NumInt
=========================

The `NumInt` class offers numerical integration functions for DFT calculations.
It is assigned to the `_numint` attribute of SCF classes. This module offers
APIs such as `get_rho` (for computing electron density in real space) and
`get_vxc` (for computing the value of the XC energy functional and XC matrices),
for DFT calculations :ref:`dft`. Both molecular and PBC DFT implementations
follow this API design.

The `pyscf.pbc.dft.multigrid` module offers the `MultiGridNumInt` class, which
is compatible with the `NumInt` class. Its methods, such as `get_vxc`, `get_fxc`,
`get_rho`, `cache_xc_kernel`, `nr_rks`, and `nr_uks`, are mostly compatible
with the corresponding methods in NumInt (with only a few additional keyword
arguments for controlling multigrid instances). These methods can be
individually invoked, like those in the `NumInt` class, to compute densities and
XC matrix elements.

The `pyscf.pbc.dft.multigrid` module also provides the `MultiGridNumInt2` class,
which further optimizes the implementations of the `MultiGridNumInt` class.
However, due to differences in the algorithm implementations, the support for
optimized k-points and non-orthogonal lattices is not as comprehensive as that
in the `MultiGridNumInt` class. Currently, the `SCF.multigrid_numint()` method
invokes the `MultiGridNumInt` class. To maximize the multigrid performance, you
can manually assign the `MultiGridNumInt2` instance to the `mf._numint`
attribute.

The two classes will be merged into one in the future release.


.. _mole_multgird:

How to Apply Multigrid in Molecular DFT
=======================================

The Multigrid algorithm is currently designed and implemented for data
structures with periodic boundary conditions. Nevertheless, it can be adapted
for molecular calculations.

First, we need to initialize molecule within a relatively large periodic lattice.
A vacuum space needs to be placed between the box and the molecule to
simulate free boundary conditions::

    cell = pyscf.M(atom='''
    O    0.000    0.118  0.
    H    0.758   -0.470  0.
    H   -0.758   -0.470  0.''',
    a=np.eye(3)*10,
    dimension=0,
    basis='gth-dzvp', pseudo='gth-pade',)

Here, we apply a 10 x 10 x 10 (unit Angstrom) box. The box does not need to be
excessively large. It only needs to ensure that the electron density of the molecule
does not leak outside the box. If there are no diffused functions, typically, a
5 Angstrom margin around the molecule is sufficient. The lattice is a virtual
box which does not need to be centered on the origin. There is no need to adjust
the molecule's coordinates to the center of the box.

Alternatively, the lattice can be automatically determined. We can just create a
`Mole` instance and then call the `Mole.to_cell()` method to convert a `Mole`
object into a Cell object::

    mol = pyscf.M(atom='''
    O    0.000    0.118  0.
    H    0.758   -0.470  0.
    H   -0.758   -0.470  0.''',
    basis='gth-dzvp', pseudo='gth-pade')
    cell = mol.to_cell(dimension=0)

When initializing the Cell, we set `dimension=0`. This setting informs the
program that the system is a 0-dimensional system (molecule), which allows the
program to more effectively utilize this property to truncate the long-range
Coulomb potential and accelerate the computation of certain integrals.
The system can also be treated as a 3-dimensional crystal, with slightly
increased computational load.

In this system, we use the GTH pseudo potential and the GTH basis set. This
basis set does not have very steep GTO functions, making it suitable for the
multigrid method.

With this `Cell` object, we can initialize the DFT instance as we would for
typical PBC DFT calculations. The following example demonstrates another
way which integrates
the Multigrid XC matrix functionality into the the molecular DFT instances.::

    from pyscf.pbc.dft.multigrid import MultiGridNumInt2
    mol = pyscf.M(atom='''
    O    0.000    0.118  0.
    H    0.758   -0.470  0.
    H   -0.758   -0.470  0.''',
    basis='gth-dzvp', pseudo='gth-pade')
    mf = mol.RKS(xc='pbe')

    cell = mol.to_cell(dimension=0)
    mf._numint = MultiGridNumInt2(cell)
    mf._numint.xc_with_j = False

    mf.run()

In this example, we use the `MultiGridNumInt` class for `mf._numint` of the
molecular DFT instance. This setup invokes the `MultiGridNumInt` algorithm to
calculate the DFT XC matrix and XC energy, while calls the standard molecular
`get_jk` functions for Coulomb energy.
Note the setting `xc_with_j=False`, which disables the computation of Coulomb
energy by the `MultiGridNumInt` class. This is necessary because the molecular
DFT program employs a different method to compute nuclear repulsion energy
compared to the method used in PBC DFT.
When using the molecular DFT program in conjunction with `MultiGridNumInt`, we
should not use `MultiGridNumInt` to compute the Coulomb energy.


Limitations of Multigrid algorithm
==================================

* The Multigrid algorithm only supports uniform grids. Currently, it cannot be used with Becke grids.

* The `MultiGridNumInt` class does not support the calculation of analytical nuclear gradients.

* The `MultiGridNumInt2` class does not support k-point
  calculations, meta-GGA functionals, or the computation of the TDDFT fxc kernel.


Examples
========
* :source:`examples/pbc/27-multigrid.py`
* :source:`examples/pbc/27-multigrid2.py`
