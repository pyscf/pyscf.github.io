# Crystal structure

**Module**: `pyscf.pbc.gto`

**Examples**: [pyscf/examples/pbc](https://github.com/pyscf/pyscf/tree/master/examples/pbc)

## Building a crystal

Periodic crystals are built using the `pbc.gto.Cell` class,
which is very similar to the `gto.Mole` class,
```python
from pyscf.pbc import gto
cell = gto.Cell()
cell.atom = '''H 0 0 0; H 1 1 1'''
cell.basis = 'gth-dzvp'
cell.pseudo = 'gth-lda'
cell.a = numpy.eye(3) * 2
cell.build()
```

The other equivalent ways to build a molecule introduced in [Molecular
structure](gto), including the shortcut functions `pbc.gto.M()` or `pyscf.M()`,
also apply here,
```python
from pyscf.pbc import gto
cell = gto.Cell()
cell.build(
  atom = '''H 0 0 0; H 1 1 1''',
  basis = 'gth-dzvp',
  pseudo = 'gth-lda',
  a = numpy.eye(3) * 2)
```
```python
from pyscf.pbc import gto
cell = gto.M(
  atom = '''H 0 0 0; H 1 1 1''',
  basis = 'gth-dzvp',
  pseudo = 'gth-lda',
  a = numpy.eye(3) * 2)
```
```python
import pyscf
cell = pyscf.M(
  atom = '''H 0 0 0; H 1 1 1''',
  basis = 'gth-dzvp',
  pseudo = 'gth-lda',
  a = numpy.eye(3) * 2)
```

## Geometry and lattice vectors

The `Cell.atom` attribute defines the positions of the atoms inside the unit cell, and
the additional parameter `Cell.a` defines the lattice vectors. 
The format of `Cell.a` is array-like,
```python
cell.a = numpy.eye(3) * 2
cell.a = [[2,0,0],[0,2,0],[0,0,2]]
```

Each row of the 3-by-3 matrix of `Cell.a` represents a lattice vector
in Cartesian coordinates, with the same unit as the input `atom` parameter
(and controllable by the `unit` attribute).

:::{warning}

The input lattice vectors should form a right-handed coordinate system, as
otherwise some integrals may be computed incorrectly in PySCF.
PySCF will print a warning if the lattice vectors do not form a right-handed
coordinate system.

:::

## Basis set and pseudopotentials

PySCF uses crystalline Gaussian-type orbitals as basis functions for periodic 
calculations.  The predefined basis sets and ECPs for molecular calculations
can be used in periodic calculations as well. 

As described more in [Density fitting for crystalline
calculations](pbc/df), many PBC calculations require the use of ECPs (or
pseudopotentials, as they are more commonly called in periodic codes).  In
addition to molecular ECPs, PySCF includes GTH pseudopotentials, which have
been parameterized for use with HF or different DFT functionals,
```python
cell.pseudo = 'gth-hf'
cell.pseudo = 'gth-lda' # an alias for 'gth-pade'
cell.pseudo = 'gth-pbe'
```
The GTH pseudopotentials should always be used with corresponding valence basis sets,
```python
cell.basis = 'gth-szv' # or gth-dzv, gth-dzvp, gth-tzvp
```

Lists of all available GTH pseudopotentials and basis sets are available in
[pbc/gto/pseudo](https://github.com/pyscf/pyscf/tree/master/pyscf/pbc/gto/pseudo)
and [pbc/gto/basis](https://github.com/pyscf/pyscf/tree/master/pyscf/pbc/gto/basis).

:::{note}

GTH basis sets and pseudopotentials are not available for all atoms of the periodic table.

:::

## K-points

The k-points to be used in solid calculations can be obtained through the 
`Cell.make_kpts` method, by specifying the number of k-points
in each lattice vector direction,
```python
kpts = cell.make_kpts([1,2,2])
print(kpts.shape)
# (4,3)
```
By default, this will return the shifted Monkhorst-Pack mesh that includes
the Gamma point. To get the non-shifted Monkhorst-Pack mesh,
```python
kpts = cell.make_kpts([1,2,2], with_gamma_point=False)
```

To get a shifted Monkhorst-pack mesh centered at a specific point,
```python
kpts = cell.make_kpts([1,2,2], scaled_center=[0.,0.25,0.25])
```
where `scaled_center` is defined in the units of reciprocal lattice vectors.

The obtained k-points are used as input for crystalline electronic structure calculations,
```python
from pyscf.pbc import scf
kpts = cell.make_kpts([2,2,2])
kmf = scf.KRHF(cell, kpts=kpts)
e_hf = kmf.kernel()
```

Calculations with k-points always return the energy per unit cell.

## Spin

Because k-point sampling formally represents a calculation on a supercell,
the attribute `Cell.spin` indicates the number of unpaired electrons in the
supercell (**not** in the unit cell). For example, in a calculation with
a 2x2x2 k-point mesh, `cell.spin = 1` leads to one unpaired electron in the
supercell (not eight).

## Low-dimensional systems

The PySCF `pbc` module also supports low-dimensional periodic systems. You can initialize
the attribute `Cell.dimension` to specify the dimension of the system,
```python
cell.dimension = 2
cell.a = numpy.eye(3) * 2
cell.build()
```

When `Cell.dimension` is smaller than 3, a vacuum of infinite size will be
applied in certain direction(s).  For example, when `Cell.dimension = 2`, 
the z-direction will be treated as infinitely large and the xy-plane
constitutes the periodic surface. When `Cell.dimension = 1`, the y and z axes
are treated as vacuum and thus the system is a wire in the x direction. 
When `Cell.dimension = 0`, all three directions are treated as vacuum, and this is
equivalent to a molecular system.

## Other parameters

The `Cell.precision` attribute determines the integral accuracy, and its
default value is `1e-8` hartree. When calling the `cell.build()` method,
some parameters are set automatically based on the value of `precision`, including

  * `mesh` - length-3 list or 1x3 array of int

    - The numbers of grid points in the FFT-mesh in each direction.

  * `ke_cutoff` - float

    - Kinetic energy cutoff of the plane waves in FFTDF

  * `rcut` - float

    - Cutoff radius (in Bohr) of the lattice summation in the integral evaluation

Other attributes of the `Mole` class such as `verbose`,
`max_memory`, etc., have the same meanings in the `Cell` class.

:::{note}

Currently, point group symmetries for crystals are not supported.

:::

## Accessing AO integrals

Periodic AO integrals can be evaluated using the `Cell.pbc_intor` function,
```python
overlap = cell.pbc_intor('int1e_ovlp')
kin = cell.pbc_intor('int1e_kin')
```

By default, the `Cell.pbc_intor` function only returns integrals at the
Gamma point.  If k-points are specified, it will return the integrals at each
k-point,
```python
kpts = cell.make_kpts([2,2,2])
overlap = cell.pbc_intor('int1e_ovlp', kpts=kpts)
```

:::{note}

The `Cell.pbc_intor` function can only be used to evaluate periodic short-range
integrals. PBC density fitting methods have to be used to compute the integrals for
long-range operators such as the electron-nuclear attraction and the electron-electron
repulsion integrals.

:::

The electron repulsion integrals can be evaluated with the periodic density fitting
methods,
```python
from pyscf.pbc import df
eri = df.DF(cell).get_eri()
```
See [Periodic density fitting](pbc/df) for more details.


## Initializing Mean-Field and Post-SCF Methods

Mean-field and post-SCF methods can be initialized directly from the `Cell` instance.

### Instantiating Mean-Field Methods
For single k-point mean-field methods, you can pass the k-point to the `cell.HF`
or `cell.KS` methods.

```python
kpt = np.array([0., 0., 0.])
mf = cell.HF(kpt=kpt)
mf = cell.KS(kpt=kpt)
```

To create mean-field methods with a k-point mesh, the HF and KS methods are
prefixed with "K". The k-points should be created and passed to these K-prefixed
mean-field methods.

```python
kpts = cell.make_kpts([3,3,3])
mf = cell.KHF(kpts=kpts)
mf = cell.KKS(kpts=kpts)
```

Based on the `cell.spin` setting, different mean-field instances are created. By default:
* Restricted closed-shell methods are used for molecules without unpaired electrons (`cell.spin=0`).
* Unrestricted methods (UHF, UKS) are used for molecules with unpaired electrons (`cell.spin > 0`).

The use of k-point symmetry in KHF  and KKS depends on the supplied `kpts`. If
`kpts` are created as a k-point symmetry-enabled (`KPoint`) instance, k-point
symmetry-adapted mean-field methods will be created.

```python
kpts = cell.make_kpts([3,3,3], space_group_symmetry=True)
mf = cell.KHF(kpts=kpts) # KsymAdaptedKRHF
mf = cell.KKS(kpts=kpts) # KsymAdaptedKRKS
```

### Explicit Mean-field Methods
To create specific mean-field instances, such as unrestricted or generalized HF
or KS methods, they can be explicitly called from the `Cell` instance. The
following rules apply:

* `cell.RHF()` and `cell.RKS()` create single k-point restricted Hartree-Fock or
  Kohn-Sham methods. When `cell.spin != 0`, ROHF and ROKS methods are
  instantiated. However, for k-point mean-field methods, `cell.KRHF()` and
  `cell.KRKS()` always create restricted closed-shell models, regardless of the
  `cell.spin` setting. This is different from the treatment in single k-point methods.

* `cell.ROHF()`, `cell.ROKS()`, `cell.KROHF()`, and `cell.KROKS()` create
  restricted open-shell HF or KS methods.

* `cell.UHF()`, `cell.UKS()`, `cell.KUHF()`, and `cell.KUKS()` create
  unrestricted Hartree-Fock or Kohn-Sham methods, regardless of the `cell.spin`
  setting.

* `cell.GHF()`, `cell.GKS()`, `cell.KGHF()` and `cell.KGKS()` create generalized
  mean-field methods.

* K-point symmetry adaptation is available for K-prefixed methods, which is
  determined by the `kpts` argument.

### Keyword Parameters for Mean-field Methods
Configuration parameters for mean-field methods can be specified through keyword
arguments of the instantiation methods. For example,

```python
mf = cell.RKS(xc='b3lyp', max_cycle=20)
mf = cell.KGKS(xc='pbe', kpts=cell.make_kpts([2,2,2]))
```

### Instantiating Post-SCF Methods
Post-SCF methods at the gamma point can be specified for the `Cell` instance.
Keyword arguments can be supplied, and they will be passed to the post-SCF
methods.

```python
mycc = cell.CCSD(frozen=1)
mytd = cell.TDA()
```

For post-SCF methods with k-point sampling, the method should be prefixed with "K".
Additionally, k-points can be supplied as keyword arguments, and they will be
automatically applied when creating the underlying mean-field object.

```python
mycc = cell.KCCSD(kpts=cell.make_kpts([2,2,2]), frozen=1)
```
