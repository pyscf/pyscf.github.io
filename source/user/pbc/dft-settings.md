# Guidelines for PBC DFT Parameter Settings

DFT calculations for extended systems with PBC involve multiple parameters
to control how integrals are evaluated. Compared to molecular DFT calculations,
these settings are significantly more complex. In molecular DFT calculations,
integral-related settings mainly affect computational efficiency, with only
minor impact on accuracy. However, in periodic DFT, integral parameters are
coupled to the choice of pseudopotential, the type of exchange-correlation
functional, basis sets, boundary conditions, and the available memory and disk
space.

To support various computation scenarios, PySCF provides multiple algorithms for
computing Coulomb integrals and the numerical integration options for DFT XC
functional. This document discusses how to select the suitable schemes for
different use cases.

## Overview of the Available Integral Schemes
### Coulomb Integral Schemes

* **FFTDF** (aka the GPW algorithm):
```python
from pyscf.pbc.df import FFTDF
kpts = cell.make_kpts(kmesh)
mf = cell.KRKS(xc='pbe0', kpts=kpts)
mf.with_df = FFTDF(cell, kpts=kpts) # Unecessary, this is the default
mf.run()
```
  This algorithm employs plane-wave (PW) functions as auxiliary functions in the
  density fitting method. Coulomb integrals are computed using real-space
  density evaluation with fast-Fourier transform. In this approach, electron
  density is evaluated on real space grids and the Poisson equation is solved in
  reciprocal space to evaluate Coulomb integrals. The fast Fourier transforms
  (FFTs) are used to transform quantities between real space and reciprocal
  space.

  With a sufficiently high energy cutoff in reciprocal space (approaching the
  complete basis limit for plane waves), the FFTDF method does not introduce any
  approximations. This algorithm primarily serves as the reference
  implementation for other integral algorithms. Various features, such as the
  nuclear gradients, stress tensor, are first developed and validated with FFTDF.

  However, FFTDF is generally considered computationally inefficient and
  intensive in memory consumption. It is not suitable for treating core
  electrons or evaluating exact exchange in hybrid or HF methods. Nevertheless,
  due to its formally exact, FFTDF remains the default integral scheme for DFT,
  HF, and post-HF calculations.

* **AFTDF** (Analytical Fourier Transform Density Fitting):
```python
from pyscf.pbc.df import AFTDF
kpts = cell.make_kpts(kmesh)
mf = cell.KRKS(xc='pbe0', kpts=kpts)
mf.with_df = AFTDF(cell, kpts=kpts)
mf.run()
```
  AFTDF extends the FFTDF approach by analytically evaluating the Fourier
  transforms of Gaussian basis functions.
  Compared to FFTDF, when using the same PW energy cutoff, AFTDF introduces
  smaller errors in Coulomb integrals, also consumes significantly less memory.
  However, AFTDF is generally slower than FFTDF. Similar to the FFTDF, AFTDF is
  generally not suitable for describing core electrons and exact exchange, due
  to the requirement of a high energy cutoff for PWs and the large number of
  auxiliary plane waves required.

  A unique feature of AFTDF, however, is its ability to perform Fourier
  transforms for non-uniform PW grids. This property allows AFTDF to efficiently
  handle low-dimensional systems (like surfaces, and wires) by introducing
  infinite vacuum regions and using non-uniform grids along non-periodic
  directions.

* **GDF** (Gaussian density fitting):
```python
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh)).density_fit(auxbasis='weigend')
mf.run()
```
  The GDF (Gaussian Density Fitting) algorithm uses Gaussian auxiliary basis
  functions to represent the product of atomic orbitals. Unlike the plane-wave
  auxiliary functions used in FFTDF or AFTDF, Gaussian auxiliary functions can
  efficiently describe both core and valence electrons with a small number of
  basis functions. This property makes GDF effective for all-electron
  calculations and hybrid functional DFT computations.

  GDF method has a significant overhead to construct the three-index tensor
  during initialization. Once the this tensor is generated, it is very efficient
  to evaluate the HF exchange matrix in GDF due to the small number of auxiliary
  functions. The available disk space and memory will limit the problem size.
  GDF may become impractical for very dense k-point sampling.

  In the PySCF CPU version, the GDF three-index tensor is stored on disk. In
  GPU4PySCF, this tensor is compressed and held completely in host memory.
  The compression reduces the size of the tensor to approximately 20% of the
  original PySCF size. However, the compression storage in GPU4PySCF requires
  additional work to decompress the tensor, which can become more expensive than
  evaluating HF exchange when many k-points are involved.

  GDF supports all types of periodicity, including fully 3D periodic systems
  and low-dimensional systems. It also allows the evaluation of four-index
  Coulomb integrals in both AO and MO representations. Other features, like the
  computation of nuclear gradients, are currently limited.

  Compared to the auxiliary PW fitting functions (FFTDF and AFTDF),
  GDF can produce larger errors in Coulomb integrals due to the incompleteness
  of the Gaussian fitting functions. Diffuse fitting functions can cause linear
  dependencies, which further increases the error. These errors are typically
  more significant than those observed in standard molecular DF calculations.

  Another known limitation of the current GDF implementation is its handling
  of k-point band calculations. Even a small number of k-points can trigger
  the generation of many three-index tensors, which often leads to huge
  memory and disk usage.

  Despite these limitations, in various scenarios, GDF remains a practical and
  preferred option for hybrid functional DFT and post-HF calculations. GDF often
  provides the best trade-off between accuracy and computational efficiency for
  all-electron and hybrid functional calculations.

  See [PBC denisty fitting](pbc_df) for more details of the GDF,
  FFTDF, and AFTDF modules.

* **RSDF** (Range-separated density fitting):
```python
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh)).rs_density_fit(auxbasis='weigend')
mf.run()
```
  RSDF is an alternative implementation of the GDF method. RSDF is designed to
  produce the same results as the GDF program. The implementation of RSDF is
  generally more efficient than GDF; however, it may introduce slightly larger
  errors. The GDF implementation is more conservative in terms of integral
  screening and various cutoff thresholds.

  The functionalities of RSDF are less comprehensive than the GDF module. It
  supports 3D calculations but does not support low-dimensional calculations.

* **MDF** Mixed density fitting):
```python
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh)).mix_density_fit(auxbasis='weigend')
mf.with_df.mesh = [7, 7, 7]
mf.run()
```
  MDF combines PWs and Gaussian functions as auxiliary basis functions. In
  principle, MDF should be more accurate than GDF, while computationally being
  cheaper than FFTDF and AFTDF.

  However, in practice, MDF often suffers from severe linear dependency between
  the plane-wave and Gaussian components of the fitting basis. Attempts to
  remove these linear dependencies can introduce additional numerical errors,
  which in some cases make MDF less accurate than GDF. From a computational
  perspective, MDF is significantly slower than GDF.

  As a result, MDF is typically not recommended. It is only considered in
  scenarios where GDF errors are unacceptably large and alternative methods such
  as RSJK are too slow.

* **ISDF** (Interpolative separable density fitting):
  This method is under development. This method is efficient for hybrid DFT.

* **RSJK** (Range-separated J/K builder):
```python
from pyscf.pbc.scf.rsjk import RangeSeparatedJKBuilder
kpts = cell.make_kpts(kmesh)
mf = cell.KRKS(xc='pbe', kpts=kpts)
mf.rsjk = RangeSeparatedJKBuilder(cell, kpts=kpts)
mf.run()
```
  The RSJK method evaluates Coulomb integrals (the short-range part) in real
  space, with the long-range integrals augmented by AFTDF. RSJK is an accurate
  integral evaluation algorithm. It can produce results identical to FFTDF (or
  AFTDF) with an infinite high energy cutoff.

  RSJK is suitable for all-electron and hybrid functional DFT calculations. It
  also supports 3D periodic systems and low-dimensional systems.

  RSJK is efficient for the evaluation of HF exchange over a large number
  k-points. However, it is not as efficient for a small number of k-points. Its
  scaling against the number of k-points is low. For hybrid DFT calculations
  with a 5x5x5 k-mesh or larger, RSJK can become more efficient than any of the
  previously mentioned density fitting methods.

### Numerical Integration Grids
In periodic DFT calculations, XC integrals are evaluated numerically on a set of
real-space mesh grids. Two primary types of grids are commonly used: uniform
grids and atomic-centered grids.

* **UniformGrids**:
```python
from pyscf.pbc.dft import UniformGrids
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh))
mf.grids = UniformGrids(cell)
mf.run()
```
  Uniform grids consist of an evenly spaced mesh that spans the unit cell. These
  mesh grids can be shared with the mesh used in the FFTDF-based Coulomb integral
  computation code. The XC and Coulomb computations can be combined to reduce the
  overall computational cost.

  Although a sufficient dense uniform mesh can accurately describe the core electron
  density, it require significantly more grids than atomic-centered grids. As a
  result, uniform grids are usually employed when pseudopotentials are used.

* **BeckeGrids**:
```python
from pyscf.pbc.dft import BeckeGrids
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh))
mf.grids = BeckeGrids(cell)
mf.run()
```
  Becke grids are classical atom-centered numerical grids adapted for periodic
  boundary conditions. These grids are suitable for all-electron calculations.
  The number of grid points in periodic BeckeGrids is significantly larger than
  those in molecular DFT. Despite this, the total number of Becke grids required
  is still substantially smaller than that needed for equivalent accuracy with
  uniform grids.

### Numerical Integration Schemes
Each type of DFT integration grid is associated with a corresponding
XC integration algorithm. Two algorithms are implemented in PySCF:

* **NumInt**:
```python
from pyscf.pbc.dft import NumInt, KNumInt
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh))
# mf.grids = KNumInt() # Unecessary, this is the default
mf.run()
```
  The NumInt module follows the implementation used in molecular DFT. It
  supports both UniformGrids and BeckeGrids. The computational cost of NumInt is
  roughly proportional to the number of grid points. The locality of orbitals
  and the sparsity of the XC matrix are not explored in the NumInt code. As a
  result, the cost scales roughly cubicly with system size. Despite this, NumInt
  remains the default XC integration scheme as it supports most functionalities
  (such as the fxc kernel integration in TDDFT).
  
  For uniform grids, NumInt can often be replaced by the more efficient
  MultiGrid algorithm.

* **MultiGridNumInt**:
```python
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh))
mf = mf.multigrid_numint()
mf.run()
```
  The MultiGrid algorithm is a numerical integration scheme optimized for Coulomb
  and XC functional evaluation. It takes advantage of the locality of orbitals,
  employing different integration grids for different orbital types.
  
  MultiGrid is designed to work only with UniformGrids and does not support
  BeckeGrids. It is commonly used for calculations with pseudopotentials.

  More details of Multigrid integration code are documented in [Multigrid Integration](pbc_dft_multigrid)
  
## Scenarios and Suggested Configurations

### Summary of Scenarios and Integration Options

The following table summarizes common PBC DFT scenarios and the available
integral or XC integration algorithms for each case:

| Scenario                              | Available Algorithms |
|---------------------------------------|----------------------|
| All-electron calculations             | GDF, RSDF, RSJK, NumInt-BeckeGrids |
| Pseudopotentials                      | FFTDF, MultiGrid, NumInt |
| Semi-local functionals (w/o HFX)      | FFTDF, GDF, RSDF, MultiGrid, NumInt-BeckeGrids |
| Hybrid DFT                            | GDF, RSDF, RSJK |
| 2D systems with truncated Coulomb     | FFTDF, GDF, RSJK, MultiGrid, NumInt-BeckeGrids |
| 1D systems with infinite vacuum       | AFTDF, GDF, RSJK, NumInt-BeckeGrids |

### Typical Scenarios

Several key dimensions influence the choice of integral algorithms and numerical
configurations in the PBC DFT calculations:

* Functional type: Semi-local (without HFX) vs. hybrid (with HFX)
* Core treatment: All-electron vs. pseudopotential
* Periodicity: Fully 3D periodic vs. low-dimensional systems (1D/2D slab models)

In the following, we analyze each combination of these dimensions, and provide
practical guidance for setting up efficient integral schemes and configurations.

* Pseudopotential + Semi-local XC + 3D Periodicity

This is one of the most common scenarios in PBC DFT calculations. The
default settings of the KRKS and KUKS classes can directly handle this type of
system. However, the default `NumInt` calculator is not computationally
efficient for large systems. For better performance, it is recommended to use
the `MultiGridNumInt` calculator. For example
```python
cell = pyscf.M(..., pseudo='gth-pbe', basis='gth-dzvp')
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh))
mf = mf.multigrid_numint()
mf.run()
```

* Pseudopotential + Hybrid functional + 3D Periodicity

Gaussian-based PBC DFT is considered more affordable for hybrid functional
calculations than plane-wave DFT program. This feature makes hybrid functional
also a common scenario in PBC DFT simulations. 

In this setting, the default FFTDF algorithm can be used
for small systems. It is inefficient for larger unit cells or
dense k-point meshes due to its high memory and computational requirements.

For larger systems, it is recommended to switch to more efficient integral
algorithms for evaluating the HF exchange integrals. The available options
include **GDF**, **RSDF**, and **RSJK**.
The size of the three-index integral tensors generated by GDF and RSDF scales
with the square of the number of k-points and the cube of the unit-cell size. In
practice, these methods are suitable when the product $N_{AO} * N_{k}$ is below
approximately 15,000. Beyond this limit, the RSJK algorithm is more appropriate
due to its lower memory requirements.

For the semi-local XC functional part of the hybrid DFT, the default
`NumInt` calculator can be safely used. It can also be replaced by the
`MultiGridNumInt` calculator. Since the XC integration typically takes a small
fraction of the total computational cost in hybrid DFT, the choice between
`NumInt` and `MultiGridNumInt` has little impact on overall performance.

```python
cell = pyscf.M(..., pseudo='gth-pbe', basis='gth-dzvp')
mf = cell.KRKS(xc='hse06', kpts=cell.make_kpts(kmesh)).density_fit()
mf = mf.multigrid_numint()
mf.run()
```

  To perform RSJK algorithm, we can execute
```python
from pyscf.pbc.scf.rsjk import RangeSeparatedJKBuilder
kpts = cell.make_kpts(kmesh)
mf = cell.KRKS(xc='pbe', kpts=kpts)
mf = mf.multigrid_numint()
mf.rsjk = RangeSeparatedJKBuilder(cell, kpts=kpts)
mf.run()
```
The order of applying `multigrid_numint()` and setting `mf.rsjk = ...`, does not
affect the final setup.


* Pseudopotential + Semi-local XC + Low-dimensional System (2D and 1D)

The setup for 2D calculations is almost the same as that for the 3D calculations.
To perform a 2D calculation, we can simply set
`cell.dimension = 2`. This will apply a truncated Coulomb potential along
the z-axis. Note, a reasonable vacuum is still required in the
z-direction as that in the 3D case. Both the `NumInt` and `MultiGridNumInt`
calculators fully support 2D calculations. The `MultiGridNumInt` is more
favorable as it is more efficient than `NumInt`.
```python
cell = pyscf.M(..., dimension=2, a=np.diag([3.6, 3.6, 15]))
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh))
mf = mf.multigrid_numint()
mf.run()
```

For 1D (wire) systems, a practical and commonly used approach is to
perform the calculation in 3D or 2D with a sufficiently large vacuum
along the non-periodic directions. For example, we can use the 2D setup with a
large vacuum along the y-axis. However, this setting may slightly break the
symmetry between the y-axis and z-axis.
```python
cell = pyscf.M(..., dimension=2, a=np.diag([3.6, 15, 15]))
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh))
mf = mf.multigrid_numint()
mf.run()
```

PySCF also offers methods for infinite vacuum boundaries in low-dimensional
systems. The `FFTDF`, `NumInt`, and `MultiGridNumInt` algorithms cannot be used
in this mode. In the infinite vacuum mode, achieving similar accuracy
requires a significantly larger computational cost than the above configurations
that use a truncated Coulomb potential with `cell.dimension=2`.
Therefore, we will not discuss the integral configurations under the infinite
vacuum mode here.

* Pseudopotential + Hybrid functional + Low-dimensional System (2D and 1D)

If the truncated Coulomb potential with `cell.dimension=2` is applied, the GDF
and RSJK algorithms can still be used for 2D calculations with hybrid
functionals. However, the RSDF scheme does not support the truncated Coulomb
potential. Settings for XC integration are the same as those for the semi-local
XC functional scenario mentioned above.

* All-electron + Semi-local XC + 3D Periodicity

When all-electron basis sets are used, the default FFTDF algorithm becomes
inappropriate because of the extremely high PW energy cutoff. Even if HF
exchange is **not** computed, it is still necessary to employ one of the GDF,
RSDF, or RSJK methods to evaluate the Coulomb matrix efficiently and accurately.

Due to the limitations in the implementation of GDF and RSDF, we cannot use a
very large k-point mesh when employing these integral methods.

For the XC numerical integration, the default `NumInt` with uniform grids is
also not suitable for all-electron calculations. The `MultiGridNumInt` algorithm
does not support all-electron calculations. In this case, one has to use
the `NumInt` calculator with atomic `BeckeGrids`. As a result, all-electron +
semi-local XC calculations are much slower than their corresponding
pseudopotential + semi-local XC counterparts.
```python
mf = cell.KRKS(xc='pbe', kpts=cell.make_kpts(kmesh)).density_fit()
mf.run()
```

* All-electron + Hybrid functional + 3D Periodicity

The settings required for all-electron calculations with hybrid functionals are
the same as those for all-electron calculations with semi-local XC functionals.
The GDF, RSDF, or RSJK integral schemes should be used for Coulomb and HF
exchange evaluations, together with NumInt-BeckeGrids for XC integration.

Although hybrid functional calculations are still slower than semi-local XC
calculations, the difference in computational cost is moderate compared to the
pseudopotential-based hybrid DFT. The all-electron hybrid DFT is about 3-5 times
slower than semi-local XC during the initialization of the GDF three-index
tensor. This slowdown is notably smaller than the performance gap in
corresponding pseudopotential calculations, where hybrid functionals can be
orders of magnitude more expensive.
```python
mf = cell.KRKS(xc='pbe0', kpts=cell.make_kpts(kmesh)).density_fit()
mf.run()
```

* All-electron + Low-dimensional Systems

The setup for all-electron low-dimensional systems is similar to that for 3D
systems. For Coulomb integrals, one should use GDF or RSJK, regardless of
whether HF exchange is required. The RSDF algorithm is not suitable in this
scenario because it does not account for the truncated Coulomb potential
(`cell.dimension=2`).

For XC integration, the only suitable option is NumInt with BeckeGrids.

```python
cell = pyscf.M(..., dimension=2, a=np.diag([3.6, 3.6, 15]))
mf = cell.KRKS(xc='pbe0', kpts=cell.make_kpts(kmesh)).density_fit()
mf.run()
```
