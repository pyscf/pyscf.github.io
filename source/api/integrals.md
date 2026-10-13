# Integral and density-fitting conventions

This page describes integral ordering, packed storage, AO-to-MO transformations,
and molecular density-fitting factors.

Relativistic and periodic integral interfaces are currently not covered.

## Chemist notation and packed integral storage

Unless stated otherwise, molecular two-electron integrals in PySCF use
chemist notation:

$$
(ij|kl) = \int\!\int
\phi_i^*(\mathbf r_1)\phi_j(\mathbf r_1)
\frac{1}{r_{12}}
\phi_k^*(\mathbf r_2)\phi_l(\mathbf r_2)
\,d\mathbf r_1\,d\mathbf r_2.
$$

For real integrals over `n` orbitals, let `npair = n * (n + 1) // 2`.
The common storage formats are:

| Format | Shape | What is stored |
| --- | --- | --- |
| `s1` | `(n, n, n, n)` | All four indices |
| `s4` | `(npair, npair)` | A matrix indexed by the pairs `(i,j)` and `(k,l)` |
| `s8` | `(npair * (npair + 1) // 2,)` | The lower triangle of the `s4` pair matrix |

Within a packed pair, indices follow `(0,0), (1,0), (1,1), (2,0), ...`.
For `i >= j`, the pair index is `i * (i + 1) // 2 + j`. There is no extra
factor of two or square-root scaling on off-diagonal entries.

Use `ao2mo.restore()` to unpack or repack real `float64` integrals:

```python
from pyscf import ao2mo

eri_full = ao2mo.restore(1, eri_packed, n)
eri_fourfold = ao2mo.restore(4, eri_full, n)
```

When packing an integral tensor, `restore()` does not average inconsistent
entries to the tensor symmetric.

Four-fold storage assumes symmetry within each pair:

$$
(ij|kl) = (ji|kl) = (ij|lk) = (ji|lk).
$$

It does **not** require `(ij|kl) = (kl|ij)`. An `s4` pair matrix can therefore
be nonsymmetric. Eight-fold storage requires this additional pair-exchange
symmetry.

## AO-to-MO output conventions

The order of the four coefficient matrices in
`ao2mo.general(eri, (C1, C2, C3, C4))` follows the four integral indices:

$$
(ij|kl)_{\mathrm{MO}} = \sum_{pqrs}
C_{1,pi}^* C_{2,qj} C_{3,rk}^* C_{4,sl} (pq|rs)_{\mathrm{AO}}.
$$

`aosym` describes the symmetry of the input AO integrals. `compact` controls
whether pairs of MO indices are packed in the output when their symmetry
allows it.

For packed AO input, the output is a two-dimensional array: the first two
MO indices form the row index, and the last two form the column index.
With `compact=False`, its shape is `(ni * nj, nk * nl)`, where `ni`, `nj`,
`nk`, and `nl` are the numbers of columns in `C1`, `C2`, `C3`, and `C4`.
This result can be reshaped to `(ni, nj, nk, nl)` in NumPy C order.

When `compact=True`, if the first and second coefficient matrices are real
and identical, their pair index keeps only the lower-triangle entries.
The same applies to the third and fourth matrices. One pair can be packed
while the other remains unpacked.

For the standard molecular file-based transformation, the HDF5 dataset
uses the two-dimensional pair layout described above.


Integrals with multiple operator components have an additional leading component
axis.

## Complex-integral support
TODO


## Mixed Cartesian and spherical one-electron overlaps

`gto.intor_cross('int1e_ovlp', mol1, mol2)` follows each molecule's `cart`
setting. Rows belong to `mol1` and columns to `mol2`. If `mol1` is Cartesian
and `mol2` is spherical, the shape is `(ncart1, nsph2)`; the reverse case
gives `(nsph1, ncart2)`.

An explicit suffix overrides this choice. For example, `'int1e_ovlp_cart'`
returns Cartesian functions on both sides, even if one molecule has
`cart=False`. Likewise, `'int1e_ovlp_sph'` requests spherical functions on
both sides.
