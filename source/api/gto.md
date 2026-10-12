# Mole and GTO API

The `pyscf.gto.Mole` class holds molecular data and settings used to initialize
calculation objects, including `stdout`, `verbose`, and `max_memory`. It also
prepares the `_atm`, `_bas`, and `_env` arrays used by the integral libraries.

## Units

Input geometry uses Angstrom by default. Set `mol.unit = 'Bohr'` to supply
coordinates in atomic units. **`mol.unit` controls geometry input only**.
It does not change the units used internally or for other inputs and results.
For example, internal coordinates are always in Bohr.
`atom_coord(i)` and `atom_coords()` return Bohr by default.

Coordinates passed to `eval_gto()` or `dft.numint.eval_ao()` must be in
Bohr, regardless of `mol.unit`.

`mol.unit` applies when geometry input is read, including when building a
molecule or updating its geometry. Changing `mol.unit` alone after
`mol.build()` does not convert stored coordinates or change the units of
energies, gradients, or Hessians.

A numeric `mol.unit` divides input coordinates by that value to obtain Bohr.
Thus `unit=1` means Bohr, while `unit=0.52917721092` interprets coordinates as
Angstrom using that Bohr-to-Angstrom factor. A numeric value can supply a
conversion constant different from `pyscf.data.nist.BOHR`.

## Building a molecule

`build()` translates the input attributes into two representations:

| Input | Internal Python representation | Integral-library representation |
| --- | --- | --- |
| `atom` | `_atom` | `_atm` and coordinates in `_env` |
| `basis` | `_basis` | `_bas` and exponents/coefficients in `_env` |
| `ecp` | `_ecp` | `_ecpbas`, effective charges in `_atm`, and `_env` |
| `pseudo` | `_pseudo` | Effective charges in `_atm`; projector data remain in `_pseudo` |

After modifying input attributes, call `build()` to regenerate the internal
data from those inputs. Changes to internal arrays are not propagated back to
the input attributes and can be overwritten by `build()`. See
[Internal Python data formats](#internal-python-data-formats) for the structures
of `_atom`, `_basis`, `_ecp`, and `_pseudo`.

`gto.M(**kwargs)` and `pyscf.M(**kwargs)` create a `Mole` and call `build()`.
Most `set_geom_()` updates and `Mole.to_cell()` also rebuild. These operations
can replace inconsistent internal data. `Cell.to_mol()` reuses the cell's
built arrays and reconstructs molecular symmetry data when needed.

## Shell ordering

`build()` groups `_bas` by atom in geometry order. For each atom type,
`format_basis()` sorts `_basis` shells by increasing angular momentum; shells
with the same angular momentum retain their input order. `_bas` follows this
shell order for every atom using that basis.

Primitive order is different: `_basis` retains the primitive rows supplied by
the parser, while `build()` stores exponents in descending order in `_env` and
reorders their contraction coefficients with them. Use `bas_exp()` and
`bas_ctr_coeff()` to inspect the built representation.

The ordering is part of the internal contract. `aoslice_by_atom()` requires
an atom's shells to be contiguous, and `decontract_basis()` requires equal
angular momenta to be grouped. AO offsets, labels, symmetry-adapted bases, and
derivative methods use the same shell-to-AO mapping. Do not reorder `_basis`,
`_bas`, or their `_env` data independently.

Code that constructs a fictitious `Mole`, such as a pseudopotential integral
builder, must preserve these relationships when using standard `Mole`
accessors.

## AO ordering

Within each shell, AO ordering follows contraction, then angular component.

- Real spherical AOs (`mol.cart = False`) have `2*l+1` components per
  contraction: `x, y, z` for p; `xy, yz, z^2, xz, x2-y2` for d; and labels
  `-l` through `+l` for higher angular momenta.
- Cartesian AOs (`mol.cart = True`) have `(l+1)*(l+2)/2` components per
  contraction, ordered by decreasing powers of x and then y. The d order is
  `xx, xy, xz, yy, yz, zz`.
- The j-adapted spinor basis groups components by `j`, then increasing `m_j`,
  as detailed in [Spinor ordering and time reversal](#spinor-ordering-and-time-reversal).

`ao_labels()` returns spherical or Cartesian labels according to `mol.cart`;
`spinor_labels()` returns spinor labels.

## Cartesian AO normalization

With `mol.cart = True`, Cartesian AOs in PySCF are not necessarily normalized
to one. The diagonal entries of the AO overlap matrix can therefore differ
from one.

## Spinor ordering and time reversal

For each contraction with angular momentum `l`, a basis entry with `kappa=0`
contains the `j=l-1/2` block followed by the `j=l+1/2` block. The first block
is absent for s functions. Positive `kappa` selects only `j=l-1/2`, and
negative `kappa` selects only `j=l+1/2`. Within each block, `m_j` increases
from `-j` to `+j`. `spinor_labels()` exposes these labels.

`mol.time_reversal_map()` returns one signed integer per spinor AO. Unlike
ordinary AO indices, its absolute values are **one-based**. For zero-based
index `i`, let `t = mol.time_reversal_map()[i]`. Then

```text
T |i> = sign(t) |abs(t) - 1>
```

The sign carries the basis phase. Time reversal also conjugates coefficients
of linear combinations; applying it twice gives the negative of the original
spinor. The map uses spinor AO ordering.

## AO offsets

`mol.ao_loc` (equivalently `ao_loc_nr()`) gives spatial AO offsets according
to `mol.cart`; `ao_loc_2c()` gives spinor offsets. Each array has length
`mol.nbas + 1`. Shell `i` occupies the slice `ao_loc[i]:ao_loc[i+1]`;
the last entry is the total AO count.

`aoslice_by_atom()` returns an integer array of shape `(mol.natm, 4)`.
Each row is `(shell_start, shell_stop, ao_start, ao_stop)` for one atom, with
zero-based, exclusive stop indices suitable for Python slices. Pass
`ao_loc=mol.ao_loc_2c()` to obtain spinor AO boundaries.

## AO label searches

`search_ao_label()` returns indices in the same AO order as `ao_labels()`
and the rows of `mo_coeff`. Numbered atom labels must match the atom index
exactly. For example, `mol.search_ao_label('0 H')` selects hydrogen AOs on
atom 0; it must not also select atoms 10 or 20. Atom indices are zero-based
by default.

General patterns such as `'H.*s'` can still use regular expressions. Exact
matching of the atom number in a numbered label does not make all searches
literal-string searches. The earlier behavior where `'0 H'` also matched
`'10 H'` was a bug, not a search convention.

## Updating geometry with `set_geom_`

`set_geom_(atoms_or_coords, unit=None, symmetry=None, inplace=True)` accepts
an `atom` string/list or a NumPy coordinate array of shape `(mol.natm, 3)`
for the existing atoms. It returns the updated molecule.

`unit=None` uses `mol.unit`. Supplying a different conversion factor updates
the returned molecule's unit. This setting describes the units of the input
coordinates; the updated internal coordinates are still in Bohr. When passing
the default output of `mol.atom_coords()` back to `set_geom_()`, explicitly
specify `unit='Bohr'` so the coordinates are interpreted correctly.

The update follows two paths:

- A NumPy array with symmetry disabled updates `atom`, `_atom`, and coordinates
  in `_env`, and clears the cached nuclear repulsion energy.

- Otherwise, `build()` refreshes the internal data.

`symmetry=None` uses the current setting. To disable symmetry on a symmetric
molecule, set `mol.symmetry = False` before the update. Passing
`symmetry=False` with a NumPy array selects the direct-update path but does not
clear the molecule's existing symmetry attributes.

`inplace=True` updates the existing molecule. `inplace=False` returns a shallow
copy with a separate `_env`; other nested data may remain shared. Neither
resets calculation objects. Use a scanner or `method.reset(new_mol)` before
running calculations at the new geometry.

## Internal Python data formats

These formats are useful for custom basis sets and interfaces. Prefer parser
and accessor functions to direct mutation of built internal attributes.

### Atoms

`mol._atom` is an ordered list of `(symbol, [x, y, z])` entries in Bohr.
Symbols are standardized but retain atom labels such as `H1` or `H@2`.
`atom_symbol(i)` retains the label; `atom_pure_symbol(i)` removes it.
Ghost atoms retain basis functions but have zero nuclear charge.

`gto.format_atom()` returns this representation. Set `unit='Bohr'` when using
it as `mol.atom` input. `atom_coords()` and `atom_charges()` read the built
integral data, including changes made directly to `_env` or `_atm`.

### Basis sets

`mol._basis` maps atom symbols or labels to lists of shells. Each shell has
one of these forms:

```text
[l, [exponent, coefficient_1, coefficient_2, ...], ...]
[l, kappa, [exponent, coefficient_1, coefficient_2, ...], ...]
```

Each row is a primitive Gaussian with a positive exponent; each coefficient
column is a contraction. The optional integer `kappa` selects the spinor block
and defaults to zero. Parsers split combined SP input into s and p shells.

`gto.basis.load(name, symbol)` and `gto.basis.parse(text)` return shell lists
for `mol.basis` dictionary values. `gto.format_basis()` formats the dictionary;
`build()` also accepts a common basis name or a `default` entry.
An exact atom-label entry, such as `'H@1'`, takes precedence over an element
entry such as `'H'`.

Take care when combining tagged atoms with `default`. An element entry does
not override the default basis for a tagged atom. For atoms labeled `H` and
`H@1`, these assignments give:

| `mol.basis` | Basis for `H` | Basis for `H@1` |
| --- | --- | --- |
| `{'H': 'cc-pvdz'}` | cc-pVDZ | cc-pVDZ |
| `{'default': 'sto-3g', 'H': 'cc-pvdz'}` | cc-pVDZ | STO-3G |
| `{'default': 'sto-3g', 'H': 'cc-pvdz', 'H@1': 'cc-pvtz'}` | cc-pVDZ | cc-pVTZ |

To override the default for a tagged atom, specify its exact label.
For example, when `default` is set, changing the `'H'` entry does not change
the basis for `'H@1'`. Set `'H@1'` separately.

Building absorbs primitive normalization into `_env` coefficients and, by
default, normalizes contractions. Use `bas_exp()` and `bas_ctr_coeff()` to
access basis data with these normalization conventions accounted for.

### Effective core potentials

The basis and ECP are separate inputs. Selecting a basis designed for an
ECP does not automatically attach that ECP. Specify both `mol.basis` and
`mol.ecp` when the intended calculation requires them. A molecule can build
successfully even when the intended ECP is missing.

`mol._ecp` maps symbols to the following structure:

```text
[number_of_core_electrons,
 [[l, [terms_for_radial_order_0, terms_for_radial_order_1, ...]], ...]]
```

`l=-1` denotes the local channel; `l>=0` denotes a projector channel.
Rows are `[exponent, scalar_coefficient]`, with an optional third SOC
coefficient. Empty lists retain missing radial orders. In the NWChem convention,
order `n` multiplies `r**(n-2) * exp(-exponent*r**2)` in the potential.

Use `gto.basis.load_ecp()` or `gto.basis.parse_ecp()` to obtain this data.
The core-electron count reduces the charge used by the built molecule and
the number of explicitly treated electrons. `atom_nelec_core(i)` reports the
removed core count.

SOC coefficients include the factor `2/(2*l+1)`. Scalar and SOC terms are
stored separately in `_ecpbas`; `has_ecp_soc()` detects SOC terms. Adding
them to an RHF/UHF Hamiltonian requires explicit customization; see the PySCF
examples `gto/20-soc_ecp.py` and `scf/44-soc_ecp.py`.

### GTH pseudopotentials

`mol._pseudo` stores the CP2K/GTH representation by atom symbol:

```text
[valence_electrons_by_channel,
 rloc, number_of_local_coefficients, local_coefficients,
 number_of_projector_channels,
 [r_l, number_of_projectors, h_matrix], ...]
```

`sum(valence_electrons_by_channel)` is the effective ionic charge. `rloc` is
the local Gaussian radius, followed by the local polynomial coefficients.
Nonlocal blocks start at `l=0` in increasing angular momentum. Each `h_matrix`
is a symmetric `(number_of_projectors, number_of_projectors)` matrix expanded
from CP2K's triangular input. Zero-projector channels may be present.

`gto.format_pseudo()` resolves named potentials supplied in a dictionary.
`pyscf.gto.basis.parse_cp2k_pp.parse()` parses CP2K potential text into one
such value. These parameters are distinct from the Gaussian ECP format.
During building, conflicting ECP and GTH assignments to the same symbol raise
`RuntimeError`. ECPs and pseudopotentials are not assigned to ghost atoms by
default.

## Serialization

Molecular data have three serialization forms:

- JSON, which stores the built molecule and most instance attributes.
- Python pickle data, which embeds the JSON representation.
- A dictionary containing selected inputs for rebuilding a molecule.

PySCF checkpoint files use the JSON representation for the `mol` record.

### JSON

`gto.dumps(mol)` and `mol.dumps()` return JSON. `gto.loads(text)` and
`Mole.loads(text)` create a molecule from it. Loading restores the internal
arrays; another `build()` is not required. `mol.loads_(text)` updates an
existing molecule.

The attributes `output`, `stdout`, `_keys`, and `_ctx_lock` are excluded.
Other non-JSON-serializable attributes may be dropped with a warning. Values
nested in lists or tuples can still make serialization fail.

The symmetry attributes `symm_orb`, `irrep_id`, and `irrep_name` are also
excluded and reconstructed during deserialization.

`atom`, `basis`, `ecp`, and `pseudo` are stored as Python representation
strings and reconstructed with `eval()`. Deserialize only trusted records.


### Pickle

Pickle serialization uses `Mole.__getstate__` and `Mole.__setstate__`, which
delegate to the JSON representation. It has the same exclusions and is not a
lossless way to preserve arbitrary custom attributes.

### Input-only dictionaries

`gto.pack(mol)` returns a Python dictionary of selected input attributes. It
excludes internal arrays and derived data. Use another serializer to store or
transfer the dictionary.

`gto.unpack(input_dict)` creates an unbuilt `Mole`. Call `build()` before using
its integrals:

```python
mol1 = gto.unpack(gto.pack(mol)).build()
```

The complete `pack()` result is not a `gto.M()` or `build()` keyword mapping:
it includes keys such as `_nelectron` and `nucprop` that those functions do not
accept.
