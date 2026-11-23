# General API patterns

This page describes API patterns shared by PySCF method objects.

## Stream objects

Most method objects inherit from `pyscf.lib.StreamObject`. Their `kernel()`
method performs the calculation and returns results that depend on the method.
These results are normally also stored as attributes on the method object.
See each method class's docstring for the relevant attributes.

Stream objects provide common operations that let you chain method calls:

- `run(*args, **kwargs)` first applies keyword arguments with `set`, calls
  `kernel(*args)`, and returns the same object rather than the kernel result.

  Keyword arguments to `run()` are settings, not arguments to `kernel()`.
  In particular, `mf.run(dm0=dm)` sets an attribute named `dm0`; it does not
  supply the initial density to SCF. Use `mf.kernel(dm0=dm)` or `mf.run(dm)`
  to start SCF from that density. Settings such as `conv_tol` can be passed
  through `run(conv_tol=...)`.

  **Intermediate `run()` calls.** Some constructors and property methods can be
  chained without calling `run()` between them. Any required earlier
  calculation is run automatically. For example, the following TDDFT gradient
  calculations are equivalent:

  ```python
  td_grad = pyscf.M(atom='H 0 0 0; H 0 0 1.1').RHF().run().TDA().run().Gradients().run()
  ```

  ```python
  td_grad = pyscf.M(atom='H 0 0 0; H 0 0 1.1').RHF().TDA().Gradients().run()
  ```

  You still need to call `run()` or `kernel()` at the end to complete the calculation.

- `set(**kwargs)` assigns attributes and returns the same object. It supports
  keyword arguments only and does not check attribute names.

- `copy()` makes a shallow copy of the object.

- `view(cls)` returns a view with the same attributes and a different class.


### Output control

`StreamObject` provides two default attributes for controlling output. See the
logger API page for logging levels and timer utilities.

* `stdout = sys.stdout`: by default, PySCF writes output to standard output,
  usually the terminal. To redirect the default output to a file, set
  `StreamObject.stdout` to that open file.

* `verbose = 0`: `logger.Logger` uses this default verbosity level. Classes
  inheriting from `StreamObject` can be passed to `Logger`. The default is
  silent. Most methods take their initial verbosity from `mol.verbose`, so
  changing the class default may not affect objects that already exist.

### Attribute spelling check with `check_sanity()`

`check_sanity()` checks public instance attributes (names not starting with `_`).

- Supported attributes are registered in each class's `_keys` set. `_keys` from
  parent classes are checked too, so a child class does not need to register
  the same keys again.

- A message such as `<class ...> does not have attributes ...` often indicates
  a typo or an attribute omitted from `_keys`, rather than a calculation error.

- It reports unknown public attributes and methods that have been replaced
  on an object. For example, an `RHF` object does not support the setting
  `spin`. Assigning it produces the following message:

```python
mf = mol.RHF()
mf.spin = 1
mf.check_sanity()
```

```text
<class 'pyscf.scf.hf.RHF'> does not have attributes  spin
```

Replacing the `get_hcore` method on an `RHF` object:

```python
mf = mol.RHF()
mf.get_hcore = lambda *args, **kwargs: None
mf.check_sanity()
```

can produce:

```
Overwritten attributes  get_hcore  of <class 'pyscf.scf.hf.RHF'>
```

### Reset method

Use `method.reset(new_mol)` to reuse a method with a new molecule or cell.
It clears saved integrals and other temporary data for the old system while
keeping settings such as convergence thresholds and frozen orbitals.

Results from the previous calculation are kept and may be reused as an
initial guess. `reset()` does not run a calculation, so stored energies and
convergence flags still refer to the previous calculation.

When scanning geometries, changing the `Mole` or `Cell` object alone may
not update every part of the calculation. Several method objects can share
the same molecule, while an attached solvent model may keep a separate copy
of the geometry. Assigning `mf.mol = new_mol` also leaves old integrals and
attached objects unchanged.

The recommended approach is to call `method.reset(new_mol)`. This also updates
the attached objects supported by the method, such as `_numint`, `with_df`,
and `with_solvent`, so the calculation does not mix data from different geometries.

SCF keeps its orbitals and occupations after `reset()`. Calling `kernel()`
without an explicit initial density reuses them as they are; it does not
adapt them to a new basis or atom order. If either has changed, supply an
initial density that matches the new system. For example,
`init_guess_by_chkfile()` can project saved orbitals into the new basis,
as described under [checkpoint restoration](scf.md#restoring-a-calculation).

[Scanners](#scanners) handle the geometry update and run the calculation in
one call. They also manage the reuse of previous results according to the
method.


## Scanners

`as_scanner()` creates a scanner that can be called repeatedly with different
geometries. Calling an energy scanner with a `Mole` or `Cell` runs the
calculation and returns the total energy. A nuclear-gradient scanner returns
`(energy, gradient)`. Results are also available on the scanner:

* Energy scanners store `e_tot` and `converged`.
* Gradient scanners store `e_tot`, `de` (nuclear gradients), and `converged`.

```python
scanner = gto.M().RHF().as_scanner()
energy = scanner(gto.M(atom='H 0 0 0; H 0 0 1.2'))
print(scanner.converged)
```

The scanner keeps the method settings and, when possible, uses the previous
result as the initial guess for the next geometry. This is useful for
potential-energy scans and geometry optimization. Each call updates the
scanner's results to those of the latest calculation.

**Side effects.** A scanner is not an independent copy of the original
method. The two can share attached objects, such as density-fitting, X2C,
or solvent objects. A scan can update these shared objects to the new
geometry while the original method's `mol` or `cell` still refers to the
old geometry. For example:

```python
mf = gto.M(atom='H 0 0 0; H 0 0 1.0').RHF().density_fit()
scanner = mf.as_scanner()
scanner(gto.M(atom='H 0 0 0; H 0 0 1.2'))
```

After this scan, `mf.mol` and `mf.with_df.mol` refer to different geometries.
Running `mf` again without resetting it can mix data from the two geometries.

Before reusing the original `mf`, reset it to the geometry you want:

```python
# Return mf and its DF object to the original geometry.
mf.reset(mf.mol)

# Or update both to the scanner's latest geometry.
mf.reset(scanner.mol)
```

### Orbital reuse in CASCI and CASSCF scanners

CASCI and CASSCF energy scanners reuse orbitals differently:

- `casci_scanner(mol)`: Uses new SCF orbitals, even if the CASCI object already
  has saved orbitals.
- `casci_scanner(mol, mo_coeff=mo)`: Uses the supplied orbitals without an SCF
  calculation or automatic projection.
- `casscf_scanner(mol)`: Adapts saved CASSCF orbitals to the new geometry by
  projection. If no orbitals have been initialized, it uses new SCF orbitals.
- `casscf_scanner(mol, mo_coeff=mo)`: Adapts the supplied orbitals by projection,
  then uses them as the CASSCF starting guess.

Assigning orbitals to `casci_scanner.mo_coeff` does not make the next scan
use them. This differs from SCF calculations, which automatically reuse
existing orbitals as the initial guess when no initial density is supplied.

A CASSCF scanner runs an SCF calculation even when starting orbitals are
supplied. Both scanners reuse their previous CI solution as a guess unless
an explicit initial guess `ci0` is supplied. After the call, the new orbitals
and CI solution are available as `scanner.mo_coeff` and `scanner.ci`.


### Serialization
TODO


### Conversion to GPU4PySCF

`to_gpu()` converts a supported PySCF object to its corresponding GPU4PySCF
object. It also converts supported attributes, including NumPy arrays and
attached method objects.

Some arrays intentionally remain NumPy arrays, including `kpt`, `kpts`,
`_kpts`, `kpts_band`, `mesh`, and `frozen`. These exceptions are listed in
`_ATTRIBUTES_IN_NPARRAY` in the `pyscf.lib.misc` module.

`from_cpu()` is the complementary conversion hook and is implemented by
GPU4PySCF modules or classes. When GPU4PySCF provides `from_cpu()`, it takes
precedence over PySCF's `to_gpu` conversion code.

PySCF and GPU4PySCF evolve separately, so the supported methods and conversion
results can vary between GPU4PySCF versions.

Not every PySCF method or scanner supports `to_gpu()`. The shared conversion
code still relies on support from individual methods.


## Dynamically composed method classes

Features such as density fitting (DF), one-electron X2C, second-order SCF
(SOSCF), and solvent models add features to an existing method by creating a
new class at runtime. The new class keeps the original method's behavior and
adds the requested feature. This is called decorating the method.

These classes are made available through `pyscf.lib` when they are created.
They are not defined in source files, but you can access them through
`pyscf.lib` like any other class.

To remove a decoration, use the corresponding `undo_*()` method. For example:

```python
mf = gto.M(atom='H 0 0 0; F 0 0 1.1').RHF().density_fit()
plain_mf = mf.undo_df()
```

Some features can also be disabled temporarily without removing the decoration.
For example, setting `mf.with_df = None` bypasses DF-based J/K construction while
leaving the composed class unchanged.

To permanently remove a supported decoration, use the corresponding undo method,
such as `undo_df()`, `undo_x2c()`, or `undo_soscf()`. These methods return a
view of the object with the corresponding mixin removed.
