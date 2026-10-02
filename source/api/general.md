# General API patterns

This page describes API patterns shared by PySCF method objects.

## Stream objects

Most method objects inherit from `pyscf.lib.StreamObject`. Their `kernel()`
method performs the calculation and has a method-specific output. These output
are normally stored as result attributes on the method object. Consult each
method class's docstring for the relevant attributes.

Stream objects provide common operations that allow methods to be chained in a
computation:

- `run(*args, **kwargs)` first applies keyword arguments with `set`, invokes
  `kernel(*args)`, and returns the same object rather than the kernel result.

  **Intermediate `run()` calls.** Some constructors and property methods can be
  chained without an intermediate `run()`. A prerequisite calculation is
  triggered automatically. For example, the following TDDFT gradient workflows
  are equivalent

  ```python
  td_grad = pyscf.M(atom='H 0 0 0; H 0 0 1.1').RHF().run().TDA().run().Gradients().run()
  ```

  ```python
  td_grad = pyscf.M(atom='H 0 0 0; H 0 0 1.1').RHF().TDA().Gradients().run()
  ```

  Calling `run()` or `kernel()` at the end to complete a calculation is still required.

- `set(**kwargs)` assigns attributes and returns the same object. It supports
  keyword arguments only and does not validate attribute names.

- `copy()` makes a shallow copy of the object.

- `view(cls)` returns a view with the same attributes and a different class.


### Output control

`StreamObject` provides two default attributes for controlling logging. See the
logger API page for logging levels and timer utilities.

* `stdout = sys.stdout`: by default, PySCF output is written to the standard
  system stdout. If you want to redirect default output to a file, you can
  set `StreamObject.stdout` to that open file.

* `verbose = 0`: `logger.Logger` uses this default verbosity level. Classes
  inheriting from `StreamObject` can be passed to `Logger`. The default is
  silent. Most methods initialize their verbosity from `mol.verbose`. So
  changing a class-level default may not affect existing method instances.

### Attribute spelling check with `check_sanity()`

`check_sanity()` checks public instance attributes (names not starting with `_`).

- Supported attributes are registered in each class's `_keys` set. `_keys` from
  parent classes are also included in the checking. Key registration does not
  need to be repeated in a child class.

- A message such as `<class ...> does not have attributes ...` often indicates
  a typo or an attribute omitted from `_keys`, rather than a calculation error.

- It reports an unknown public attribute or an overwritten of a class method.
  For example, `RHF` instance does not support configuration `spin`, setting
  this attribute leads to

```python
mf = mol.RHF()
mf.spin = 1
mf.check_sanity()

<class 'pyscf.scf.hf.RHF'> does not have attributes  spin
```

  Overwriting `the get_hcore` method of the `RHF` instance,

```python
mf = mol.RHF()
mf.get_hcore = lambda *args, **kwargs: None
mf.check_sanity()
```

can report:

```
Overwritten attributes  get_hcore  of <class 'pyscf.scf.hf.RHF'>
```

### Reset method

`reset()` is implemented by individual method classes. It updates the method's
associated system and clears cached data tied to the old system while preserving
method settings such as convergence thresholds and frozen orbitals.

It does not clear results from a previous calculation. This can allow old
results to be used as an initial guess for a new system.

It is common to run the same calculation settings for molecules with slightly
different geometries (or other configurations). A `Mole` (or `Cell`) object may
be referenced by several method objects. Updating it in place can appear to work
in some cases, but is not generally safe. A nested object such as an
implicit-solvent model may keep its own copy of the molecule geometry
separately.

`method.reset(new_mol)` can update the system recursively. Cached intermediates
are cleared during this update to avoid mixing quantities pre-computed from
different systems.

`Scanner <Scanners>_` is another tool designed for this workflow. It
automatically calls the available `reset` method and performs additional
method-specific checks between systems.


### Serialization
TODO


### Conversion to GPU4PySCF

`to_gpu()` converts a supported PySCF object to its corresponding GPU4PySCF
object. It recursively converts supported attributes, including NumPy arrays
and nested method attributes.

Some arrays intentionally remain NumPy arrays, including `kpt`, `kpts`,
`_kpts`, `kpts_band`, `mesh`, and `frozen`. This setting is currently provided
by a global list, `_ATTRIBUTES_IN_NPARRAY` in the `pyscf.lib.misc` module.

`from_cpu()` is the complementary conversion hook and is implemented by
GPU4PySCF modules or classes. When GPU4PySCF provides `from_cpu()`, it takes
precedence over PySCF's `to_gpu` conversion code.

PySCF and GPU4PySCF evolve separately, support and return details can vary
between GPU4PySCF versions.

Not every PySCF method or scanner supports `to_gpu()`. Although a global
conversion is implmented, the to gpu conversion still relies on each individual
method for many methods. 


## Dynamically composed method classes

Features such as density fitting (DF), one-electron X2C, second-order SCF
(SOSCF), and solvent models decorate an existing method by constructing a new
class at runtime. The resulting class retains the behavior of the original class
while adding the requested functionality.

These dynamically composed classes are created and registered in the
`pyscf.lib` module at runtime. They do not exist in the source tree, but can be
accessed like any other class through `pyscf.lib`.

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


## Scanners

`as_scanner()` creates a callable, stateful derivative of a method object for
evaluating related geometries. Calling an energy scanner with a `Mole`, `Cell`
returns total energy. A nuclear-gradient scanner returns `(energy, gradient)`.
Additional attributes are accessible via the scanner:
* Energy scanners expose `e_tot` and `converged`.
* Gradient scanners expose `e_tot`, `de` (nuclear gradients), and `converged`.

```python
scanner = gto.M().RHF().as_scanner()
energy = scanner(gto.M(atom='H 0 0 0; H 0 0 1.2'))
print(scanner.converged)
```

The scanner retains method options and reuses the previous result as an initial
guess when compatible, which is useful for potential-energy scans and geometry
optimization. After a call, its normal result attributes are updated, as the
original method class.

**Side effects.** Scanners share nested attributes with the method from which
they were created, such as `_scf`, `with_df`, `with_x2c`, or solvent models.
A scan can leave these nested attributes changed for the last geometry while
the parent method still refers to its original `mol` or `cell`. For example:

```python
mf = gto.M(atom='H 0 0 0; H 0 0 1.0').RHF().density_fit()
scanner = mf.as_scanner()
scanner(gto.M(atom='H 0 0 0; H 0 0 1.2'))
```

After this scan, `mf.mol` and `mf.with_df.mol` refer to different systems.
Running `mf` directly can then mix data associated with the original molecule
and the scanned molecule.

If you intend to reuse the original `mf`, restore consistent system state first:

```python
# Restore the parent and its nested DF object to the parent's original geometry.
mf.reset(mf.mol)

# Or update the parent and its nested DF object to the scanner's latest geometry.
mf.reset(scanner.mol)
```
