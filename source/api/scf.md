# SCF customization API

This page describes how to customize molecular SCF calculations by replacing
the Hamiltonian or parts of the SCF procedure. It does not list all SCF settings.

## SCF driver contract

`SCF.kernel()` and `SCF.scf()` are aliases. They run the same driver, update
the SCF object, and return the final total energy. The lower-level
`pyscf.scf.hf.kernel()` function instead returns
`(converged, e_tot, mo_energy, mo_coeff, mo_occ)` and does not provide the
usual setup and final steps handled by the SCF object.

For each calculation, the SCF object performs these steps:

1. Report settings with `dump_flags()` and initialize caches with `build()`.
2. Select an initial density matrix with `get_init_guess()` or reuse the
   existing orbitals when possible.
3. Obtain the AO overlap and one-electron Hamiltonian from `get_ovlp()` and
   `get_hcore()`.
4. Build the effective potential with `get_veff()` and evaluate the total
   energy with `energy_tot()`.
5. In each iteration, call `get_fock()`, `eig()`, `get_occ()`, `make_rdm1()`,
   `get_veff()`, and `energy_tot()`.
6. Test convergence using `get_grad()` and either the default thresholds or
   `check_convergence`.

The driver stores its final state in `converged`, `e_tot`, `mo_energy`,
`mo_coeff`, `mo_occ`, and `cycles`. A finite energy does not imply
convergence. Check `converged` before using the result as a converged reference.

`max_cycle <= 0` evaluates the energy without SCF iterations. If the object has
not been initialized, the driver still constructs an initial wavefunction. If
orbitals already exist, this mode preserves them and updates only `e_tot`.


## SCF results

After a calculation, the SCF object holds the total energy in `e_tot` and
orbitals in `mo_energy`, `mo_coeff`, and `mo_occ`. `cycles` records the number
of iterations. `converged` indicates whether the calculation converged. Results
are stored even when convergence is not reached.

For methods such as RHF, `mo_coeff` is usually a square matrix. Each column
represents an MO, and each row corresponds to an AO. The `mo_energy` and
`mo_occ` vectors give the energy and occupation of each column. When linearly
dependent AO functions are removed, `mo_coeff` has fewer columns.
The `mo_energy` and `mo_occ` arrays shrink accordingly.

UHF and UKS store alpha and beta orbitals separately. Their `mo_coeff` has shape
`(2, nao, nmo)`, and `mo_energy` and `mo_occ` have shape `(2, nmo)`.

GHF and GKS store a single set of orbitals. General GHF/GKS orbitals need not be
eigenvectors of the spin operator $S_z$. An orbital can mix alpha and beta
spin. For generalized orbitals, `scf.ghf.guess_orbspin(mo_coeff)` returns 0 for a
pure alpha orbital, 1 for a pure beta orbital, and -1 when both spin components
are present.

Converting an RHF/RKS or UHF/UKS object with `to_ghf()` or `to_gks()` initially
preserves pure alpha and beta orbitals. The conversion records their spin
labels in `mo_coeff.orbspin`. For closed-shell RHF/RKS, each spatial orbital
becomes an adjacent alpha--beta pair with the same energy. For UHF/UKS,
occupied orbitals come first, followed by unoccupied orbitals. Each group is
sorted by energy rather than arranged in alpha--beta pairs.

The `reset(mol)` method updates the molecule and clears stored integrals,
but keeps these result attributes. Calling `kernel()` again
without an explicit initial density uses the stored orbitals and occupations as
the initial guess. They are reused without projection, even if the geometry,
basis, or atom order has changed. If the basis or atom order changes,
`init_guess_by_chkfile()` can construct an initial density by projecting
saved orbitals onto the current basis.

The `scf_summary` dictionary records energy components from the SCF
iterations. Its common keys include `e1` (one-electron energy), `e2`
(two-electron energy), `coul` (Coulomb energy), `exc` (exchange energy for HF or
exchange-correlation energy for DFT), `gap` (HOMO-LUMO gap), and `nuc` (nuclear
repulsion). Depending on the calculation, it can also contain `dispersion`
and `e_solvent`. Energy components are in Hartree, while `gap` is in eV.

The `reset()` method clears this dictionary.


## Defining a custom Hamiltonian

For a model Hamiltonian in an AO-like basis, set up the following:

- Set the electron count on the `Mole` object. For an unrestricted calculation,
  `mf.nelec = (nalpha, nbeta)` can be used to specify the two spin populations.
- Override `get_hcore(mol)` and `get_ovlp(mol)`. They normally return square
  `(nao, nao)` matrices. UHF also accepts
  [different alpha and beta core Hamiltonians](#different-alpha-and-beta-core-hamiltonians).
- Supply the two-electron interaction through `_eri` or override `get_jk()`.
  Real integrals with four-fold or eight-fold permutation symmetry can use
  the corresponding packed format. Complex integrals must use a full
  four-index tensor.
  See the [limitations](#limitations) for GHF, complex integrals, and method
  composition.
- Alternatively, supply the two-electron interaction through
  [three-index density-fitting integrals](#density-fitting-integrals).
- Supply an initial density through `kernel(dm0=...)` when the molecule does
  not provide a meaningful atomic initial guess.
- Override `energy_nuc()` when the constant part of the model Hamiltonian is
  not `mol.enuc`.

`ao2mo.restore(8, eri, nao)` converts a full four-index tensor of real,
permutation-symmetric integrals to the eight-fold packed format. The
number of basis functions implied by `_eri`, `get_hcore()`, `get_ovlp()`, and
the initial density must agree.

Replacing methods on an individual object is convenient for a single
calculation:

```python
import numpy
from pyscf import ao2mo, gto, scf

mol = gto.M()
mol.nelectron = 2

h1 = numpy.array([[0.0, -1.0], [-1.0, 0.0]])
eri = numpy.zeros((2, 2, 2, 2))
eri[0, 0, 0, 0] = eri[1, 1, 1, 1] = 2.0

mf = scf.RHF(mol)
mf.get_hcore = lambda mol=None: h1
mf.get_ovlp = lambda mol=None: numpy.eye(2)
mf._eri = ao2mo.restore(8, eri, 2)
mf.kernel(dm0=numpy.eye(2))
```

`check_sanity()` reports replaced methods such as `get_hcore` and `get_ovlp`
as overwritten attributes. This message is expected when customizing an
object this way. For reusable code, normally create a subclass of the
appropriate RHF, ROHF, or UHF class instead.

### Different alpha and beta core Hamiltonians

The standard UHF solver uses the same `(nao, nao)` core Hamiltonian for both
spins. To supply different alpha and beta Hamiltonians, return an array of
shape `(2, nao, nao)` from `get_hcore()`, with alpha first and beta second.
The UHF Fock matrices and energy then include the corresponding spin-dependent
one-electron terms. Both spins still use the same `(nao, nao)` overlap
matrix; separate alpha and beta overlap matrices are not supported by the
standard UHF solver.

### Density fitting integrals

Ordinary `density_fit()` does not fit the custom tensor in `mf._eri`.
Instead, a DF-SCF object uses three-index integral tensors stored in
`mf.with_df._cderi`. To supply custom integrals, first create a DF-SCF object
with `density_fit()`, then assign the three-index tensor to that attribute.

For the standard molecular DF solver, supply a real `float64` array with
shape `(naux, nao * (nao + 1) // 2)` in C-contiguous order. Each row contains
the lower triangle of a symmetric AO matrix, ordered as
`(0,0), (1,0), (1,1), (2,0), ...`.

These tensors are factors of the two-electron interaction, not raw
three-center integrals. They are related to `_eri` in the four-fold packed
format by:

```python
_eri = einsum('pi,pj->ij', _cderi, _cderi)
```

`_cderi` can be a NumPy array held in memory or a path to an HDF5 file.
The HDF5 file must contain the three-index tensor in a dataset named `j3c`.

Assigning `mf._eri` to a DF-SCF object does not change the integrals used
by its DF J/K builder.

### Limitations

#### Two-electron integrals for GHF

GHF uses a spin-orbital basis for its one-electron quantities, but its
`_eri` attribute accepts only integrals in the **spatial AO basis**.
A full `_eri` tensor with `nao**4` elements is accepted. For real integrals, the
four-fold or eight-fold packed formats can also be used.

Assigning two-component spinor integrals or an arbitrary spin-dependent
two-electron interaction to GHF `_eri` is not supported. Such an interaction
requires a custom `get_jk()` that accepts and returns matrices in the
spin-orbital basis.

#### Customized hooks

SOSCF can use a custom Hamiltonian supplied through `get_hcore()`,
`get_ovlp()`, and `_eri` (or `_cderi` for DF objects).

Overriding iteration hooks such as `get_fock()`, `get_occ()`, or `get_grad()`
is a different case. SOSCF may not call these hooks in the same
way as the ordinary SCF driver. A calculation may run and produce results
without following the behavior you intended.

#### Complex Hamiltonians and post-HF methods

RHF, UHF, and GHF can use complex Hermitian one-electron Hamiltonians and
complex orbitals. A complex spatial AO two-electron interaction must be
supplied as a full tensor.

Support at SCF level does not guarantee support in a post-HF solver or
its AO-to-MO transformation. Some post-HF methods, such as
`pyscf.cc.rccsd.RCCSD`, support complex integrals.
Other CC methods, such as UCCSD and EOM-CC, do not generally support
complex integrals or orbitals.

## Use by post-HF methods

After replacing `get_hcore()` and `_eri`, post-HF methods such as MP, CI,
CC, and MCSCF can use these integrals from the customized SCF calculation.

Some integral transformations use a disk-based calculation with integrals
from `mol` when `_eri` is too large for the memory limit. For a model
Hamiltonian, this would replace the custom integrals with the molecule's
ordinary integrals. Set

```python
mol.incore_anyway = True
```

to force supported post-HF methods to use `mf._eri`. This can require a
large amount of memory.

Method composition, including density fitting, X2C, solvent models, and SOSCF,
may replace customized methods. Apply these features before customizing
the methods they replace. For SOSCF, follow the reference setup described
under [customized hooks](#customized-hooks).

Analytical gradients and Hessians require derivatives of the customized
Hamiltonian, overlap, and two-electron interaction with respect to nuclear
coordinates. Supplying their values for the ground-state calculation does
not supply these derivatives. Standard molecular gradient and Hessian
calculations may therefore give incorrect results for a custom Hamiltonian
unless the derivative methods are customized too.

Electronic response calculations need the change in the effective potential
caused by a small change in density. `gen_response()` returns a function
that converts a density change into the corresponding potential change. Standard
HF response routines use J/K builders, so customizing `get_jk()` can also serve
response calculations. Overriding `get_veff()` alone does not update
this response kernel.


## Method type

Each SCF method has a class, such as `RKS`, `UKS`, `GKS`, `ROHF`, or `DHF`.
Features such as density fitting, solvent models, and QM/MM are built on
these classes through [dynamic composition](#dynamically-composed-methods).
For example, `mol.RKS().density_fit().PCM().soscf()` creates an object whose
class is named `SecondOrderPCMDFRKS`.

To check the underlying method, use `isinstance()` with the relevant class:

```python
isinstance(mf, rks.RKS)
isinstance(mf, gks.GKS)
```

Alternatively, use the SCF object's `istype()` method:

```python
mf.istype('RKS')
mf.istype('GKS')
```


## J/K builders

`get_jk()` returns `(vj, vk)`, the Coulomb and exchange matrices for the supplied
density. It accepts a single `(nao, nao)` matrix or a batch with shape
`(..., nao, nao)`. The output J and K matrices have the same shape as the input
density matrices.

For UHF/UKS, a ground-state density has shape `(2, nao, nao)`, with alpha first
and beta second. A batch of unrestricted response densities can have shape
`(2, npert, nao, nao)`, where `npert` counts perturbations. `get_jk()` contracts
each density separately. Periodic k-point methods have additional k-point
conventions: the spin axis still comes first, followed by the k-point axis.

When a custom `get_jk()` method is provided, the inherited `get_j()` and
`get_k()` also use the custom builder.


### Incremental potential builds

The driver passes the current density and the previous density/potential pair
to `get_veff(mol, dm, dm_last, vhf_last)`. The method must return the **full
potential for `dm`**, not just its change since the previous iteration.

For HF, the potential is linear in the density. An incremental build therefore
evaluates the potential of `dm - dm_last` and adds it to `vhf_last`. The previous
potential must correspond to `dm_last` for the same Hamiltonian and AO basis.
With `direct_scf` enabled and no stored `_eri`, the standard molecular HF
implementation uses this approach. With stored integrals, or when direct SCF
is disabled, it rebuilds the potential from the full density.

A custom `get_veff()` may ignore `dm_last` and `vhf_last` and always rebuild
from `dm`. It must still accept these arguments because the driver supplies
them. Density differences are not ordinary physical densities: a custom J/K
builder must handle them without assuming positive occupations or a fixed
electron count. For nonlinear potentials such as exchange-correlation, the HF
increment formula cannot be applied to the entire potential; use a full build
or a method-specific update.


## Hooks for SCF iteration

The SCF driver provides the following hooks. They receive an `envs`
dictionary containing the driver's local variables, obtained with Python's
`locals()`.

- `check_convergence(envs)` replaces the default energy-and-orbital-gradient
  test and must return a boolean. Returning `True` requests convergence; the
  optional final convergence check can call it again.

- `pre_kernel(envs)` runs before the first iteration, after the initial
  density, potential, energy, and optional DIIS object have been prepared.

- `callback(envs)` runs at the end of each regular iteration, after the
  convergence test and checkpoint write. Its return value is ignored. A common
  use is to monitor or save iteration data.

- `post_kernel(envs)` runs after the iteration loop and any final convergence
  check, before the driver returns.

## Checkpoint files

Use `mf.chkfile` to choose an HDF5 checkpoint file. By default, an SCF
object creates a temporary file that is removed when its temporary-file
object is released. Assign a path to keep the checkpoint, or set `mf.chkfile = None` to disable
checkpoint writes.

The checkpoint is an HDF5 file. Its `scf/` group stores the energy and orbital
arrays, rather than the entire SCF object:

- `scf/e_tot`: the current total energy;
- `scf/mo_energy`: orbital energies;
- `scf/mo_coeff`: orbital coefficients;
- `scf/mo_occ`: orbital occupations.

The dataset names match the SCF result attributes, so saved results can be
assigned directly:

```python
from pyscf.scf import chkfile

mf.__dict__.update(chkfile.load('chkfile_name', 'scf'))
```

The `mol` dataset stores the `Mole` object in serialized form. Use
`lib.chkfile.load_mol(filename)` to restore it.

### Restoring a calculation

Loading the `scf` results from a PySCF checkpoint file does not restore all SCF
settings. Despite its name, a PySCF checkpoint does not contain everything
needed to repeat a calculation. The standard SCF checkpoint stores the molecule
and wavefunction results. Settings such as the DFT functional (`xc`), auxiliary
basis, relativistic corrections, solvent model, and convergence thresholds
are not saved. These settings and method combinations
must be supplied again in a new calculation.

SCF can use existing `mo_coeff` and `mo_occ` values to construct an initial
guess. It assumes that the orbitals use the same AO order and basis as the
current system. Supplying `mo_coeff` does not automatically adapt the
orbitals to the current basis or atom order.

SCF also does not check that the supplied orbitals are orthonormal. When
reusing saved orbitals, for example as an initial guess for SOSCF, make sure
they are orthonormal with respect to the current AO overlap matrix.

## Dynamically composed methods

Density fitting, X2C, solvent models, and SOSCF can be added to an existing
method object.

The order matters. Some features can be applied in either order because
they affect different parts of the calculation, but this is not always true.
The recommended order is:

X2C → density fitting → PCM → SOSCF

- X2C (`mf.x2c()`), density fitting (`mf.density_fit()`), and solvent models
  (`mf.PCM()`) are independent at the mean-field level and can be applied
  in any order for the SCF calculation itself.

- Density fitting and SOSCF depend on their order. Applying
  `density_fit().soscf()` uses SOSCF to optimize the density-fitted mean-field
  calculation. When converged, it should give the same results as
  `density_fit().run()`. Reversing the order with `mf.soscf().density_fit()`
  uses density fitting only to speed up the orbital Hessian calculation in
  SOSCF. When converged, this should give the same results as the original
  `mf.run()` calculation.

- Derivative calculations currently support `density_fit().PCM()`, but not
  the reverse order. Although the order does not affect SCF results, use
  `density_fit().PCM()` if derivatives will be needed.

- The QM/MM module adds the field from MM charges to the core Hamiltonian
  of a mean-field calculation. Combining QM/MM with PCM is not supported
  and can give undefined results. When using both QM/MM and X2C, apply X2C
  first; applying X2C after QM/MM is not supported.
