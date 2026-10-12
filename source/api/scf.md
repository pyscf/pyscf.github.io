# SCF customization API

This page describes the interfaces used by the molecular SCF driver. It is
intended for implementations that replace the Hamiltonian or part of the SCF
procedure. It does not list all SCF settings.

## SCF driver contract

`SCF.kernel()` and `SCF.scf()` are aliases. They run the same driver, update
the SCF object, and return the final total energy. The lower-level
`pyscf.scf.hf.kernel()` function instead returns
`(converged, e_tot, mo_energy, mo_coeff, mo_occ)` and does not provide the
normal object-level setup and finalization by itself.

For each calculation, the object-level driver performs these operations:

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
convergence; callers that require a converged reference should inspect
`converged`.

`max_cycle <= 0` requests a non-iterative energy evaluation. If the object has
not been initialized, the driver still constructs an initial wavefunction. If
orbitals already exist, this mode preserves them and updates only `e_tot`.


## SCF results

After a calculation, the SCF object holds the total energy in `e_tot` and
orbitals in `mo_energy`, `mo_coeff`, and `mo_occ`. `cycles` records the number
of iterations. `converged` indicates whether the calculation converged. Results
are stored even when convergence is not reached.

Typically, (for instance, the RHF method), the shape of `mo_coeff` is a square
matrix. Each column of `mo_coeff` represents an MO, and each row corresponds to
an AO. `mo_energy` and `mo_occ` vectors give the energy and occupation of each
column. When AO functions are linearly dependent, `mo_coeff` has fewer columns.
The `mo_energy` and `mo_occ` arrays shrink accordingly.

UHF and UKS store alpha and beta orbitals separately. Their `mo_coeff` has shape
`(2, nao, nmo)`, and `mo_energy` and `mo_occ` have shape `(2, nmo)`.

GHF and GKS store a single set of orbitals. General GHF/GKS orbitals need not be
eigenvectors of the spin operator $S_z$. The orbital can mix alpha and beta
spin. For generalized orbitals, `scf.ghf.guess_orbspin(mo_coeff)` returns 0 for a
pure alpha orbital, 1 for a pure beta orbital, and -1 when both spin components
are present.

Converting an RHF/RKS or UHF/UKS object with `to_ghf()` or `to_gks()` initially
preserves pure alpha and beta orbitals. The conversion records their spin
labels in `mo_coeff.orbspin`. For closed-shell RHF/RKS, each spatial orbital
becomes an adjacent alpha--beta pair with the same energy. For UHF/UKS,
occupied orbitals come first, followed by unoccupied orbitals. Each group is
sorted by energy rather than arranged in alpha--beta pairs.

The `reset(mol)` method of SCF classes replaces the molecular system and clears
integral caches, but keeps these result attributes. Calling `kernel()` again
without an explicit initial density uses the stored orbitals and occupations as
the initial guess. Note, they are reused without projection, even if the
geometry, basis, or atom ordering has changed. (If a basis change or atom
reordering, `init_guess_by_chkfile()` can construct one by projecting checkpoint
orbitals onto the current basis).

In addition to the result attributes, the `scf_summary` attribute of SCF class
records the energy components created by SCF iterations. This attribute is a
dictionary. Its common keys include `e1` (one-electron energy), `e2`
(two-electron energy), `coul` (Coulomb energy), `exc` (exchange energy for HF or
exchange-correlation energy for DFT), `gap` (HOMO-LUMO gap), and `nuc` (nuclear
repulsion). Depending on the calculation, it can also contain `dispersion`,
`e_solvent`. Note, energy components are in Hartree, while `gap` is in eV.

The `reset()` method clears this dictionary.


## Defining a custom Hamiltonian

For a model Hamiltonian in an AO-like basis, configure the following pieces:

- Set the electron count on the `Mole` object. For an unrestricted calculation,
  `mf.nelec = (nalpha, nbeta)` can be used to specify the two spin populations.
- Override `get_hcore(mol)` and `get_ovlp(mol)`. They normally return square
  `(nao, nao)` matrices. [For UHF with different alpha-beta core hamiltonian,
  and add a hyperlink to the subsection which provides detailed description]
- Supply the two-electron interaction through `_eri` or override `get_jk()`.
  [For AO integrals in spinor basis or complex-valued orbitals]
- Supply an initial density through `kernel(dm0=...)` when the molecule does
  not provide a meaningful atomic initial guess.
- Override `energy_nuc()` when the constant part of the model Hamiltonian is
  not `mol.enuc`.

For real, permutation-symmetric integrals, `_eri` should use the four-fold or
eight-fold packed format accepted by `pyscf.ao2mo`. `ao2mo.restore(8, eri,
nao)` converts a full four-index tensor to the common eight-fold format. The
number of basis functions implied by `_eri`, `get_hcore()`, `get_ovlp()`, and
the initial density must agree.

Assigning methods on one object is supported and convenient for a single
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

`check_sanity()` reports instance methods such as `get_hcore` and `get_ovlp`
as overwritten attributes. This diagnostic is expected for this customization
style. A reusable implementation should normally subclass the appropriate
RHF, ROHF, or UHF class instead.

### add a description for different alpha and beta core Hamiltonians.
    - the UHF get_hcore produces alpha-beta degenerated core Hamiltonian (N,N) matrix.
    - A spin-dependent UHF core Hamiltonian may have shape
      `(2, nao, nao)`; a common two-dimensional Hamiltonian is broadcast to both
      spins.
    - get_fock and energy_elec can use spin-resolved hcore
      while the overlap matrix is shared between the two eigensolves.


### two-electron integrals for GHF
The custom integrals supplied to the ._eri attribute must be integrals for
scalar basis.
  - integrals with 4-fold and 8-fold permutation symmetry, and N^4 without
    permutation symmetry, which of them are supported
  - Is complex integrals supported.
  - integrals of two-component spinor bases are not supported. To handle such
    case, `get_jk` needs to be customized accordingly. However, the subsequent
    MO integral transformation in post-HF calculations does not automatically
    recogonize this extension, it will still break the post-HF.

## Use by post-HF methods

After overwriting the Hcore and `_eri`, post-HF methods (MP, CI, CC, MCSCF, etc)
can use these integrals from a customized SCF calculation.

Some transformation paths choose an out-of-core calculation from `mol` when
the stored `_eri` is considered too large for the memory limit. For a model
Hamiltonian, that fallback would generate the molecule's ordinary integrals
instead. Set

```python
mol.incore_anyway = True
```

to force supported post-HF paths to use `mf._eri`. This can require substantial memory.

Method composition, including density fitting, X2C, solvent models, and SOSCF,
may replace customized methods. It's recommend apply these method composition
first, then apply customization.

Analytical gradients and Hessians require derivatives of the customized
Hamiltonian, overlap, and two-electron interaction with respect to nuclear
coordinates. Overriding their ground-state values does not supply those
derivatives. The ordinary molecular derivative routines can therefore be
inconsistent with a custom Hamiltonian unless their derivative interfaces are
customized as well.

Electronic response has a different requirement: it needs the change in the
effective potential caused by a density perturbation. `gen_response()` returns
a callable that maps a perturbed density to that induced potential. Standard
HF response routines use J/K builders, so customizing `get_jk()` can also serve
response calculations. Overriding `get_veff()` alone does not update
this response kernel.


## Method type

Each SCF method is associated to a class, such as `RKS`, `UKS`, `GKS`, `ROHF`,
and `DHF`, etc.
Derived methods such as density fitting, solvent models, and QM/MM are build on
these classes through [dynamic composition](#dynamically-composed-methods).
For example, `mol.RKS().density_fit().PCM().soscf()` creates an object whose
class is named `SecondOrderPCMDFRKS`.

To check the underlying method, one method is to use `isinstance()` with the
relevant class, e.g.:
```python
isinstance(mf, rks.RKS)
isinstance(mf, gks.GKS)
```
Another approach is to use the `istype()` method of SCF classes, e.g.:
```python
mf.istype('RKS')
mf.istype('GKS')
```


## J/K builders

`get_jk()` returns `(vj, vk)`, the Coulomb and exchange matrices for the supplied
density. It accepts a single `(nao, nao)` matrix or a batch with shape
`(..., nao, nao)`. The output J and K matrices has the same shape as the input
density matrices.

For UHF/UKS, a ground-state density has shape `(2, nao, nao)`, with alpha first
and beta second. A batch of unrestricted response densities can have shape
`(2, npert, nao, nao)`, where `npert` counts perturbations. `get_jk()` contracts
each density separately. Periodic k-point methods have additional k-point
conventions. The spin-channels still placed at the leading dimension. The next
dimension is k-point indices.

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

Following hooks are provided in SCF driver. They can access all variables in the
SCF driver via a Python `locals()` dict.

- `check_convergence(envs)` replaces the default energy-and-orbital-gradient
  test and must return a boolean. Returning `True` requests convergence; the
  optional final convergence check can call it again.

- `pre_kernel(envs)` runs before the first iteration, after the initial density, potential, energy, and
  optional DIIS object have been constructed.

- `callback(envs)` runs at the end of each regular iteration, after the
  convergence test and checkpoint write. Its return value is ignored. One common
  scenario is to use it to observe or save iteration data.

- `post_kernel(envs)` runs after the iteration loop and any final convergence
  check, before the driver returns.

## Checkpoint files

Use `mf.chkfile` to specify an HDF5 file as checkpoint. By default, an SCF
object creates a temporary file that is removed with the temporary-file object.
Assign a path to retain the checkpoint, or set `mf.chkfile = None` to disable
checkpoint writes.

The checkpoint is an HDF5 file. Its `scf/` group stores the energy and orbital
arrays, rather than a serialization of the SCF object:

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

The `mol` dataset stores a serialized `Mole` object.
`mol` can be deserialized from a chkfile using `lib.chkfile.load_mol(filename)`.

## Dynamically composed methods
Density fitting, X2C, solvent models, and SOSCF are provided via the method class composition.

The order matters.
Although some composition are exchangable, as they modify different part of the code, but no guarantee.
It's recommended to perform the following order
X2C -> density_fitting -> PCM -> SOSCF

- X2C (mf.x2c()), density fitting (mf.density_fit()), and solvent models (mf.PCM()) are decoupled at mean-field level. They can be applied in arbitrary order.

- Density fitting and SOSCF are coupled. Their order matters. Applying `density_fit().soscf()` will use SOSCF for the density fitting mean-field calculation. It should produce the same results as `density_fit().run()` when converged. Applying in a flipped order `mf.soscf().density_fit()` will use density fitting approximation to speed up integral computation for orbital hessian in SOSCF, this would produce the same results as the basic `mf.run()` when converged

- in the derivative computations, the current implementation only supports `density_fit().PCM()` type of method composition. Although the order of `density_fit()` and `PCM()` does not affect SCF results, it's recommended to perform `density_fit().PCM()` initialization than the flipped way due to the limitation of derivatives functionality.

- QM/MM module also provide an initializer for mean-field objects, which can apply MM charges as external field to the Hcore.
  QM/MM cannot be used with PCM. Their composition are not properly implemented, will lead to undefined results.
  QM/MM and X2c order matters. when applying both QM/MM and X2C, X2C must be applied first. performing X2c on top of QM/MM are not handled in current implementation










## CCSD recomputes a Fock matrix for the supplied orbitals

**Status:** current convention. **Sources:** [Discussion #2650](https://github.com/pyscf/pyscf/discussions/2650).

**Question or historical answer.** Supplying mo_energy from another program does not necessarily control the CC denominators.

**Current code.** _ChemistsERIs._common_init_ builds a density, calls the reference get_veff/get_fock, transforms the Fock matrix and takes its real diagonal as eris.mo_energy. This step is not a full SCF optimization. A custom CC ao2mo implementation can control its ERI object, but its integrals, Fock matrix and reference must agree.

**Documentation task.** Document the CC ERI-object boundary and the distinction between reference mo_energy, eris.mo_energy and the full Fock matrix. Treat the subclass snippet in the thread as an advanced customization, not a standalone guarantee of physical consistency.

**Code evidence:** [`pyscf/cc/ccsd.py:1393`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/cc/ccsd.py#L1393).

## Decorators and custom optimizer hooks

**Status:** current convention with caveats. **Sources:** [Issue #2212](https://github.com/pyscf/pyscf/issues/2212), [Issue #2788](https://github.com/pyscf/pyscf/issues/2788), [Issue #874](https://github.com/pyscf/pyscf/issues/874), [Discussion #2802](https://github.com/pyscf/pyscf/discussions/2802), [Discussion #3242](https://github.com/pyscf/pyscf/discussions/3242), [Discussion #3034](https://github.com/pyscf/pyscf/discussions/3034).

**Question or historical answer.** SOSCF, density fitting and method conversion can change which customized hooks are used.

**Current code.** DF copies instance state but clears _eri and supplies its own JK and derivative methods. SOSCF uses an underlying _scf reference, so assigning a method only to an outer wrapper is not a universal customization recipe. Repeated CASSCF density_fit now checks for an existing DF mixin; a local repeat call succeeded. The old conversion/DF-loss issue has corresponding conversion changes, but custom subclasses still need individual checks.

**Documentation task.** Extend scf.md's composition section with hook ownership and a supported-combination matrix. Cross-link the ML-XC discussion rather than implying a single override survives every wrapper. Omit the old repeated-DF MRO workaround.

**Code evidence:** [`pyscf/df/df_jk.py:109`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/df/df_jk.py#L109), [`pyscf/soscf/newton_ah.py:670`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/soscf/newton_ah.py#L670), [`pyscf/mcscf/df.py:30`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/mcscf/df.py#L30), [`pyscf/dft/uks.py:207`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/dft/uks.py#L207), [`pyscf/scf/hf.py:2430`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/scf/hf.py#L2430).

## Complex Hamiltonian support is solver and path specific

**Status:** current limitation. **Sources:** [Discussion #2628](https://github.com/pyscf/pyscf/discussions/2628), [Discussion #2659](https://github.com/pyscf/pyscf/discussions/2659), [Discussion #3184](https://github.com/pyscf/pyscf/discussions/3184).

**Question or historical answer.** A successful complex SCF or RCCSD run is taken to imply unrestricted CCSD support.

**Current code.** The molecular UCCSD implementation still uses optimized real/symmetric integral contractions; the shared _contract_s4vvvv_t2 asserts float64 amplitudes. GCCSD has explicit complex-capable paths, but its integral representation and orbital conventions differ. Removing a dtype assertion is not enough to establish correct complex conjugations throughout UCCSD.

**Documentation task.** Add a solver/input/path support matrix, with complex MO coefficients distinguished from arbitrary complex custom AO ERIs. Keep the reported local patch and suggested GCCSD comparison as validation leads, not as upstream support claims.

**Code evidence:** [`pyscf/cc/uccsd.py:454`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/cc/uccsd.py#L454), [`pyscf/cc/ccsd.py:487`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/cc/ccsd.py#L487), [`pyscf/cc/gccsd.py:344`](https://github.com/pyscf/pyscf/blob/9ceafa41b217bd3afbe6a0399d875d805ff81306/pyscf/cc/gccsd.py#L344).
