# Periodic modules conventions

This page describes electron counts, smearing settings, and orbital arrays
for periodic calculations. For calculation setup, see the
[periodic SCF user guide](../user/pbc/scf.rst).

## Cell charge and spin convention

For a calculation with `nkpts` equally weighted k-points, `cell.charge` and
`cell.spin` specify the charge and spin of the **sampled supercell**, containing
`nkpts` copies of the input unit cell. As with `Mole.spin`, `cell.spin` is
the electron-number difference `Nalpha - Nbeta`, not the multiplicity `2S+1`.

For symmetry-reduced k-point sampling, `nkpts` counts the k-points in the
full mesh, including those omitted by symmetry.

- Setting `cell.charge = 1` removes one electron from the sampled supercell.
  The charge per unit cell is then `1 / nkpts`. To keep the same charge per
  unit cell when changing the mesh, scale `cell.charge` with `nkpts`.
- Setting `cell.spin = 2` requests two more alpha than beta electrons in
  the sampled supercell. As in molecular calculations, `spin` means
  `Nalpha - Nbeta`, not the multiplicity `2S+1`.

The `cell.nelectron` attribute uses a different convention: it counts the
electrons in a single cell. For a charged system (`cell.charge != 0`), this
count treats the charge as belonging to one cell. Do not use this attribute
as the total electron count for k-point sampled calculations.
Assigning `cell.nelectron` modifies the electron count **per unit cell**.

The `cell.nelec` attribute is a pair `(nalpha, nbeta)`. Like `cell.nelectron`,
it refers to a single cell, giving the alpha and beta electron counts.

For k-point sampled calculations, `KUHF.nelec` and `KUKS.nelec` give the
alpha and beta counts for the sampled supercell. For a neutral system with
`cell.spin=0`, each population is `nkpts` times the corresponding population
in `cell.nelec`.

Unless `mf.nelec` is explicitly overridden, the following relations can be
expected:

```python
nalpha, nbeta = mf.nelec
nalpha + nbeta == cell.tot_electrons(nkpts)
nalpha - nbeta == cell.spin
```

TODO: Update the electron-number rules after merging PR 3447 (fractional
electron numbers for `Mole` and `Cell`).

## Smearing width, electron targets and free energies

The following conventions apply to
`mf.smearing(sigma, method, mu0=None, fix_spin=False)`.

### Width and units

`sigma` is an energy width in **Hartree** for both smearing methods
(`method='fermi'` and `method='gaussian'`).

For `method='fermi'`, convert a temperature in Kelvin or a width in eV to
`sigma` using the following conversion factors:

```python
from pyscf.data import nist

sigma = temperature_in_K * nist.BOLTZMANN / nist.HARTREE2J  # ~= temperature_in_K * 3.167e-6
sigma = width_eV / nist.HARTREE2EV  # ~= width_eV / 27.211
```

### Electron targets

With `mu0=None`, the electron target applies across the full k-point mesh,
not separately at each k-point.
Default smearing for KUHF/KUKS can relax the spin-population constraint.
For example, it can give unequal total alpha and beta populations even
with `cell.spin=0`.

Use the `fix_spin` parameter to control the spin-population constraint.

- `fix_spin=False` (default): Neither `cell.spin` nor an assigned `mf.nelec` fixes the separate spin populations.
- `fix_spin=True`: Alpha and beta populations match `mf.nelec` separately.

Supplying `mu0` fixes the chemical potential rather than providing an
initial guess. This takes precedence over the electron targets: the
resulting electron count can differ from the count requested through
`cell.charge`, `cell.nelectron`, or `mf.nelec`, even with `fix_spin=True`.

### Energy attributes

Enabling smearing does not change `mf.e_tot` or the return value of
`kernel()` into a free energy. Both remain the HF/DFT total energy at the
fractional-occupation density. After the calculation, the energy choices
are:

| Attribute | Energy |
| --- | --- |
| `mf.e_tot` | $E$, without an entropy correction |
| `mf.e_free` | $E - \sigma S$, including the entropy correction |
| `mf.e_zero` | $E - \tfrac{1}{2}\sigma S$, an estimate of the zero-smearing energy |

Here $S$ is `mf.entropy`. All three energies are in Hartree **per unit
cell**, even though the electron target refers to the sampled supercell.
The entropy is also per unit cell.

## PBC orbital coefficients

K-point orbital arrays place the k-point axis before the AO and orbital
axes. Unrestricted arrays place spin first. Normally, `mo_coeff` has the
following shape:

- KRHF/KRKS: `(nkpts, nao, nao)`
- KUHF/KUKS: `(2, nkpts, nao, nao)`

Removing linear dependencies preserves these shapes, even if different
k-points retain different numbers of orbitals.
When removing basis functions reduces the number of MOs at a k-point,
the coefficient array is padded with zero columns.
The corresponding `mo_energy` entries are `1e30` Hartree, and the
`mo_occ` entries are zero.

These placeholder columns must be handled appropriately or excluded during
post-SCF processing.

Single-k-point RHF/RKS and UHF/UKS omit the k-point axis and have coefficient
shapes `(nao, nmo)` and `(2, nao, nmo)`.

Code consuming k-point orbitals must support complex coefficients, even
when a Gamma-only mesh produces real orbitals.

For single-k-point calculations:

| Calculation | Expected coefficients |
| --- | --- |
| RHF/UHF or RKS/UKS at Gamma with a real Hamiltonian | Real |
| Any method at a general nonzero k-point | Complex |
| GHF/GKS at Gamma | May be complex |
