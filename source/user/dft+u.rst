.. _dftplusu:

DFT+U
*****

The DFT+U functional is an effective method to address self-interaction errors
in strongly localized electrons such as the d or f electrons in transition metals.
In these systems, the local or semi-local DFT functionals often over delocalizes
these electrons, leading to wrong electronic structures.

In PySCF, Hubbard U can be combined with all XC functionals, including
LDA, GGA, meta-GGA, even the hybrid DFT functionals to for the DFT+U functional.

DFT+U has been implemented separately for molecular and periodic
boundary condition (PBC) calculations. For periodic setups, users can invoke
the `Cell.KRKSpU` and `Cell.KUKSpU` methods to create the DFT+U instance. These
two classes are provided by the `pbc.dft.krkspu` and `pbc.dft.kukspu` modules.
For molecular systems, the `Mole.RKSpU` and `Mole.UKSpU` methods are available.
There is no specialized implementation for periodic systems at the gamma point.
The DFT+U calculations at the gamma point should be handled using the k-point
sampling version `KRKSpU` and `KUKSpU`.

The DFT+U functionals can be integrated with various other modules, such as the
X2C method for relativistic effect corrections, the PCM methods to include
solvation effects.
In addition to calculating the ground state energy, analytical gradients have
been implemented for DFT+U methods. The optimization of molecular geometries and
crystal structures using DFT+U functionals is availabel.
However, please note that the orbital-Hessian for DFT+U functionals
is not available. DFT+U cannot be used for methods that require
orbital Hessians, such as analytical nuclear Hessian, and TDDFT calculations.


How to perform DFT+U
====================

For molecular systems, we use the RKSpU and UKSpU classes to create a DFT+U
mean-field object. To create a DFT+U instance, you need to provide the `U_idx`
and `U_val` parameters to specify which orbitals the U term will be applied to.
For example::

    import pyscf
    mol = pyscf.M(
        atom='Mn 0. 0. 0.; O 1.6 0 0',
        basis='ccpvdz-dk',
        spin=5,
        verbose=4
    )
    # Add U term (in eV) to 3d orbitals
    U_idx = ["Mn 3d"]
    U_val = [2.8]
    mf = mol.UKSpU(xc='pbe', U_idx=U_idx, U_val=U_val)
    mf = mf.x2c()
    mf.run()

In this example, `U_idx` specifies that the U term is applied to the 3d orbitals
of the Mn atom, and the `U_val` parameter sets the value of the effective `J-U`
term to 2.8 eV.

The DFT+U code provides various options to configure `U_idx` and `U_val`, which
allows us to use individual settings for atoms or orbitals within the system.
More details of the `U_idx` and `U_val` configurations will be explained in the
following sections.

`U_idx` and `U_val` are saved as attributes of the DFT+U mean-field instance.
The `U_val` in the input arguments, as well as the `U_val` attribute of the
mean-field instance, are provided in eV. During the SCF iteration, the U
parameters are converted to atomic units (AU) to adapt to the existing
DFT Fock matrix implementation. When customizing the DFT+U code, the differences
in units should be noted.

Various features for mean-field instances can be combined with this DFT+U
mean-field object, such as relativistic corrections, the implicit solvent model,
and density fitting approximations. Relativistic corrections are included in
this example. The second-order SCF method can also be applied to the DFT+U
instance. However, the effect of U is not taken into account in the orbital
Hessian computation in the current PySCF implementation. As a result, one may
experience slower convergence for the DFT+U calculations than standard Kohn-Sham DFT.

The DFT+U calculations for crystal systems are similar to those for molecular
systems. Calling the RKSpU and UKSpU classes will create a DFT+U instance at the
gamma point. For DFT+U calculations with k-point sampling, we use the KRKSpU and
KUKSpU classes to instantiate a DFT+U mean-field object. For example::

    import pyscf
    cell = pyscf.M(
        atom='C 0.,  0.,  0.; C 0.8917,  0.8917,  0.8917',
        a='''0.      1.7834  1.7834
             1.7834  0.      1.7834
             1.7834  1.7834  0.    ''',
        basis='gth-dzvp',
        pseudo='gth-pade',
        verbose=4)
    kmesh = [2, 2, 2]
    kpts = cell.make_kpts(kmesh)
    # Add U term (in eV) to the 2p orbital of the second Carbon atom
    U_idx = ['1 C 2p']
    U_val = [5.0]
    mf = cell.KRKSpU(kpts=kpts, U_idx=U_idx, U_val=U_val, minao_ref='gth-szv')
    mf.run()

The DFT+U method in PySCF by default employs MINAO as the reference basis to
localize the orbitals. MINAO is an all-electron basis set that includes only core
and valence orbitals. 

Pseudopotentials are commonly used in periodic boundary condition (PBC)
calculations. Using all-electron orbitals as a reference for orbital
localization may not yield optimal valence orbitals. In this example, we use the
GTH-SZV basis as the localization reference, which is compatible with the
pseudopotential. This setting is provided to the `minao_ref` argument.

Configurable Options
--------------------
The DFT+U class provides four attributes that can be configured for DFT+U calculations.

* `U_idx`: The orbital indices to which the U modifications are applied.
* `U_val`: The effective U-J term [in eV].
* `C_ao_lo`: The orbital coefficients of custom local orbitals.
* `alpha`: the perturbation [in eV] used to compute U using the linear-response
  DFT+U calculations. Reasonable values are 0.01 - 0.05. A smaller `alpha` value
  could lead to large errors due to its impact to the SCF convergence accuracy.

Input Formats for U_idx and U_val
---------------------------------
The DFT+U module supports multiple input methods for `U_idx` and `U_val`, which
enables us to apply apply different U values to various atoms and orbitals
within the system.

`U_idx` is a list where each element specifies the orbitals to which U will be
applied. `U_val` is also a list. Its length must match that of the `U_idx` list.
Each element of the `U_val` list records a U value for the corresponding element
specified in `U_idx`.

The elements of `U_idx` support two input formats. They can be a list of
integers or a string. When an element of `U_idx` is a list of integers, it
specifies the indices of the AOs for applying U. If a string is specified for
the `U_idx` element, it represents a collection of AOs that match the string
representation. For example, `'Mn 3d'` matches all the 3d type orbitals for all
Mn atoms; `'1 C 2pz'` refers to the `2pz` orbital of the second atom, which
points to a single orbital. The string input is parsed by the
`mol.search_ao_label()` method to locate the indices of the corresponding
orbitals. You can use this function to validate the input string.

`U_idx` also supports the mix of both input formats. For example, you can input
an `U_idx` as `['Mn 3d', [6, 7, 8, 15, 16, 17]]`.

It is important to note that the elements of `U_idx` should be core and valence orbitals.
You should avoid including high-lying AOs. This is because the current DFT+U
implementation utilizes a set of reference orbitals to localize the orbitals to
which U is applied. The high-lying AOs might not be available in the reference
orbitals. Sepcifying high-lying AOs in `U_idx` can lead to poor orbital
localization or may even result in a complete failure of the calculation.

Orbital Localization
--------------------
In the current DFT+U implementations, U is applied to a set of localized
orbitals. The DFT+U code uses a set of reference orbitals to identify the core
and valence orbitals, and localizes them through Lowdin orthogonalization. By
default, the reference orbitals are based on the MINAO basis set. Other
reference orbitals can be set through the input keyword `minao_ref` or the
`minao_ref` attribute of the DFT+U instance.

The format for `minao_ref` is the same as the format for the `mol.basis` input
in PySCF. For example, a truncated ANO for O 3s 2p can be input as `{'O': 'ANO@3S2P'}`.
If you have generated a set of custom atomic orbital references, you can convert
the AO contraction coefficients in the same way as the PySCF basis set input and
assign them to `minao_ref`.

The default configuration for `minao_ref`, which is the MINAO basis set, is
suitable for most all-electron calculations. However, in some scenarios,
`minao_ref` can be configured for additional computational requirements.

* The MINAO basis set does not include data of fourth-period or heavier elements.
  If a system contains these heavy elements, we can use other AO
  references, such as the ANO basis set as the `minao_ref`.

* The MINAO basis set only includes parameters for core and valence orbitals. If
  you need to apply a U correction to high-lying AOs, you will need to configure a
  reference basis set that includes the high-lying AOs.

* The calculations with ECP or pseudopotential is another scenario that requires
  modifying `minao_ref`. In these cases, the valence orbitals differ from the
  valence orbitals in all-electron calculations, particularly near the
  nucleus. Using the MINAO basis set as a reference may lead to poor valence
  orbitals. In such calculations, we can assign an appropriate reference basis
  set to `minao_ref`. For example, for the GTH pseudopotential, we can use the
  GTH-svp basis as the reference basis set.

In addition to using a reference basis set to localize orbitals, the DFT+U
program also supports customized local orbitals. The local orbitals can be
obtained using the `pyscf.lo` module, and the orbital coefficients can be
assigned to the keyword argument `C_ao_lo` or the `C_ao_lo` attribute of the
DFT+U instance. 

It is important to note that when the custom local orbitals is assigned, the
program will automatically skip the orbital localization and orthogonalization
process. The function for Fock matrix in DFT+U will assume that these local
orbitals are orthogonal. If the custom local orbitals are not orthogonal, the U
contribution may be double counted.

Another limitation for `C_ao_lo` customization is that the current program does
not support the calculations of nuclear gradients for custom orbitals.
This feature cannot be used for tasks such as geometry optimization or crystal
structure optimization.

Geometry Optimization and Crystal Structure Optimization
========================================================
The nuclear gradients of the DFT+U method are supported by the PySCF package. We
can use the `geomeTRIC` package with the DFT+U gradients to optimize the
molecular geometry::

    mol = pyscf.M(
        atom='Mn 0. 0. 0.; O 1.6 0 0',
        basis='ccpvdz',
        spin=5)
    U_idx = ["Mn 3d"]
    U_val = [2.8]
    mf = mol.UKSpU(xc='pbe', U_idx=U_idx, U_val=U_val)
    mol_eq = mf.Gradients().optimizer(solver='geomeTRIC').kernel()

For crystal calculations, geometry optimization can be performed using the
PySCF-ASE interface::

    from pyscf.pbc.tools.pyscf_ase import PySCF, pyscf_to_ase_atoms
    from ase.optimize import BFGS
    cell = pyscf.M(
        atom='C 0.,  0.,  0.; C 0.8917,  0.8917,  0.8917',
        a='''0.      1.7834  1.7834
             1.7834  0.      1.7834
             1.7834  1.7834  0.    ''',
        basis='gth-dzvp',
        pseudo='gth-pade')
    U_idx = ['1 C 2p']
    U_val = [5.0]
    kmesh = [2, 2, 2]
    kpts = cell.make_kpts(kmesh)
    mf = cell.KRKSpU(xc='pbe', kpts=kpts, U_idx=U_idx, U_val=U_val, minao_ref='gth-szv')

    # Setup the ASE calculator from PySCF objects
    atoms = pyscf_to_ase_atoms(cell)
    atoms.set_calculator(PySCF(method=mf))
    opt = BFGS(atoms).run()
    print(atoms.get_positions())

The DFT+U module also provides the capability to compute the stress tensor. This
functionality is also integrated into the PySCF-ASE interface. We can utilize the
crystal structure optimization features provided by ASE to optimize the crystal
structure::

    from ase.filters import StrainFilter
    opt = BFGS(StrainFilter(atoms)).run()
    print(atoms.cell)

or simultaneously optimize the atom positions and the crystal structure::

    from ase.filters import UnitCellFilter
    opt = BFGS(UnitCellFilter(atoms)).run()
    print(atoms.get_positions())
    print(atoms.cell)


Determining U with Linear-response Calculations
===============================================

One challenge in DFT+U calculations is the determination of the appropriate value of the
Hubbard U parameter for a specific system. The linear-response U method offers a
systematic approach to generate the Hubbard U for individual systems.

For the theoretical background of the linear-response U method,
the following references provide comprehensive insights:
* M. Cococcioni and S. de Gironcoli, Phys. Rev. B 71, 035105 (2005)
* H. J. Kulik, M. Cococcioni, D. A. Scherlis, and N. Marzari, Phys. Rev. Lett. 97, 103001 (2006)
* A tutorial by H. J. Kulik: [Calculating Hubbard U](https://hjkgrp.mit.edu/tutorials/2011-05-31-calculating-hubbard-u/)

A linear-response calculation of the effective Hubbard parameter U involves
multiple self-consistent DFT+U perturbation calculations. The procedure
consists of the following steps:

1. A standard DFT+U calculation without any external perturbation. This step
obtains the reference electronic structure and occupation matrices :math:`n_I` for
each correlated site :math:`I`.

2. Apply a series of small on-site potential shifts :math:`a_I` as a perturber

.. math::
    \alpha_I \sum{m} \hat{n}_{Im}

where :math:`n_{Im}` is the occupation operator for orbital m on site I.

3. For each perturbation,

  * Record the occupations before the first SCF iteration, :math:`n_I^0`,
    representing the bare (non-interacting) response.

  * Record the occupations after full self-consistency, :math:`n_I`,
    representing the interacting response.

4. Linear regression for :math:`\chi_I^0` :math:`\chi_I` against the
   :math:`alpha_I` perturbation.

.. math::
   \begin{align}
   \chi_I &= \frac{\partial n_I}{\partial \alpha_I} \\
   \chi_I^0 &= \frac{\partial n_I^0}{\partial \alpha_I}
   \end{align}

5. The effective Hubbard interaction parameter

.. math::
    U_I = \frac{1}{\chi_I^0} - \frac{1}{\chi_I}

PySCF provides the `linear_response_u()` API to perform the calculation of the
Hubbard U parameter. This function returns an optimized U parameter (in eV unit).
It can be directly utilized in subsequent DFT+U calculations. The following
example demonstrates the use of this API to determine the U parameter for the
sextet state of the MnO molecule::

    from pyscf.dft.ukspu import linear_response_u
    mol = pyscf.M(atom='Mn 0 0 0; O 0 0 1.6',
                  basis='ccpvdz-dk',
                  spin=5, # sextet, 5 un-paired electrons
                  verbose=4)
    # Hubbard U on Mn 3d shells only.
    # Scalar relativisitic corrections are considered in this system.
    mf = mol.UKSpU(xc='pbe', U_idx=['Mn 3d'], U_val=[3.0]).x2c().run()
    u0 = linear_response_u(mf)
    print(f'linear response U = {u0} eV')

In this example, the initial U value is set to 3.0 eV. The X2C scalar relativistic
corrections are utilized here. However, their impact is minimal in this case.
The linear regression analysis and the corresponding U parameter can be found at
the end of the output.::

    Line fitting chi0 = -0.228159 x + 5.340062
    Line fitting chif = -0.129000 x + 5.340046
    Uresp = 3.369056, chi0 = -0.228159, chif = -0.129000
    linear response U = 3.3690559047534574 eV

The linear-response U calculations can be repeatedly invoked to
self-consistently optimize the U parameter. The self-consistent U typically
converges quickly. Within two to three iterations, typically, the change of U
value would drop to under 0.1 eV.


Examples
========

* :source:`examples/dft/22-dft+u_linear_response_u.py`
* :source:`examples/pbc/22-dft+u.py`
