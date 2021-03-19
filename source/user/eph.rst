.. _user_eph:

**********************
Electron Phonon Matrix
**********************

*Modules*: :mod:`eph`, :mod:`pbc.eph`

Introduction
============
Electron-Phonon (eph) Interaction describes the coupling between the electrons and the nuclear degree of freedom. First order eph matrix elements are most commonly used to compute electron/phonon transport properties. PySCF supports first order eph matrix calculation by either Coupled Perturbed Kohn Sham/Hartree Fock (CPKS/CPHF) or finite difference. In periodic system, only Gamma-point calculation with FFTDF is supported.

A minimal example of analytical eph matrix calculation is as follows::

    from pyscf import gto, scf, eph
    mol = gto.M(
    atom = [['O', [0.000000000000,  -0.000000000775,   0.923671924285]],
            ['H', [-0.00000000000,  -1.432564848017,   2.125164039823]],
            ['H', [0.000000000000,   1.432564848792,   2.125164035930]]],
    basis = 'sto3g',
    unit = 'Bohr'
    )
    # Note that this is a pre-computed relaxed molecule
    mf = scf.RHF(mol).run()

    myeph = eph.EPH(mf)
    ephmat, omega = myeph.kernel()
    print('eph matrix', ephmat)
    print('phonon frequencies', omega)

Alternatively, one can compute the matrix elements using finite difference::

    from pyscf import gto, scf, eph
    mol = gto.M(
    atom = [['O', [0.000000000000,  -0.000000000775,   0.923671924285]],
            ['H', [-0.000000000000,  -1.432564848017,   2.125164039823]],
            ['H', [0.000000000000,   1.432564848792,   2.125164035930]]],
    basis = 'sto3g',
    unit = 'Bohr'
    )
    mf = scf.RHF(mol).run()

    ephmat, omega = eph.eph_df.kernel(mf, disp=1e-4)
    print('eph matrix', ephmat)
    print('phonon frequencies', omega)

Spin symmetry
=============
The eph module in PySCF supports mean field methods including Hartree Fock and DFT using either spin-restricted or spin-unrestricted reference.

Analytical eph matrix elements
------------------------------

The module-level ``eph.EPH(mf)`` constructor can infer the correct method based
on the level of symmetry-breaking in the mean-field argument.  For more explicit
control or inspection, the respective classes and functions can be found in
``rhf.py``, ``uhf.py``, ``rks.py`` and ``uks.py``.

For example, a spin-unrestricted calculation can be performed
as follows::

    mf = scf.UHF(mol).run()
    myeph = eph.uhf.EPH(mf)
    ephmat, omega = myeph.kernel()
    print('eph matrix', ephmat)
    print('phonon frequencies', omega)

Finite difference
-----------------

The finite difference kernel function ``eph.eph_fd.kernel`` can handle all the supported mean field methods above.

Output control
==============

Filtering noisy frequencies
---------------------------

The diagonalization of the mass-weighted hessian can sometimes yield unphysical vibrational modes that correspond to imaginary frequencies or to non-zero frequencies which should strictly vanish. These unphysical modes may arise either from the use of not fully relaxed geometries, the use of insufficiently accurate density functional quadrature, or simply numerical noise. By default, modes below 80 cm-1 are filtered out in PySCF; note that this also includes any negative frequencies.

One can specify a different cutoff frequency (in cm-1) when constructing the EPH object as follows::

    myeph = eph.EPH(mf, cutoff_frequency=60)
    ephmat, omega = myeph.kernel()

Similar to keeping imaginary frequencies, one can set keep_imag_frequency to True when initializing the object::

    myeph = eph.EPH(mf, keep_imag_frequency=True)
    ephmat, omega = myeph.kernel()

Matrix element representation
-----------------------------

The eph matrix are computed in the atomic orbital (AO) basis by default. One can also request them in the MO basis as follows::

    myeph = eph.EPH(mf)
    ephmat, omega = myeph.kernel(mo_rep=True)
