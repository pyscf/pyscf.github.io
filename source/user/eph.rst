.. _user_eph:

************************
Electron-phonon coupling
************************

*Modules*: :mod:`eph`, :mod:`pbc.eph`

Introduction
============
Because nuclei are much heavier than the electrons, one usually assumes the Born-Oppenheimer (BO) approximation, which
defines an adiabatic electronic wavefunction for fixed nuclear positions. However, this approximation breaks down when the electronic states both change very rapidly with nuclear position and approach each other in energy.

To move beyond the BO approximation, one can assume that the electrons experience the moving nuclei as a perturbation of the external potential. Carrying out an expansion to linear order in nuclear positions and assuming quadratic nuclear motion around the minimum leads to an approximate Hamiltonian for the electron and nuclear motion 
known as the linear electron-phonon coupling approximation. When deriving this Hamiltonian, the electronic degrees of freedom are assumed to be treated within a mean-field theory. This approximation is widely used in solid state calculations to model phenomena arising from coupled electron/nuclear motion (see :cite:`giustino2017eph` for more discussion).

PySCF supports the calculation of the first-order electron-phonon (e-ph) coupling matrix elements using either coupled perturbed Kohn-Sham/Hartree-Fock (CPKS/CPHF) or finite differences. In periodic systems, only FFTDF-based Gamma-point calculations are supported by the finite difference implementation. In both cases, a geometrically relaxed structure is required for the consistency of the theory. 

A minimal example of an analytical e-ph matrix calculation is as follows::

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

Alternatively, one can compute the matrix elements using finite differences::

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
The ``eph`` module in PySCF supports mean field methods including Hartree Fock and DFT using either spin-restricted or spin-unrestricted reference.

Analytical eph matrix elements
------------------------------

The module-level ``eph.EPH(mf)`` constructor can infer the correct method based
on the type of the mf argument.  For more explicit
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

Similarly, to keep the imaginary frequencies, one can set keep_imag_frequency to True when initializing the object::

    myeph = eph.EPH(mf, keep_imag_frequency=True)
    ephmat, omega = myeph.kernel()

Matrix element representation
-----------------------------

The eph matrix is computed in the atomic orbital (AO) basis by default. One can also request it in the MO basis as follows::

    myeph = eph.EPH(mf)
    ephmat, omega = myeph.kernel(mo_rep=True)

References
==========

.. bibliography:: ref_eph.bib
  :style: unsrt
