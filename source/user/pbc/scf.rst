.. _user_pbc_scf:

*******************
SCF and DFT methods
*******************

*Modules*: :py:mod:`pyscf.pbc.scf`, :py:mod:`pyscf.pbc.dft`

Introduction
============

PySCF supports periodic Hartree-Fock and density functional theory calculations
with Brillouin zone sampling.  The results of these calculations serve
as input to a variety of periodic post-HF and post-DFT calculations.
A minimal example of a periodic HF calculation on diamond with
a 2x2x2 sampling of the Brillouin zone is shown below.  Note that the
:attr:`kpts` keyword argument must be in physical units (inverse bohr),
which can be easily achieved using the ``Cell.make_kpts()`` method::

    from pyscf.pbc import gto, scf
    cell = gto.M(
        atom = '''C 0.0000 0.0000 0.0000
                  C 0.8917 0.8917 0.8917''',
        a = '''0.0000 1.7834 1.7834
               1.7834 0.0000 1.7834
               1.7834 1.7834 0.0000''',
        pseudo = 'gth-pade',
        basis = 'gth-szv'
        )

    kmf = scf.KRHF(cell, kpts=cell.make_kpts([2,2,2])).run()
    # converged SCF energy = -10.9308552994574

See :ref:`user_pbc_gto` and :ref:`user_pbc` for more details about
the :class:`Cell` object and Brillouin zone sampling.

Various integral schemes have been developed for periodic calculations in PySCF.
A regular periodic calculation can be started with some type of density fitting.
For large systems, multi-grid algorithm (see :ref:`user_pbc_multigrid`) or
range-separation (RS) integral algorithm (see :ref:`user_pbc_rsjk`) can be
considered to reduce the cost of Fock build.
See :ref:`user_pbc_df_comparison` for the comparison of different integral
schemes in periodic calculations.

Density fitting
===============

The default behavior is to use plane-wave density fitting (FFTDF).
The number of plane-waves used as an auxiliary basis is controlled by kinetic
energy cutoff, which is specified by the :attr:`Cell.ke_cutoff` parameter.  The
default value of this parameter is chosen to provide many digits of precision in
the ERIs and the subsequent energies.  If reduced precision is tolerable, this
parameter can be set manually, resulting in significant speedups and memory
savings.  For the diamond example above, lowering the kinetic energy cutoff to 100
hartree changes the SCF energy by about 5 microhartree::

    cell.ke_cutoff = 100.0 # in Hartree
    cell.build()
    kmf = scf.KRHF(...).run()
    # converged SCF energy = -10.9308602914696

For many systems, Gaussian density fitting (GDF) is more economical, although it
incurs larger errors than FFTDF.  Periodic GDF can be activated in the same way
as for molecules::

    kmf = scf.KRHF().density_fit().run()

If a corresponding auxiliary basis is found for the chosen atomic orbital basis,
it will be used.  Otherwise, an even-tempered Gaussian basis will be used.

See :ref:`user_pbc_df` for more details about periodic density fitting.


Finite-size effects
===================

The long-ranged nature of the Coulomb interaction is responsible for a number of
divergent contributions to the energy.  For charge-neutral unit cells, the
divergence of the nuclear repulsion energy, the electron-nuclear attraction
energy, and the electronic Hartree energy cancel one another.

The nonlocal exact exchange energy exhibits an integrable divergence at
:math:`G=0` that is responsible for the leading-order finite-size error of HF
and hybrid DFT calculations :cite:`Paier2005,Broqvist2009,Sundararaman2013`.  In
PySCF, this exchange divergence can be addressed in a number of ways using the
:attr:`exxdiv` keyword argument to the mean-field class, with the following
possible values.

* ``'ewald'`` (default)

  The :math:`G=0` value of the Coulomb potential is the supercell Madelung
  constant :cite:`Paier2005,Broqvist2009,Sundararaman2013`, which is evaluated by
  Ewald summation.  The finite-size error of the exchange energy decays as
  :math:`N_k^{-1}`, where :math:`N_k` is the number of k-points sampled in the
  Brillouin zone.

* ``None``

  The :math:`G=0` value of the Coulomb potential is set to zero.  The
  finite-size error of the exchange energy decays slowly as :math:`N_k^{-1/3}`.

* ``'vcut_sph'``

  The Coulomb potential is spherically truncated in real space at a radius equal
  to half of the supercell side length :cite:`Spencer2008`.  The finite-size
  error of the exchange energy decays as :math:`\exp(-aN_k^{1/3})`.  Only
  supported with plane-wave density fitting (FFTDF).

* ``'vcut_ws'``

  The Coulomb potential is truncated outside of the Wigner-Seitz supercell
  :cite:`Sundararaman2013`, which is more appropriate than spherical truncation
  for anisotropic cells.  The finite-size error of the exchange energy decays as
  :math:`\exp(-aN_k^{1/3})`.  Only supported with plane-wave density fitting
  (FFTDF).

An example calculation with exchange treated with the spherically truncated
Coulomb potential is shown here::

    kmf = scf.KRHF(cell, kpts=kpts, exxdiv='vcut_sph').run()

Band structure calculations
===========================

After an SCF calculation has been performed, the band structure can be
calculated non-self-consistently along a k-point path using the
``SCF.get_bands(kpts)`` function, where ``kpts`` is a list of k-points along
which the band structure is desired.

.. warning::
   For Hartree-Fock or hybrid DFT, the discontinuity of the exchange potential
   at :math:`G=0` is problematic for band structure calculations.  Using
   :attr:`exxdiv='vcut_sph'` with FFTDF is recommended instead.  Alternatively,
   the SCF procedure can be repeated at each k-point, which is much more
   expensive but allows the use of any :attr:`exxdiv` or density fitting.

See :source:`examples/pbc/09-band_ase.py` for an example DFT band structure
calculation.

Add-ons
=======

All molecular SCF add-ons are also available for periodic SCF but
must be accessed through the molecular :mod:`pyscf.scf.addons` module.
Here we highlight a few of the most useful add-ons.

Linear dependencies
-------------------

The dense packing of atoms in solids combined with the use of diffuse
atom-centered basis functions is responsible for frequent linear dependencies.
The linear dependency can be eliminated by Cholesky orthogonalization::

    from pyscf import scf as mol_scf
    kmf = scf.KRHF(cell, kpts=kpts)
    kmf = mol_scf.addons.remove_linear_dep_(kmf).run()

Smearing
--------

For metals or small band gap semiconductors, it can be useful to smear the
orbital occupation numbers away from integer values.  This can improve SCF
convergence and can expedite convergence to the thermodynamic limit with k-point
sampling.  Because this approach assumes a finite electronic temperature, it
yields an entropy and free energy::

    kmf = scf.KRHF(cell, kpts=kpts)
    kmf = scf.addons.smearing_(kmf, sigma=0.01, method='fermi').run()
    print('Entropy = %s' % kmf.entropy)
    print('Free energy = %s' % kmf.e_free)
    print('Zero temperature energy = %s' % ((kmf.e_tot+kmf.e_free)/2))

Fermi-Dirac smearing (:attr:`method='fermi'`) and Gaussian smearing
(:attr:`method='gauss'`) are supported.

.. warning::
   Because most functions in PySCF assume integer occupations, they may fail if
   combined with a mean-field calculation that was performed with smearing.

Stability
---------

Periodic SCF solutions can be checked with stability analysis::

    kmf = scf.KRHF(cell).run()
    kmf.stability()

.. _user_pbc_multigrid:

Multigrid
---------
For pure DFT calculations, multi-grid algorithm becomes efficient if the size of
unit cell is relatively large (e.g. more than 10 atoms) or the kinetic energy
cutoff is relatively high (e.g. more than 100k plane waves). This algorithm is
implemented in the :mod:`pyscf.pbc.dft.multigrid` module::

  from pyscf.pbc.dft import multigrid
  mf = multigrid.multigrid(dft.KRKS(cell))
  mf.run()

More examples can be found

.. literalinclude:: ../../../examples/pbc/27-multigrid.py

The multi-grid algorithm does not support exact exchange.

.. _user_pbc_rsjk:

Range-separation integration
----------------------------
This algorithm computes most of the four-center integrals in real space. For small
systems, it is less efficient than most density fitting algorithms. When you
need to handle a huge unit cell, or to compute exact exchange for many k points,
or to use all-electron basis sets, you can consider to invoke this algorithm.
For periodic SCF models, you can call the method `jk_method`::

  kmf = scf.KRHF(cell, kpts).jk_method('RS')
  kmf.kernel()

This algorithm has much smaller memory footprint than density fitting algorithm.
It is also a good choice if your computer does not have enough memory to call
density fitting methods.

.. note:: The implementation in current release does not support band structure calculations


.. _user_pbc_df_comparison:

How to choose integral scheme
=============================

In the tables below, we provide a very rough estimation of the capability and
characters for each integral algorithm in typical periodic mean-field
calculations.

1. Number of basis functions for Gamma point calculations

============= ================ ================
Scheme        Pure DFT         hybrid DFT
============= ================ ================
FFTDF         < 5000           < 300
GDF           < 5000           < 5000
MDF           < 5000           < 2000
Multi-grid    100 - 100000     N/A
RS            100 - 10000      100 - 10000
============= ================ ================

2. Assuming 100 AOs per unit cell, the number of k-points

============= ================ ================
Scheme        Pure DFT         hybrid DFT
============= ================ ================
FFTDF         < 2000           < 10
GDF           < 2000           < 100
MDF           < 2000           < 50
Multi-grid    5000             N/A
RS            1000             1000
============= ================ ================

3. Relative performance for systems of 100 AOs per unit cell, 10 k-points

============= ================ ================
Scheme        Pure DFT         hybrid DFT
============= ================ ================
FFTDF         1x               1000x
GDF           10x              100x
MDF           30x              300x
Multi-grid    1x               N/A
RS            100x             100x
============= ================ ================

4. Other metrics

============= ================ ============== ========== =============== ============
Scheme        All-electron     Warm-up cost   Accuracy   Memory cost     IO cost
============= ================ ============== ========== =============== ============
FFTDF         Extremely slow   Low            High       Extremely high  None
GDF           Well supported   High           Low        Moderate        Huge
MDF           Well supported   High           Medium     Moderate        Huge
Multi-grid    Limited cases    Low            High       High            None
RS            Well supported   Low            High       Low             None
============= ================ ============== ========== =============== ============
