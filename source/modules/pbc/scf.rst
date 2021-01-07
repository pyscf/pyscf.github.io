.. _pbc_scf:

:mod:`pbc.scf` --- Self-consistent field with periodic boundary conditions
**************************************************************************

.. module:: pbc.scf
   :synopsis: Hartree-Fock methods with periodic boundary conditions

This module is constructed to be analogous to the molecular :mod:`scf` module to handle
mean-field calculations with periodic boundary conditions.

Gamma point and single k-point calculation
==========================================
The usage of the Gamma point Hartree-Fock program is very close to that of the
molecular program.  In a PBC Gamma point calculation, one needs to initialize the
:class:`Cell` object and the corresponding :class:`pyscf.pbc.scf.hf.RHF` class::

    from pyscf.pbc import gto, scf
    cell = gto.M(
        atom = '''H     0.      0.      0.    
                  H     0.8917  0.8917  0.8917''',
        basis = 'sto3g',
        h = '''
        0       1.7834  1.7834
        1.7834  0       1.7834
        1.7834  1.7834  0     ''',
        gs = [10]*3,
        verbose = 4,
    )
    mf = scf.RHF(cell).run()

Compared to the :class:`pyscf.scf.hf.RHF` object for molecular calculations,
the PBC-HF calculation with :class:`pyscf.pbc.scf.hf.RHF` or
:class:`pyscf.pbc.scf.uhf.UHF` has three differences

* :class:`pyscf.pbc.scf.hf.RHF` is the single k-point PBC HF class.  By default,
  it creates a Gamma point instance.  You can change to other (single) k-points by
  setting the :attr:`kpt` attribute::

    mf = scf.RHF(cell)
    mf.kpt = cell.get_abs_kpts([.25,.25,.25])  # convert from scaled kpts
    mf.kernel()

* The exchange integrals of the PBC Hartree-Fock method show a slow convergence
  with respect to the number of k-points.  A proper treatment for the divergent
  part of the exchange integral can improve the convergence.  The attribute
  :attr:`exxdiv` is used to control the method to handle the divergent exchange
  term.  The default ``exxdiv='ewald'`` can be used for most scenarios.  However,
  if molecular post-HF methods are mixed with the Gamma point mean-field methods (see
  :ref:`mix_to_mol`, you will need to explicitly use ``exxdiv=None`` to obtain a consistent total
  energy (see :ref:`exxdiv`). If PBC post-HF methods are used, the :mod:`exxdiv` attribute
  is always treated as ``None`` in the correlated calculation to obtain a consistent total energy.

* In PBC calculations there are different choices for how to 
  evaluate 2-electron integrals.  The default integral scheme (FFTDF) is accurate   for pseudo-potentials and fast for DFT calculations but
  slow in other scenarios. A second integral
  scheme which is a good balance of speed and accuracy for
  Hartree-Fock, all-electron, and post-Hartree-Fock calculations, is
  Gaussian density fitting. This can be used by setting the :attr:`with_df` attribute (see :ref:`pbc_df`) or, conveniently, by using the :func:`density_fit` function (see :ref:`pbc_density_fit`). For
higher accuracy, you may wish to use
  mixed density fitting (MDF) (see :ref:`with_df`).  Here is an example to update :attr:`with_df`

.. literalinclude:: /../examples/pbc/11-gamma_point_all_electron_scf.py

  

.. _mix_to_mol:

Mixing with molecular program for post-HF methods
-------------------------------------------------
The Gamma point HF code adopts the same code structure, function and
method names and argument' conventions as the molecular SCF code.
This desgin allows one to mix PBC HF objects with the existing molecular post-HF
code for PBC electron correlation treatments.  A typical molecular post-HF
calculation starts from the finite-size HF method with the :class:`Mole`
object::

    from pyscf import gto, scf
    mol = gto.M(atom='H 0 0 0; H 0 0 1', basis='ccpvdz')
    mf = scf.RHF(mol).run()

    from pyscf import cc
    cc.CCSD(mf).run()

The PBC Gamma point post-HF calculation requires the :class:`Cell` object and
PBC HF object::

    from pyscf.pbc import gto, scf
    cell = gto.M(atom='H 0 0 0; H 0 0 1', basis='ccpvdz',
                 h=numpy.eye(3)*2, gs=[10,10,10])
    mf = scf.RHF(cell).run()

    from pyscf import cc
    cc.CCSD(mf).run()

The differences are the ``mol`` or ``cell`` object to create and the
``scf`` module to import.  With the system-specific mean-field object, one
can carry out various post-HF methods (MP2, Coupled cluster, CISD, TDHF,
TDDFT, ...) using the same code for finite-size and extended systems.
See :ref:`mix_mol` for more details of the interface between PBC and molecular
modules.

k-point sampling
================

An example demonstrating the use of k-points for a KRHF, KRKS calculation, and using the Newton (second-order SCF) solver.

.. literalinclude:: /../examples/pbc/20-k_points_scf.py

Smearing
--------

In many periodic systems, the HOMO-LUMO gap may become very small. In such a case, one needs to use occupation number smearing
to converge the SCF calculation, as illustrated here::

    import numpy
    from pyscf.pbc import gto, scf
    
    cell = gto.Cell()
    cell.atom = '''
    He 0 0 1
    He 1 0 1
    '''
    cell.basis = 'ccpvdz'
    cell.a = numpy.eye(3) * 4
    cell.verbose = 4
    cell.build()
    
    #
    # Use scf.addons.smearing_ function to modify the
    # PBC (gamma-point or k-points) SCF object
    #
    nks = [2,1,1]
    mf = scf.KRHF(cell, cell.make_kpts(nks))
    mf = scf.addons.smearing_(mf, sigma=.1, method='fermi')
    mf.kernel()
    

.. _exxdiv:

Exchange divergence treatment
=============================
The attribute :attr:`exxdiv` controls the handling of the slow convergence
of the HF exchange integrals. 
Generally, :attr:`exxdiv` leads to a shift in the total energy that
can be thought of as a HF finite size-correction. It also modifies the spectrum of
the occupied orbitals. When a Gamma-point PBC mean-field calculation is mixed with a molecular post-HF implementation, setting this attribute to anything other than ``None`` will lead to an inconsistency in the total energy.

Possible values are ``None``, ``vcut_sph`` (spherical cutoff), ``vcut_ws`` (Wigner-Seitz cutoff), ``ewald`` (probe-charge Ewald correction).
``vcut_sph`` and ``vcut_ws`` are only available when using FFTDF. When using hybrid functionals, the exchange treatment will use the :attr:`exxdiv` attribute.

Note that when calling :func:`get_bands`, the choice of :attr:`exxdiv` affects the band structure. In particular, ``vcut_sph`` (spherical cutoff), ``vcut_ws`` (Wigner-Seitz cutoff) should be used to obtain smooth bands, otherwise, the band structure should be computed by twisting the
SCF calculation (i.e. redoing self-consistency at each k point).


.. _pbc_density_fit

Gaussian density fitting
========================
This example uses the :func:`density_fit` function to enable Gaussian density fitting::

    import numpy as np
    from pyscf import gto as mol_gto
    from pyscf.pbc import gto, scf, cc, df
    
    cell = gto.Cell()
    cell.atom='''
    C 0.000000000000   0.000000000000   0.000000000000
    C 1.685068664391   1.685068664391   1.685068664391
    '''
    cell.basis = 'gth-dzv'
    cell.pseudo = 'gth-pade'
    cell.a = '''
    0.000000000, 3.370137329, 3.370137329
    3.370137329, 0.000000000, 3.370137329
    3.370137329, 3.370137329, 0.000000000'''
    cell.unit = 'B'
    cell.verbose = 5
    cell.build()
    
    #
    # Default DF auxiliary basis is a set of even-tempered gaussian basis (with
    # exponents alpha * beta**i, i = 1,..,N).  The even-tempered parameter alpha
    # is determined automatically based on the orbital basis.  beta is set to 2.0
    #
    mf = scf.RHF(cell).density_fit()
    mf.kernel()

Program reference
=================

.. automodule:: pyscf.pbc.scf.hf
   :members:

.. automodule:: pyscf.pbc.scf.uhf
   :members:

.. automodule:: pyscf.pbc.scf.khf
   :members:

.. automodule:: pyscf.pbc.scf.kuhf
   :members:

