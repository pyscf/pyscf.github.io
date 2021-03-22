.. _mcscf:

:mod:`mcscf` --- Multi-configurational self-consistent field
************************************************************

The :mod:`mcscf` implements orbital optimization for
MCSCF and CASSCF. 1-step (combined orbital and wavefunction
optimization) and 2-step algorithms (alternating orbital and wavefunction
optimization) are available. Different kinds of active space solvers can
be used with this module.

For example, a simple CASCI calculation can be run as::

    import pyscf
    
    mol = pyscf.M(
        atom = 'O 0 0 0; O 0 0 1.2',
        basis = 'ccpvdz',
        spin = 2)
    
    myhf = mol.RHF().run()
    
    # 6 orbitals, 8 electrons
    mycas = myhf.CASCI(6, 8).run()

and a simple CASSCF can be run as::
  
   import pyscf

    mol = pyscf.M(
        atom = 'O 0 0 0; O 0 0 1.2',
        basis = 'ccpvdz',
        spin = 2)
    
    myhf = mol.RHF().run()
    
    # 6 orbitals, 8 electrons
    mycas = myhf.CASSCF(6, 8).run()
    
The CASSCF orbital optimization is general and can be combined
with many different solvers, such as DMRG and selected CI solvers. 
Optimized orbitals are stored in the attribute :attr:`mycas.mo_coeff`.

Examples
========

* :source:`examples/mcscf/00-simple_casci.py`
* :source:`examples/mcscf/00-simple_casscf.py`
* :source:`examples/mcscf/01-for_expensive_fci.py`
* :source:`examples/mcscf/03-natural_orbital.py`
* :source:`examples/mcscf/04-density_matrix.py`
* :source:`examples/mcscf/10-define_cas_space.py`
* :source:`examples/mcscf/11-casscf_with_uhf_uks.py`
* :source:`examples/mcscf/12-c2_triplet_from_singlet_hf.py`
* :source:`examples/mcscf/13-load_chkfile.py`
* :source:`examples/mcscf/13-restart.py`
* :source:`examples/mcscf/14-project_init_guess.py`
* :source:`examples/mcscf/15-state_average.py`
* :source:`examples/mcscf/15-state_specific.py`
* :source:`examples/mcscf/15-transition_dm.py`
* :source:`examples/mcscf/16-density_fitting.py`
* :source:`examples/mcscf/17-approx_orbital_hessian.py`
* :source:`examples/mcscf/18-o2_spatial_spin_symmetry.py`
* :source:`examples/mcscf/18-spatial_spin_symmetry.py`
* :source:`examples/mcscf/19-frozen_core.py`
* :source:`examples/mcscf/20-change_symmetry.py`
* :source:`examples/mcscf/21-active_space_symmetry.py`
* :source:`examples/mcscf/21-nosymhf_then_symcasscf.py`
* :source:`examples/mcscf/22-x2c.py`
* :source:`examples/mcscf/23-local_spin.py`
* :source:`examples/mcscf/24-callback.py`
* :source:`examples/mcscf/34-init_guess_localization.py`
* :source:`examples/mcscf/40-customizing_hamiltonian.py`
* :source:`examples/mcscf/41-mcscf_custom_df_hamiltonian.py`
* :source:`examples/mcscf/41-state_average.py`
* :source:`examples/mcscf/41-state_average_custom_fci.py`
* :source:`examples/mcscf/42-compare_cas_space.py`
* :source:`examples/mcscf/43-PiOS-C8H10.py`
* :source:`examples/mcscf/43-avas.py`
* :source:`examples/mcscf/43-dmet_cas.py`
* :source:`examples/mcscf/44-mcscf_active_space_hamiltonian.py`
* :source:`examples/mcscf/50-casscf_then_dmrgscf.py`
* :source:`examples/mcscf/50-casscf_with_selected_ci.py`
* :source:`examples/mcscf/50-cornell_shci_casscf.py`
* :source:`examples/mcscf/50-dmrgscf_with_block.py`
* :source:`examples/mcscf/51-o2_triplet_by_various_fci.py`
* :source:`examples/mcscf/60-uhf_based_ucasscf.py`
* :source:`examples/mcscf/70-casscf_optimize_scheduler.py`

CASSCF active space solver
==========================
.. automodule:: pyscf.mcscf

DMRG solver
-----------

FCIQMC solver
-------------

State-average FCI solver
------------------------

State-average with mixed solver
-------------------------------


Symmetry broken
===============
.. note OMP threads cause numerical instability and break spin symmetry


Initial guess
=============


Canonical orbitals
==================
Orbital canonicalization are controlled by parameters :attr:`mc.canonicalization`
and :attr:`mc.natorb` (assuming the MCSCF object is ``mc``). The order of
orbitals are affected by the parameter :attr:`mc.sorting_mo_energy`.

* canonicalization:
  This flag canonicalizes orbitals in core and external space using general Fock matrix.
* natorb:
  Transforms active orbitals using 1-particle density matrices.
* sorting_mo_energy:
  Sort orbitals based on the diagonal elements of the general Fock matrix.
* Enabling natorb or sorting_mo_energy may slightly change the total energy
  of DMRG solver or selected CI solver.

General Fock matrix is defined as

.. math::
  \mathbf{F} &= \mathbf{h}_{core} + \mathbf{J} - \mathbf{K} \\
  J_{pq} &= \sum_{rs} (pq|rs) \gamma_{sr} \\
  K_{pq} &= \sum_{qr} (pq|rs) \gamma_{qr} \\

:math:`\gamma` is the total density matrix which includes the doubly occupied
core density matrix and correlated density matrix in active space.

If :attr:`mc.canonicalization` is enabled, CASCI/CASSCF will call
the :func:`mc.canonicalize` function to diagonalize orbitals in
**core space and external space**. Orbitals in active space are not changed if
merely setting :attr:`mc.canonicalization`. In the attribute
:attr:`mc.mo_energy`, eigenvalues of general Fock matrix for core and external
subspaces are stored in the corresponding sub-sectors. The sector of active
space in :attr:`mc.mo_energy` stores the expectation value of general Fock
:math:`\langle \phi|F|\phi\rangle`.
By default, :attr:`mc.canonicalization` is enabled because the canonicalized
MCSCF orbitals can simplify the implementations of MRPT (NEVPT2) methods.  

:attr:`mc.natorb` controls whether the CASCI/CASSCF active space orbitals are
transformed to natural orbitals w.r.t. the correlated density matrix.  When this
parameter is enabled, the natural orbitals will be stored in the active part of
the attribute :attr:`mc.mo_coeff` and the CI coefficients :attr:`mc.ci` (if
applicable) will be transformed accordingly.
By default :attr:`mc.natorb` is disabled because natural orbitals may not be
favored by total energy for an arbitrary CI solver. We make this default value
to ensure that a CASCI calculation followed by a CASSCF calculation (e.g.
DMRG-CASSCF then DMRG-CASCI) produces results same to the CASSCF results.
The CASCI calculation may produce different
The value of :attr:`mc.natorb` does not affect (the default) FCI solver. But
this is not true for external large active space solvers such
as DMRG, selected CI methods. It is recommended to disable :attr:`mc.natorb` in
these calculations.

Following presents what the :attr:`mc.mo_coeff` would be like for different
combinations of :attr:`mc.canonicalization` and :attr:`mc.natorb` in a
CASCI calculation:

* ``mc.canonicalization = False`` and ``mc.natorb = False``:
The resultant orbitals :attr:`mc.mo_coeff` are identical to the input orbitals.
If the CASCI was initialized with a RHF calculation,  :attr:`mc.mo_coeff` points
to RHF orbitals.
 
* ``mc.canonicalization = True`` and ``mc.natorb = False``:
Core part and external part of :attr:`mc.mo_coeff` are canonicalized orbitals,
which diagonalize the core and external blocks of general Fock matrix.  The
orbitals in active space are identical to the active orbitals in the input.
 
* ``mc.canonicalization = False`` and ``mc.natorb = True``
Core and external part of :attr:`mc.mo_coeff` are identical to the core and
external part of the input orbitals.  Active space orbitals are transformed to
the natural orbitals of the correlated density matrix.

* ``mc.canonicalization = True`` and ``mc.natorb = True``
:attr:`mc.mo_coeff` are completely different to the input orbitals.

Please note that elements of :attr:`mc.mo_energy` may not be sorted ascendantly.
Parameter :attr:`mc.sorting_mo_energy` can affect the ordering of MCSCF orbitals
**when :attr:`mc.canonicalization` or :attr:`mc.natorb` is enabled**.

By default, canonical orbitals in the core and external space
are sorted by the orbital energies (from low to high) and the natural orbitals
in the active space are sorted by natural occupations (from large to small).
If point group symmetry is enabled in the calculation, canonical orbitals are
sorted within each symmetry sector only (rather than the entire core or external
space). Irreducible representation labels (can be accessed via
:attr:`mc.mo_coeff.orbsym`) are assigned to orbitals in the initial
guess and they will not be changed in the MCSCF optimization and the subsequent
canonicalization procedure.  Setting :attr:`mc.sorting_mo_energy` (though not
recommended) can force the orbitals to be sorted against energy (or occupations
in active space) regardless whether the point group symmetry is used.

In certain scenario, you may want to enable both :attr:`mc.natorb` and
:attr:`mc.sorting_mo_energy`.
``examples/dmrg/31-cr2_scan/cr2-scan.py`` provides one example that needs both
parameters. In that example, the dissociation curve of Cr dimer
was scanned using heat-bath selected-CI method in which the active space of
selected-CI-CASSCF was gradually enlarged in a series of CASSCF calculations.
Since the selected-CI algorithm depends on the initial single determinant, the
orbital ordering does matter to the final CASSCF results.
:attr:`mc.natorb` and :attr:`mc.sorting_mo_energy` have to be enabled to make
sure that the each selected-CI starts from the similar initial reference at
each point on the dissociation curve. Without these settings, the differences
in the orbital ordering can lead to discontinuous potential energy curve.


Program reference
=================

CASCI
-----

.. automodule:: pyscf.mcscf.casci
   :members:

.. automodule:: pyscf.mcscf.casci_symm
   :members:

.. automodule:: pyscf.mcscf.ucasci
   :members:


CASSCF
------

.. automodule:: pyscf.mcscf.mc1step
   :members:

.. automodule:: pyscf.mcscf.mc1step_symm
   :members:

.. automodule:: pyscf.mcscf.umc1step
   :members:

.. automodule:: pyscf.mcscf.mc_ao2mo
   :members:

.. automodule:: pyscf.mcscf.umc_ao2mo
   :members:


addons
------

.. automodule:: pyscf.mcscf.addons
   :members:

