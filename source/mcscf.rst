mcscf --- Multi-configurational self-consistent field
*****************************************************

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

Examples
========

Relevant examples
:file:`examples/mcscf/00-simple_casci.py`
:file:`examples/mcscf/00-simple_casscf.py`
:file:`examples/mcscf/01-for_expensive_fci.py`
:file:`examples/mcscf/03-natural_orbital.py`
:file:`examples/mcscf/04-density_matrix.py`
:file:`examples/mcscf/10-define_cas_space.py`
:file:`examples/mcscf/11-casscf_with_uhf_uks.py`
:file:`examples/mcscf/12-c2_triplet_from_singlet_hf.py`
:file:`examples/mcscf/13-load_chkfile.py`
:file:`examples/mcscf/13-restart.py`
:file:`examples/mcscf/14-project_init_guess.py`
:file:`examples/mcscf/15-state_average.py`
:file:`examples/mcscf/15-state_specific.py`
:file:`examples/mcscf/15-transition_dm.py`
:file:`examples/mcscf/16-density_fitting.py`
:file:`examples/mcscf/17-approx_orbital_hessian.py`
:file:`examples/mcscf/18-o2_spatial_spin_symmetry.py`
:file:`examples/mcscf/18-spatial_spin_symmetry.py`
:file:`examples/mcscf/19-frozen_core.py`
:file:`examples/mcscf/20-change_symmetry.py`
:file:`examples/mcscf/21-active_space_symmetry.py`
:file:`examples/mcscf/21-nosymhf_then_symcasscf.py`
:file:`examples/mcscf/22-x2c.py`
:file:`examples/mcscf/23-local_spin.py`
:file:`examples/mcscf/33-make_init_guess`
:file:`examples/mcscf/34-init_guess_localization.py`
:file:`examples/mcscf/40-customizing_hamiltonian.py`
:file:`examples/mcscf/41-mcscf_custom_df_hamiltonian.py`
:file:`examples/mcscf/41-state_average.py`
:file:`examples/mcscf/42-compare_cas_space.py`
:file:`examples/mcscf/43-avas.py`
:file:`examples/mcscf/43-dmet_cas.py`
:file:`examples/mcscf/44-mcscf_active_space_hamiltonian.py`
:file:`examples/mcscf/50-casscf_then_dmrgscf.py`
:file:`examples/mcscf/50-casscf_with_selected_ci.py`
:file:`examples/mcscf/50-cornell_shci_casscf.py`
:file:`examples/mcscf/50-dmrgscf_with_block.py`
:file:`examples/mcscf/51-o2_triplet_by_various_fci.py`
:file:`examples/mcscf/60-uhf_based_ucasscf.py`
:file:`examples/mcscf/61-rcas_vs_ucas`
:file:`examples/mcscf/70-casscf_hot_tuning.py`
:file:`examples/mcscf/70-casscf_optimize_scheduler.py`



.. automodule:: pyscf.mcscf



CASSCF active space solver
==========================

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
There are two relevant parameters for orbital canonicalization.  They are
:attr:`mc.canonicalization` and :attr:`mc.natorb` (assuming the MCSCF object is ``mc``).
In the CASCI/CASSCF calculations, the resultant orbitals are stored in the
attribute :attr:`mc.mo_coeff`.  These orbitals may be identical or partially
identical to the initial orbitals, depending on the values of
:attr:`mc.canonicalization` and :attr:`mc.natorb`.

:attr:`mc.canonicalization` controls whether the resultant CASCI/CASSCF orbitals
are canonicalized with respect to the general Fock matrix.  General Fock matrix
is defined as

.. math::
  \mathbf{F} &= \mathbf{h}_{core} + \mathbf{J} - \mathbf{K} \\
  J_{pq} &= \sum_{rs} (pq|rs) \gamma_{sr} \\
  K_{pq} &= \sum_{qr} (pq|rs) \gamma_{qr} \\

:math:`\gamma` is the total density matrix which is the summation of doubly
occupied core density matrix and correlated density matrix in active space.
If :attr:`mc.canonicalization` is enabled, the CASCI/CASSCF program will call
the :func:`mc.canonicalize` function to diagonalize the core and external space
wrt the general Fock matrix.  The eigenvalues of the core and external subspace
are stored in attribute :attr:`mc.mo_energy`.
By default, :attr:`mc.canonicalization` is enabled because the canonicalized
MCSCF orbitals can simplify the algorithm of MRPT methods.  

:attr:`mc.natorb` controls whether the CASCI/CASSCF active space orbitals are
transformed to natural orbitals of the correlated density matrix.  When this
parameter is enabled, the natural orbitals will be stored in the active part of
the attribute :attr:`mc.mo_coeff` and the CI coefficients :attr:`mc.ci` (if
applicable) will be transformed accordingly.
By default :attr:`mc.natorb` is disabled and it is important for the MCSCF solver.
Generally, the value of :attr:`mc.natorb` does not affect (the default) FCI
solver because an independent CASCI calculation following a previous MCSCF
calculation should give the same solutions no matter :attr:`mc.natorb` is enabled
or not.  But this is not true for some external large active space solvers such
as DMRG, selected CI methods.  The CASCI calculation may produce different
answers depending on the value of :attr:`mc.natorb`.  Therefore, it is
recommended to disable :attr:`mc.natorb` in your calculation.

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

There is another parameter :attr:`mc.sorting_mo_energy` which may affect the
ordering of MCSCF orbitals when :attr:`mc.canonicalization` or :attr:`mc.natorb`
is enabled.  Generally, the canonical orbitals in the core and external space
are sorted by the orbital energies (from low to high) and the natural orbitals
in the active space are sorted by natural occupations (from large to small).
This ordering may not be held if point group symmetry is enabled in the
calculation.  When a system has high spatial symmetry and point group symmetry
is enabled, each SCF orbital will be assigned to an irreducible representation
label.  In the MCSCF calculation and the canonicalization, the irreducible
representation label of the orbitals will not be changed.  They are always the
same to the symmetry labels of the input orbitals.  Although the orbitals
are still ordered within each irreducible representation, the orbital energies
(or occupancies) for all orbitals are not strictly sorted.  Setting
``mc.sorting_mo_energy = Trye`` (though not recommended) can force the orbitals
to be sorted regardless whether the point group symmetry is enabled.  In certain
scenario, you may want to enable :attr:`mc.natorb` and :attr:`mc.sorting_mo_energy`.
``examples/dmrg/31-cr2_scan/cr2-scan.py`` provides one example that you need to
enable the two parameters.  In that example, the dissociation curve of Cr dimer
was scanned by heat-bath selected-CI method in which the active space of
selected-CI-CASSCF was gradually enlarged in a series of CASSCF calculations.
Since the selected-CI algorithm depends on the initial single determinant, the
orbital ordering do have matters to the final CASSCF results. Thus
:attr:`mc.natorb` and :attr:`mc.sorting_mo_energy` have to be enabled to make
sure that the each selected-CI starts from the similar initial reference for
each point on the dissociation curve.  At some critical points, the difference
in the orbital ordering in the active space can lead to discontinuous potential
energy curve.
 

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

