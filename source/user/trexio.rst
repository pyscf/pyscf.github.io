.. _user_trexio:

TREXIO support
**************

*Modules*: :py:mod:`pyscf.tools.trexio`

Introduction
============

The `TREXIO <https://trex-coe.github.io/trexio/>`_ file format
is an exchange format for wavefunction data among quantum-chemistry packages
(`PySCF <https://github.com/pyscf/pyscf>`_,
`Quantum Package <https://github.com/QuantumPackage/qp2>`_,
`CP2K <https://github.com/cp2k/cp2k>`_, ...) and downstream consumers, in
particular ab initio Quantum Monte Carlo (QMC) packages such as
`TurboRVB <https://github.com/sissaschool/turborvb>`_,
`CHAMP <https://github.com/filippi-claudia/champ>`_,
`QMC=Chem <https://github.com/TREX-CoE/qmcchem2>`_,
`jQMC <https://github.com/jqmc-project/jQMC>`_, and the
`QMCkl <https://github.com/TREX-CoE/qmckl>`_ kernel library.

For the full list of data that can be stored in a TREXIO file, see the
`TREXIO documentation <https://trex-coe.github.io/trexio/>`_.

Usage
=====

Exporting a PySCF object to a TREXIO file is done by calling :func:`~pyscf.tools.trexio.to_trexio`.
The function dispatches on the type of the object passed in (``Mole``/``Cell``, ``SCF``, or MCSCF)
and writes the appropriate quantities.  The following examples are taken from
:source:`examples/tools/13-trexio.py`.

1. Export only molecular geometry and basis (no SCF)
----------------------------------------------------

Passing a :class:`~pyscf.gto.mole.Mole` (or a :class:`~pyscf.pbc.gto.cell.Cell`)
writes only the nuclear geometry, basis set, ECP, and symmetry information.

.. code-block:: python

    from pyscf import gto
    from pyscf.tools import trexio
    mol = gto.M(atom='H 0 0 0; F 0 0 1.8', basis='cc-pvdz', verbose=0)
    trexio.to_trexio(mol, 'hf_mol.h5')

2. Export SCF results without integrals or density matrices
-----------------------------------------------------------

Passing a converged mean-field object additionally stores the MO coefficients,
orbital energies, occupations, and spin labels.  Integrals and density matrices
can also be written.

.. code-block:: python

    from pyscf import gto, scf
    from pyscf.tools import trexio
    mol = gto.M(atom='H 0 0 0; F 0 0 1.8', basis='cc-pvdz', verbose=0)
    mf = scf.RHF(mol).run()
    trexio.to_trexio(mf, 'hf_scf.h5',
        write_ao_eri=False, write_mo_eri=False, eri_sym='s1', write_mo_rdm=False)

3. Export SCF results with MO integrals and density matrices
------------------------------------------------------------

Setting ``write_mo_eri=True`` and ``write_mo_rdm=True`` writes the MO-basis
one- and two-electron integrals and the one- and two-body reduced density 
matrices.

.. code-block:: python

    from pyscf import gto, scf
    from pyscf.tools import trexio
    mol = gto.M(atom='H 0 0 0; F 0 0 1.8', basis='cc-pvdz', verbose=0)
    mf = scf.RHF(mol).run()
    trexio.to_trexio(
        mf, 'hf_full.h5',
        write_ao_eri=False, write_mo_eri=True, eri_sym='s4',
        write_mo_rdm=True,
    )

4. Export CASSCF results with active-space integrals and density matrices
-------------------------------------------------------------------------

Passing an MCSCF object stores the SCF data plus the CI determinants, natural
occupations, and (by default) the active-space effective integrals and reduced
density matrices.  ``ci_threshold`` discards determinants whose CI coefficient
is below the given magnitude.

.. code-block:: python

    from pyscf import gto, scf, mcscf
    from pyscf.tools import trexio
    mol = gto.M(atom='H 0 0 0; F 0 0 1.8', basis='cc-pvdz', verbose=0)
    mf = scf.RHF(mol).run()
    mc = mcscf.CASSCF(mf, 6, 6).run()
    trexio.to_trexio(
        mc, 'cas.h5',
        write_mcscf_eri=True,
        write_mcscf_rdm=True,
        ci_threshold=1e-3,
    )

Supported objects
=================

The table below summarizes which PySCF objects can be exported with
:func:`~pyscf.tools.trexio.to_trexio` and which quantities are written for each.

.. list-table::
   :header-rows: 1
   :widths: 26 12 12 12 12 12 14
   :stub-columns: 1

   * - PySCF class
     - Geometry / basis / ECP
     - MO coefficients
     - AO integrals
     - MO integrals
     - 1-/2-RDM
     - CI determinants
   * - ``Mole`` / ``Cell``
     - ✓
     - —
     - —
     - —
     - —
     - —
   * - ``RHF/RKS``, ``UHF/UKS``, ``ROHF/ROKS`` (and PBC Γ-point)
     - ✓
     - ✓
     - ○
     - ○
     - ○
     - —
   * - ``CASCI/CASSCF``, ``UCASCI/UCASSCF``
     - ✓
     - ✓
     - —
     - ○
     - ○
     - ✓
   * - ``MP2`` / ``CCSD`` / ``CISD`` and other post-HF
     - —
     - —
     - —
     - —
     - —
     - —

Legend:

* ✓ — always written.
* ○ — optional; controlled by a keyword argument.
* — — not supported.

.. note::

   * Periodic (PBC) calculations must be Γ-point only; k-point objects are not
     supported.
   * Complex-valued integrals are not supported.
   * Post-Hartree-Fock objects (``MP2``, ``CCSD``, ``CISD``, ...) are not yet
     supported.
