.. _user_lo:

**********************
Localized Orbitals
**********************

*Modules*: :mod:`lo`

Introduction
============
A molecular orbital is usually delocal, i.e. it has non-negligible amplitude all over the whole system rather than only around a specific region localized around some atom(s) or bond(s).
However, one can choose a proper unitary rotation :math:`U`

.. math::

    \phi = \psi U

such that the resulting
orbitals :math:`\phi` are as spatially localized as possible. This is typically achieved by one of two classes of
methods. 
The first is to project the orbitals onto a predefined local set of orbitals, which can be e.g. real atomic orbitals or pseudo-atomic orbitals.
The second is to optimize a cost function :math:`f`, which measures the locality of the molecular orbitals.
Because there is no unambiguous choice for the localization criterion, several criteria have been suggested.
Boys localization minimizes the spread of the orbital

.. math::

    f(U) = \sum_{i} \langle\psi_i|r^2|\psi_i\rangle - \langle\psi_i|r|\psi_i\rangle^2

Boys localization :cite:`Foster60` in periodic systems is typically termed as maximally localized
wannier orbitals (MLWF) :cite:`Marzari97`.

Pipek-Mezey (PM) localization :cite:`Pipek98` maximizes the orbital charges on the atoms

.. math::

    f(U) = \sum^{\mathrm{atoms}}_{I} \sum_{i} \left|q^{I}_{i} \right|^2

Note that the PM localization depends on the choice of orbitals used for
population analysis. Several choices of populations are available, e.g. 
Mulliken, (meta-) L\"owdin orbitals. 
Intrinsic bond orbitals (IBOs) can be viewed as a special case of PM
localization using intrinsic atomic orbitals (IAOs) as population method.
See Ref. :cite:`Lehtola14PM` for a summary of choices of orbitals. 
PM localization has a good feature that enables the
distinguishment of :math:`\sigma` and :math:`\pi` orbitals.

Edmiston-Ruedenberg (ER) localization :cite:`Edmiston63` maximizes the summation of orbital coulomb integrals,

.. math::

    f(U) = \sum_{i} (ii|ii)

ER localization, however, is computationally more expensive than Boys or PM approaches.


The :mod:`lo` module implements various orbital localizations, including:

=========================== ============== ==================== ======== =====
Method                       optimization   cost function        PBC     ref
(meta-) L\"owdin                 No            -                 yes     :cite:`Lowdin50,Sun14qmmm`
Natural atomic orbitals          No            -                 gamma   :cite:`Reed85` 
Intrinsic atomic orbitals        No            -                 yes     :cite:`Knizia13IAO`
Boys                             yes         dipole              no      :cite:`Foster60`
Pipek-Mezey                      yes         local charges       gamma   :cite:`Pipek98`
Intrinsic bond orbitals          yes         IAO charges         gamma   :cite:`Knizia13IAO`
Edmiston-Ruedenberg              yes         coulomb integral    gamma   :cite:`Edmiston63`
=========================== ============== ==================== ======== =====

For example, to obtain the natural atomic orbital coefficients (in terms
of the original atomic orbitals)::

    import numpy
    from pyscf import gto, scf, lo
    
    x = .63
    mol = gto.M(atom=[['C', (0, 0, 0)],
                      ['H', (x ,  x,  x)],
                      ['H', (-x, -x,  x)],
                      ['H', (-x,  x, -x)],
                      ['H', ( x, -x, -x)]],
                basis='ccpvtz')
    mf = scf.RHF(mol).run()
    
    # C matrix stores the AO to localized orbital coefficients
    C = lo.orth_ao(mf, 'nao')

References
==========

.. bibliography:: ref_lo.bib
  :style: unsrt
