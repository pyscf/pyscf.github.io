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
orbitals :math:`\phi` are as local as possible. This is typically achieved by two kinds of
methods: One is by projection onto a pre-defined local orbital set
(e.g., AO, pseudopotential orbitals); The other is by numerical optimization on
a cost function :math:`f`, which can be chosen as different localization
criterion, e.g., in Boys localization, the dipole (spread) of a molecule,

.. math::

    f(U) = \langle\psi|r^2|\psi\rangle - \langle\psi|r|\psi\rangle^2


The :mod:`lo` module implements various orbital localizations, including:

=========================== ============== ==================== ======== =====
Method                       optimization   cost function        PBC     ref
(meta-) Lowdin                   No            -                 yes     :cite:`Lowdin50,Sun14qmmm`
Natural atomic orbitals          No            -                 no      :cite:`Reed85` 
Intrinsic Atomic orbitals        No            -                 yes     :cite:`Knizia13IAO`
Intrinsic Bond orbitals          yes         IAO charges         gamma   :cite:`Knizia13IAO`
Boys                             yes         dipole              no      :cite:`Foster60`
Pipek-Mezey                      yes         local charges       gamma   :cite:`Pipek98`
Edmiston-Ruedenberg              yes         coloumb integral    gamma   :cite:`Edmiston63`
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
