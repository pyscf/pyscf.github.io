.. _user_lo:

******************
Localized orbitals
******************

*Modules*: :mod:`lo`

Introduction
============
A molecular orbital is usually delocalized, i.e. it has non-negligible amplitude over the whole system rather than only around some atom(s) or bond(s).
However, one can choose a unitary rotation :math:`U`

.. math::

    \phi = \psi U

such that the resulting
orbitals :math:`\phi` are as spatially localized as possible. This is typically achieved by one of two classes of
methods. 
The first is to project the orbitals onto a predefined local set of orbitals, which can be e.g. atomic orbitals or pseudo-atomic orbitals.
The second is to optimize a cost function :math:`f`, which measures the locality of the molecular orbitals.
Because there is no unambiguous choice for the localization criterion, several criteria have been suggested.
Boys localization minimizes the spread of the orbital

.. math::

    f(U) = \sum_{i} \langle\psi_i|r^2|\psi_i\rangle - \langle\psi_i|r|\psi_i\rangle^2

Boys localized orbitals :cite:`Foster60` in periodic systems are typically termed maximally localized
Wannier orbitals (MLWF) :cite:`Marzari97`.

Pipek-Mezey (PM) localization :cite:`Pipek98` maximizes the population charges on the atoms

.. math::

    f(U) = \sum^{\mathrm{atoms}}_{I} \sum_{i} \left|q^{I}_{i} \right|^2

Note that PM localization depends on the choice of atomic orbitals used for
the population analysis. Several choices of populations are available, e.g. 
Mulliken or based on (meta-) L\"owdin orbitals. 
Intrinsic bond orbitals (IBOs) can be viewed as a special case of PM
localization using intrinsic atomic orbitals (IAOs) as population method.
See Ref. :cite:`Lehtola14PM` for a summary of choices of orbitals. 
Note that PM localization preserves the separation between  :math:`\sigma` and :math:`\pi` orbitals.

Edmiston-Ruedenberg (ER) localization :cite:`Edmiston63` maximizes the orbital Coulomb self-repulsion,

.. math::

    f(U) = \sum_{i} (ii|ii)

ER localization, however, is computationally more expensive than the Boys or PM approaches.

Localized orbitals can be calculated via the pivoted Cholesky factorization of a density-like
matrix :math:`\mathbf{D} = \mathbf{C} \mathbf{C}^\dagger`. :cite:`Aquilante06` Since :math:`\mathbf{C}` is 
generally a rectangular matrix containing only the subset of :math:`N` orbitals intended for localization,
the matrix :math:`\mathbf{D}` is positive-semidefinite. It can be factored using a Cholesky decomposition
with full column pivoting,

.. math::
    \mathbf{P}^\dagger \mathbf{D} \mathbf{P} = \mathbf{L} \mathbf{L}^\dagger ,

where :math:`\mathbf{L}` is a lower triangular matrix and :math:`\mathbf{P}` is a permutation matrix.
In the end, the :math:`N` leftmost columns of :math:`\mathbf{P L}` are taken as the localized orbitals.
While Cholesky orbitals are usually not as localized as, for example, PM or Boys orbitals, the procedure
is non-iterative and produces unique result, except possibly for the impact of degeneracies.
Cholesky orbitals can serve as an excellent guess for iterative localization procedures.


A summary of the functionality of the :mod:`lo` module is given below:

=========================== ============== ==================== ======== =====
Method                       optimization   cost function        PBC     ref
(meta-) L\"owdin                 No            -                 yes     :cite:`Lowdin50,Sun14qmmm`
Natural atomic orbitals          No            -                 gamma   :cite:`Reed85` 
Intrinsic atomic orbitals        No            -                 yes     :cite:`Knizia13IAO`
Cholesky orbitals                No            -                 no      :cite:`Aquilante06`
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
