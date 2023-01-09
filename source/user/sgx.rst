.. _user_sgx:

Seminumerical exchange (SGX)
****************************

*Modules*: :mod:`sgx`

The SGX module implements seminumerical computation of the exchange matrix.

Introduction
============

Direct computation of the Hartree-Fock exchange matrix in the atomic orbital basis scales poorly with system size.
To achieve better scaling, one three-dimensional integral in the 6-dimensional two-electron integrals can be computed analytically, while the other can be evaluated on a real-space grid, as proposed by Friesner :cite:`Friesner1985`.
The PySCF implementation resembles the chain-of-spheres (COSX) algorithm of Neese et al. :cite:`Neese2009`, but uses more conservative grids and a slightly different P-junction screening function.
Overlap fitting is used to reduce aliasing errors :cite:`Izsak2011`.
SGX scales as :math:`O(N^2)` with system size, as opposed to the :math:`O(N^4)` scaling of analytical exchange :cite:`Neese2009`.

Usage and Example
=================

Any :attr:`scf.hf.SCF` object :attr:`mf` can be converted to an equivalent object that computes the Coulomb and exchange matrices with SGX instead of analytical integration by calling :code:`sgx.sgx_fit(mf)`.

* :source:`examples/sgx/00-simple_sgx.py`

.. literalinclude:: ../../examples/sgx/00-simple_sgx.py

In this case, the error of DFJ+SGX compared to analytical exchange is about 0.03 mEh. The line

.. code::

    mf.with_df.dfj = True

specifies to compute the Coulomb matrix using density fitting (DF-J) while using SGX for the exchange matrix.

The :attr:`pjs` option turns on P-junction screening, which means that the three-center integrals are screened by the density matrix.
If a three-center integral gets contracted only with negligibly small density matrix elements, it is not computed.
This feature can only be used with :attr:`dfj=True`.
If :attr:`dfj=False` and :attr:`pjs=True`, :attr:`dfj` is set to :attr:`True`, and a warning is issued.
The P-junction screening threshold is determined by :attr:`mf.direct_scf_tol`.

Direct evaluation of the J matrix can be used by setting :attr:`mf.with_df.direct_j = True`, but this is much slower than SGX-J or DF-J and only intended for testing purposes.

Adjustable Parameters
=====================

Calling the :attr:`sgx_fit` function on an :attr:`scf.hf.SCF` object returns an equivalent :attr:`SGXHF` object with a :attr:`with_df` attribute that handles SGX integration. To use a non-default auxiliary basis (for :attr:`dfj=True`), :attr:`auxbasis` can be specified in the call to :attr:`sgx_fit`. In addition, there are various adjustable parameters for SGX:

:attr:`with_df` attributes:

* :attr:`grids_level_i`: The grid level to use for initial SCF iterations.
* :attr:`grids_level_f`: The grid level to use for final SCF iterations.
* :attr:`grids_thrd`: The grid points where no atomic orbital is significant (has a value greater than this threshold) are removed from consideration.
* :attr:`grids_switch_thrd`: The threshold for the magnitude of the change in the density matrix that triggers the switch from the initial grid specified by :attr:`grids_level_i` to the final grid specified by :attr:`grids_level_f`.
* :attr:`blockdim`: The maximum number of grid points to loop over at once. The number of grid points per batch is the minimum of :attr:`blockdim` and the maximum number of points allowed by the memory available for the calculation. The maximum memory can be adjusted by setting the :attr:`max_memory` attribute, which is initially set to :attr:`mol.max_memory`, the max memory setting for the Mole object.
* :attr:`dfj`: If :attr:`True`, density fitting is used for the J-matrix. If :attr:`False`, SGX is used.
* :attr:`direct_j`: If :attr:`True`, direct evaluation of the J-matrix is used (slow, for testing only).
* :attr:`pjs`: If :attr:`True`, P-junction screening is used. Threshold determined by :attr:`direct_scf_tol`.

:attr:`SGXHF` attribute:

* :attr:`direct_scf_sgx`: Whether to use direct SCF within the SGX module, meaning that the J and K matrices are evaluated from the difference in the density matrix from the previous iteration.
* :attr:`rebuild_nsteps`: Rebuild the SGX JK matrix from scratch every :attr:`rebuild_nsteps` steps (default 5). Set to 0 to turn off rebuilding the JK matrix (Warning: This can cause creeping numerical error).

References
==========

.. bibliography:: ref_sgx.bib
   :style: unsrt
