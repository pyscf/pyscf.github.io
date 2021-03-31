.. _user_sgx:

Seminumerical Exchange (SGX)
****************************

*Modules*: :mod:`sgx`

The SGX module implements seminumerical computation of the exchange matrix.

Introduction
============

Direct computation of the Hartree-Fock Exchange matrix in the atomic orbital basis scales poorly with system size. To achieve better scaling, one coordinate of the two-electron integrals is computed analytically, while the other is evaluated on a real-space grid, as proposed by Friesner :cite:`Friesner1985`. The PySCF implementation most closely resembles the chain-of-spheres (COSX) algorithm of Neese et al. :cite:`Neese2009`, but with more conservative grids and without P-junction screening. Overlap fitting is used to reduce aliasing errors :cite:`Izsak2011`. SGX scales as :math:`O(N^2)` with system size, as opposed to the :math:`O(N^4)` scaling of analytical exchange :cite:`Neese2009`.

Usage and Example
=================

Any :attr:`scf.hf.SCF` object :attr:`mf` can be converted to an equivalent object that computes the Coulomb and Exchange matrices with SGX instead of analytical integration by calling :code:`sgx.sgx_fit(mf)`.

* :source:`examples/sgx/00-simple_sgx.py`

.. literalinclude:: ../../examples/sgx/00-simple_sgx.py

.. code::

  converged SCF energy = -76.0267374704185
  converged SCF energy = -76.0267978618974

In this case, the error of SGX compared to analytical exchange is about 0.06 mEh. The line

.. code::

    mf.with_df.dfj = True

specifies to compute the Coulomb matrix using Density Fitting (DF) while using SGX for the Exchange matrix.

Adjustable Parameters
=====================

Calling the :attr:`sgx_fit` function on an :attr:`scf.hf.SCF` object returns an equivalent :attr:`SGXHF` object with a :attr:`with_df` attribute that handles SGX integration. To use a non-default auxiliary basis (for :attr:`dfj=True`), :attr:`auxbasis` can be specified in the call to :attr:`sgx_fit`. In addition, there are five main adjustable parameters for this object:

* :attr:`grids_level_i`: The grid level to use for initial SCF iterations.
* :attr:`grids_level_f`: The grid level to use for final SCF iterations.
* :attr:`grids_thrd`: The grid points at which all atomic orbitals have a value below this threshold are removed from the integration grid.
* :attr:`grids_switch_thrd`: The threshold for the magnitude of the change in the density matrix that triggers the switch from the initial grid specified by :attr:`grids_level_i` to the final grid specified by :attr:`grids_level_f`.
* :attr:`blockdim`: The maximum number of grid points to loop over at once. The number of grid points per batch is the minimum of :attr:`blockmin` and the maximum number of points allowed by the memory available for the calculation. The maximum memory can be adjusted by setting the :attr:`max_memory` attribute, which is initially set to :attr:`mol.max_memory`, the max memory setting for the Mole object.

References
==========

.. bibliography:: ref_sgx.bib
   :style: unsrt
