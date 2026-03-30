.. _user_sgx:

Seminumerical exchange (SGX)
****************************

*Modules*: :py:mod:`pyscf.sgx`

The SGX module implements seminumerical computation of the exchange matrix (and also the Coulomb matrix).

Introduction
============

Direct computation of the Hartree-Fock exchange matrix in the atomic orbital basis scales poorly with system size.
To achieve better scaling, one three-dimensional integral in the 6-dimensional two-electron integrals can be computed analytically, while the other can be evaluated on a real-space grid, as proposed by Friesner :cite:`Friesner1985`.
The PySCF implementation resembles the chain-of-spheres (COSX) algorithm of Neese et al. :cite:`Neese2009`, but uses slightly different default grids and a modified P-junction screening algorithm.
Overlap fitting is used to reduce aliasing errors :cite:`Izsak2011`.
Before screening of negligible integrals, SGX scales as :math:`O(N^3)` with system size, as opposed to the :math:`O(N^4)` scaling of analytical exchange :cite:`Neese2009`. This makes its performance favorable for large basis sets and large systems.
In principle, both SGX and analytical exchange scale asymptotically as :math:`O(N^2)` if negligible integrals are screened out, and as :math:`O(N)` if the system has a gap and negligible density matrix components are screened.
However, the linear scaling is difficult to achieve in practice and only seen for some very large systems.

Usage and Example
=================

Any :mod:`scf.hf.SCF` object :code:`mf` can be converted to an equivalent object that computes the Coulomb and exchange matrices with SGX instead of analytical integration by calling :code:`sgx.sgx_fit(mf)`.

* :source:`examples/sgx/00-simple_sgx.py`

.. literalinclude:: ../../examples/sgx/00-simple_sgx.py

In this case, the error of DF-J+SGX-K compared to analytical exchange is about 0.03 mEh. The line

.. code::

    mf.with_df.dfj = True

specifies to compute the Coulomb matrix using density fitting (DF-J) while using SGX for the exchange matrix.

Passing :code:`pjs=True` to :code:`sgx_fit` sets :code:`dfj=True` and also turns on P-junction screening (:code:`optk=True`), which means that the three-center integrals are screened by the density matrix.
If a three-center integral gets contracted only with negligibly small density matrix elements, it is not computed. Because the J-matrix cannot be screened in this manner, :code:`optk=True` only activates P-junction screening when used with :code:`dfj=True`.
If :code:`dfj=False`, :code:`optk` is ignored.

When :code:`optk=True`, the P-junction screening threshold is determined by two parameters: :code:`mf.with_df.sgx_tol_energy` and :code:`mf.with_df.sgx_tol_potential`.
These parameters correspond to approximate upper bounds on the error of the exchange energy and K-matrix due to P-junction screening, respectively.
By default, :code:`sgx_tol_energy=sgx_tol_potential="auto"`, in which case :code:`sgx_tol_energy` is set to :code:`mf.direct_scf_tol` and :code:`sgx_tol_potential=sgx_tol_energy**0.5`.
For more aggressive screening that is still numerically stable and sufficiently precise for most cases, we recommend :code:`sgx_tol_energy=mf.conv_tol` or :code:`0.1*mf.conv_tol` and :code:`sgx_tol_potential="auto"`.
The first thing to try if you have numerical issues is setting :code:`mf.rebuild_nsteps=1`.
This turns off incremental building of the K-matrix, which is faster but can accumulate errors between SCF steps.

Direct evaluation of the J matrix can be used by setting :code:`mf.with_df.direct_j = True`, but this is much slower than SGX-J or DF-J and only intended for testing purposes.

Adjustable Parameters
=====================

Calling the :code:`sgx_fit` function on an :mod:`scf.hf.SCF` object returns an equivalent :code:`SGXHF` object with a :code:`with_df` attribute that handles SGX integration. To use a non-default auxiliary basis (for :code:`dfj=True`), :code:`auxbasis` can be specified in the call to :code:`sgx_fit`. The other input option to :code:`sgx_fit` is :code:`pjs`. If :code:`pjs=False` (default), :code:`dfj=optk=False`, otherwise :code:`dfj=optk=True`. In addition, there are various adjustable parameters for SGX that can be set after initialization.

:code:`with_df` attributes:

* :code:`grids_level_i`: The grid level to use for initial SCF iterations.
* :code:`grids_level_f`: The grid level to use for final SCF iterations.
* :code:`dfj`: If :code:`True`, density fitting is used for the J-matrix. If :code:`False`, SGX is used.
* :code:`optk`: If :code:`True`, P-junction screening is used. Threshold determined by :code:`sgx_tol_energy` and :code:`sgx_tol_potential`. Only applies if :code:`dfj=True`.
* :code:`sgx_tol_energy` and :code:`sgx_tol_potential`: Sets the approximate upper bound on the energy and K-matrix errors due to density matrix screening. Ignored if :code:`optk=False` or :code:`dfj=False`. See the Usage and Example section above for instructions on how to set these, as well as the docstrings for more details.
* :code:`use_opt_grids`: Whether to use optimized grids for SGX based on COSX in ORCA :cite:`Helmich2021`. Default is True.
* :code:`fit_ovlp`: Whether to numerically fit the overlap matrix to improve numerical precision. Default is True.
* :code:`grids_thrd`: When the values of all atomic orbitals, multiplied by the qudrature weight, are less than this threshold for a given grid point, it is removed from the SGX numerical grid.
* :code:`grids_switch_thrd`: The threshold for the magnitude of the change in the density matrix that triggers the switch from the initial grid specified by :code:`grids_level_i` to the final grid specified by :code:`grids_level_f`.
* :code:`bound_aglo`: Determines how to estimate the maximum value of the three-center SGX integrals. Options are :code:`"ovlp"`, :code:`"sample"`, and :code:`"sample_pos"` (default). :code:`"ovlp"` assumes the maximum value of each integral is roughly equal to the maximum overlap between orbitals in a shell (with a value of 1 for shells on the same atom to account for orthogonality). :code:`"sample"` constructs a coarse grid for each atom pair and computes the maximum value of the integrals over these grids. :code:`"sample_pos"` does the same but also contructs an approximate estimate for the integral inside each grid block (position-dependent estimate), which can produce a small speedup for global exchange and a larger speedup for short-range exchange. Usually the default value :code:`"sample_pos"` should be fine.
* :code:`direct_j`: If :code:`True`, direct evaluation of the J-matrix is used (slow, for testing only).


:code:`SGXHF` attributes:

* :code:`auxbasis`: The auxiliary basis for J-matrix density fitting. Used if :code:`dfj=True`.
* :code:`rebuild_nsteps`: By default, when :code:`mf.direct_scf=True`, the K-matrix is constructed incrementally from the change in the density matrix between SCF steps (as is the case for other exchange algorithms in PySCF). The SGX J and K matrices are then rebuilt from scratch every :code:`rebuild_nsteps` steps (default 5). Set to 1 to rebuild from scratch at every step, and set to greater than the max number of cycles to turn off rebuilding the JK matrix (Warning: This can cause creeping numerical error). The SGX matrix is built from scratch when the numerical grid changes, regardless of the value of :code:`rebuild_nsteps`, because the same density matrix will have different energies for different grid sizes.
