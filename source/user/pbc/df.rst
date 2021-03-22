.. _user_pbc_df:

Density fitting for crystalline calculations
********************************************

*Modules*: :mod:`df`, :mod:`pbc.df`

Introduction
============

Like for molecular calculations, density fitting (DF) is useful to reduce the computational cost of manipulating the electron repulsion integrals (ERIs) in periodic calculations. Currently, PySCF provides three periodic DF methods that are different in the auxiliary basis functions being used.

* The fast Fourier transform density fitting (:class:`FFTDF`) uses plane waves (PWs) as the auxiliary basis functions. FFTDF is also known as the Gaussian and plane wave (GPW) approach in literature and other software packages. :cite:`VandeVondele05CPC,Kuhne20JCP`

* The Gaussian density fitting (:class:`GDF`) uses Gaussian-type orbitals (GTOs) as the auxiliary basis functions. :cite:`Sun17JCP`

* The mixed density fitting (:class:`MDF`) uses both PWs and GTOs as the auxiliary basis functions. :cite:`Sun17JCP`

This section of the user document covers the basic usage of these DF methods in periodic calculations. The recommended choice of DF methods for different applications is also provided. See :ref:`theory_df` for the use of DF in molecular calculations.


Using DF in periodic calculations
=================================

Unlike for molecular calculations, DF is used by default in periodic SCF calculations. For example, initialize a KRHF object uses FFTDF by default::

    import numpy as np
    from pyscf.pbc import gto, scf
    cell = gto.M(atom="He 0 0 0", a=np.eye(3)*2, basis="6-31g")
    kpts = cell.make_kpts([2,2,2])
    mf = scf.KRHF(cell, kpts)
    print(mf.with_df)  # <pyscf.pbc.df.fft.FFTDF at 0x7feb78373c88>

The periodic SCF class also provides methods for using GDF::

    mf = scf.KRHF(cell, kpts).density_fit()

and MDF::

    mf = scf.KRHF(cell, kpts).mix_density_fit()

In addition to using the APIs above, the user can also initialize a DF method first and overwrite the attribute :attr:`with_df` of the SCF object. For example::

    >>> from pyscf.pbc import df
    >>> mydf = df.GDF(cell, kpts).build()
    >>> mf.with_df = mydf

Once the periodic mean-field calculation is done, the subsequent correlated calculations will automatically use the same DF method to handle ERIs. For example, a MP2 calculation using GDF can be performed as follows::

    >>> ... # initialize cell and kpts
    >>> from pyscf.pbc import scf, mp
    >>> mf = scf.KRHF(mf, kpts).density_fit()
    >>> mf.kernel()         # perform HF with GDF
    >>> mmp = mp.KMP2(mf)
    >>> mmp.kernel()        # perform MP2 with GDF


Controlling DF error
====================

The DF error is introduced by the incompleteness of the finite auxiliary basis used to expand the atomic orbital pair densities. The DF error can often be reduced by increasing the number of auxiliary basis functions being used.

.. _df_err_fftdf:

FFTDF
-----

FFTDF uses plane waves (PWs) as the auxiliary basis, whose size is determined by :attr:`FFTDF.mesh`, which is set to :attr:`Cell.mesh` upon initialization. :attr:`Cell.mesh` is a 1d array-like object of three integer numbers, ``[nx, ny, nz]``, that defines the number of PWs (or the real-space grid points in the unit cell) for :math:`x`, :math:`y` and :math:`z` directions, respectively. The total number of PWs being used for FFTDF is hence ``nx * ny * nz``. PySCF determines :attr:`Cell.mesh` from :attr:`Cell.ke_cutoff` - the kinetic energy cutoff. By default, :attr:`Cell.ke_cutoff` is determined by requiring the PW expansion of the most compact atomic orbitals in the basis set to achieve :attr:`Cell.precision`.

To use a PW basis of size other than the default, the user can either overwrite :attr:`FFTDF.mesh` directly, or change it by changing :attr:`Cell.ke_cutoff`. An example is provided as follows::

    import numpy as np
    from pyscf.pbc import gto, df

    def print_mesh(mesh):
        print("mesh = [%d, %d, %d]  (%d PWs)" % (*mesh, np.prod(mesh)))

    cell = gto.M(atom="He 0 0 0", a=np.eye(3)*2, basis="gth-dzvp", pseudo="gth-pade")
    kpts = cell.make_kpts([2,2,2])
    mydf = df.FFTDF(cell, kpts)
    print_mesh(mydf.mesh)
    # output: mesh = [42, 42, 42]  (74088 PWs)
    mydf.mesh = [17,17,17]
    print_mesh(mydf.mesh)
    # output: mesh = [17, 17, 17]  (4913 PWs)
    cell.ke_cutoff = 60   # unit: Rydberg
    cell.build()          # rebuild cell to update cell.mesh
    mydf = df.FFTDF(cell, kpts)
    print_mesh(mydf.mesh)
    # output: mesh = [14, 14, 14]  (2744 PWs)

Note that PySCF's default for :attr:`Cell.precision` is relatively conservative (:math:`10^{-8} E_h`). This often leads to a :attr:`Cell.ke_cutoff` higher than the default used by other packages using FFTDF. :cite:`Kuhne20JCP`
For a more cost-effective choice, the user can scan over a series of :attr:`mesh` size (or :attr:`ke_cutoff` values) and monitor the convergence of e.g., the converged SCF energies.


GDF
---

GDF uses Gaussian-type orbitals (GTOs) as the auxiliary basis and parallels the :mod:`df` module for molecular calculations. We guide the readers to :ref:`choice_of_auxbasis` for more details on how to choose the auxiliary basis sets for GDF. A GDF-specific example can be found in :source:`examples/pbc/35-gaussian_density_fit.py`.


MDF
---

MDF uses mixed GTOs and PWs as the fitting basis. The GTO part of the auxiliary basis can be set in the same way as for GDF (again, see :ref:`choice_of_auxbasis`), while the PW part is similar to FFTDF, i.e., setting either :attr:`MDF.mesh` or :attr:`Cell.ke_cutoff`. The default size of the PW basis is also relatively conservative, and the user is recommended to conduct a convergence test as mentioned in :ref:`df_err_fftdf`.


Saving and reusing DF integrals
===============================

While FFTDF is implemented in the so-called integral-direct manner and needs no "initialization", both GDF and MDF pre-compute required integrals and dump them to disk for later use. The APIs for saving and reusing GDF/MDF integrals are the same as in the molecular :mod:`df` module; we guide the user to :ref:`save_reuse_df_integrals` for a detailed description. An example is provided in :source:`examples/pbc/35-gaussian_density_fit.py`.


Choice of DF methods
====================

The choice of DF methods depends on the type of applications, the required accuracy, and the available computational resources.

Type of calculations
--------------------

* For **all-electron** calculations, only GDF and MDF can be used, while FFTDF would require an impractically large PW basis to describe the core orbitals accurately (hydrogen and helium are two exceptions since they don't have core orbitals).

* For **pseudopotential (PP)**-based calculations, all three DF methods can be used.

* For calculations on **low-dimension** systems (0D, 1D, and 2D), only GDF and MDF can be used. The user needs to specify the dimension by setting :attr:`Cell.dimension`.

Required accuracy
-----------------

* FFTDF can be regarded "exact" for PP-based calculations within the given AO basis if a sufficiently large PW basis is used.

* GDF has a typical error of :math:`10^{-5} \sim 10^{-4} E_h` when using PySCF's default auxiliary basis. :cite:`Sun17JCP,Ye20arXiv` This error can be reduced by using a larger auxiliary basis set. :cite:`Ye20arXiv`

* MDF is in general more accurate than GDF and comparable to FFTDF if a sufficiently large PW basis is used. The typical error of MDF is :math:`10^{-6} E_h` or lower with the default parameters. :cite:`Sun17JCP`

Computational resources
-----------------------

* GDF requires **enough disk space** to hold the pre-computed DF integrals. The size of these integrals grows quickly with the system size and scales as :math:`O(N_k^2 n_{\mathrm{AO}}^2 n_{\mathrm{aux}})`, where :math:`N_k` is the number of k-points, :math:`n_{\mathrm{AO}}` is the number of AOs per unit cell, and :math:`n_{\mathrm{aux}}` is the number of auxiliary basis functions per unit cell. Note that for DFT calculations using pure exchange-correlation functionals (LDA and GGA), the storage requirement reduces to :math:`O(N_k n_{\mathrm{AO}}^2 n_{\mathrm{aux}})`, which is much more modest.

* FFTDF uses very little disk space but requires :math:`O(N_k n_{\mathrm{AO}}^2 N_{\mathrm{PW}})` **memory**, where :math:`N_{\mathrm{PW}}` is the size of the PW basis. Despite the modest, linear dependence on :math:`N_k`, the memory requirement could be high for systems that require a relatively large PW basis.

* MDF requires both a :math:`O(N_k^2 n_{\mathrm{AO}}^2 n_{\mathrm{aux}})` disk space to store pre-computed integrals of the GTO part of the auxiliary basis and a :math:`O(N_k n_{\mathrm{AO}}^2 N_{\mathrm{PW}})` memory for the PW part. However, both :math:`n_{\mathrm{aux}}` and :math:`N_{\mathrm{PW}}` here are smaller than that required by GDF and FFTDF, respectively.


References
==========

.. bibliography:: ../ref_df.bib
   :style: unsrt
