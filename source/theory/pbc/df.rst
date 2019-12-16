.. _theory_pbc_df:

Density fitting
***************

Introduction
============

Density fitting (DF) techniques are useful to reduce the computational cost of
two-electron repulsion integrals (ERIs), especially when periodic boundary conditions (PBCs) are applied.
In PySCF, there are four types of DF methods available for PBC calculations:
FFTDF (plane-wave density fitting with fast Fourier transformation, which is the default scheme),
AFTDF (plane-wave density fitting with analytical Fourier transformation),
GDF (Gaussian density fitting) and
MDF (mixed Gaussian and plane-wave density fitting).
They are implemented in the :mod:`pbc.df` module.
The characters of these DF methods are summarized in the following table.

========================= =========== =========== ========== ==============
Subject                   FFTDF       AFTDF       GDF        MDF
------------------------- ----------- ----------- ---------- --------------
Initialization            No          No          Slow       Slow
HF Coulomb matrix (J)     Fast        Slow        Fast       Moderate
HF exchange matrix (K)    Slow        Slow        Fast       Moderate
Building ERIs             Slow        Slow        Fast       Moderate
All-electron calculation  Huge error  Large error Accurate   Most accurate
Low-dimension system      N/A         0D,1D,2D    0D,1D,2D   0D,1D,2D
========================= =========== =========== ========== ==============


.. _fftdf:

FFTDF --- Fast Fourier transform based density fitting
======================================================

FFTDF represents the method to compute ERIs in
reciprocal space with the Fourier transformed Coulomb kernel by using
numerical fast Fourier transform (FFT), which is implmented in the
PySCF class :class:`FFTDF`.
An :class:`FFTDF` object can be initialized as follows::

    >>> import numpy as np
    >>> from pyscf.pbc import gto, df, scf
    >>> cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    >>> fftdf = df.FFTDF(cell)
    >>> print(fftdf)
    <pyscf.pbc.df.fft.FFTDF object at 0x111e278d0>
    >>> mf = scf.RHF(cell)
    >>> print(mf.with_df)
    <pyscf.pbc.df.fft.FFTDF object at 0x1206b0f28>

As the default DF scheme of PBC calculations,
an :class:`FFTDF` object is created when initializing the PBC mean-field object and held in the attribute :attr:`with_df`.

4-index ERI tensor and integral transformation
----------------------------------------------
For a general 4-index ERI, we have

.. math::

   (i_{\mathbf{k}_i} j_{\mathbf{k}_j}|k_{\mathbf{k}_k} l_{\mathbf{k}_l}) =
   \Omega \sum_{\mathbf{G}\neq \mathbf{k}_i-\mathbf{k}_j} \rho_{i_{\mathbf{k}_i},j_{\mathbf{k}_j}}(\mathbf{G})
   \frac{4\pi}{| \mathbf{k}_j-\mathbf{k}_i+\mathbf{G}|^2}
   \rho_{k_{\mathbf{k}_k},l_{\mathbf{k}_l}}(\mathbf{G}_{ikjl}-\mathbf{G}) \;,

where :math:`\mathbf{G}` is a reciprocal lattice vector,

.. math::

   \mathbf{G}_{ikjl} = \mathbf{k}_i + \mathbf{k}_k - \mathbf{k}_j - \mathbf{k}_l \;,

and :math:`\rho_{i_{\mathbf{k}_i},j_{\mathbf{k}_j}}(\mathbf{G})` is the Fourier transformed orbital pair density:

.. math::

   \rho_{i_{\mathbf{k}_i},j_{\mathbf{k}_j}}(\mathbf{G})
   = \frac{1}{\Omega} \int_{\Omega} d\mathbf{r} \phi_{i_{\mathbf{k}_i}}^{*}(\mathbf{r}) \phi_{j_{\mathbf{k}_j}}(\mathbf{r})
   e^{-i(\mathbf{k}_j - \mathbf{k}_i + \mathbf{G})\cdot\mathbf{r}} \;.

Note that the four k points (corresponding to the four AO indices) should follow the momentum
conservation law:

.. math::
    (\mathbf{k}_j - \mathbf{k}_i + \mathbf{k}_l - \mathbf{k}_k) = \mathbf{G}.

To evaluate these 4-index ERIs, :class:`FFTDF` provides the function :func:`FFTDF.get_eri`.
By default, four :math:`\Gamma` points are assigned to the four AO indices.
As the format of molecular ERI tensor, the PBC ERI tensor is reshaped to a 2D
array::

    >>> kpts = cell.make_kpts([2,2,2])
    >>> eri = fftdf.get_eri() # \Gamma points only
    >>> print(eri.shape)
    (3, 3)
    >>> eri = fftdf.get_eri([kpts[0],kpts[0],kpts[1],kpts[1]])
    >>> print(eri.shape)
    (4, 4)

In addition, one can perform AO to MO transformations using the function :func:`FFTDF.ao2mo`.
Similar to :func:`FFTDF.get_eri`, the
returned integral tensor is reshaped to a 2D array::

    >>> orbs = np.random.random((4,2,2)) # MO coefficients
    >>> eri_mo = fftdf.ao2mo(orbs, [kpts[0],kpts[0],kpts[1],kpts[1]])
    >>> print(eri_mo.shape)
    (4, 4)

Coulomb and exchange integrals
------------------------------
The :class:`FFTDF` class provides a method :func:`FFTDF.get_jk` to compute
Hartree-Fock Coulomb matrix (J) and exchange matrix (K).  This method can take
one density matrix or a list of density matrices as input and return the J and K
matrices for each density matrix::

    >>> dm = np.random.random((2,2))
    >>> j, k = fftdf.get_jk(dm)
    >>> print(j.shape)
    (2, 2)
    >>> dm = np.random.random((3,2,2))
    >>> j, k = fftdf.get_jk(dm)
    >>> print(j.shape)
    (3, 2, 2)

When k points are specified, the input density matrices should have the correct
shape that matches the number of k points::

    >>> kpts = cell.make_kpts([1,1,3])
    >>> dm = np.random.random((3,2,2))
    >>> j, k = fftdf.get_jk(dm, kpts=kpts)
    >>> print(j.shape)
    (3, 2, 2)
    >>> dm = np.random.random((5,3,2,2))
    >>> j, k = fftdf.get_jk(dm, kpts=kpts)
    >>> print(j.shape)
    (5, 3, 2, 2)


Nuclear type integrals
----------------------

PBC nuclear-electron interaction and pseudo-potential (PP) integrals can be
computed with the FFTDF methods :func:`FFTDF.get_nuc` and :func:`FFTDF.get_pp`.
:func:`FFTDF.get_nuc` function only evaluates the integral of the point charge.
If PP was specified in the cell object, :func:`FFTDF.get_nuc` produces the
integrals of the point nuclei with the effective charges.  If PP was not
defined in the cell object, :func:`FFTDF.get_pp` and :func:`FFTDF.get_nuc`
produce the same integrals.  Depending on the input k-point(s),
the two functions can produce the nuclear-type integrals for a single k-point or
a list of nuclear-type integrals for the k-points.  By default, they compute the
nuclear-type integrals of Gamma point::

    >>> vnuc = fftdf.get_pp()
    >>> print(vnuc.shape)
    (2, 2)
    >>> kpts = cell.make_kpts([2,2,2])
    >>> vnuc = fftdf.get_pp(kpts)
    >>> print(vnuc.shape)
    (8, 2, 2)
    >>> vnuc = fftdf.get_pp()
    >>> print(vnuc.shape)
    (2, 2)


Kinetic energy cutoff
---------------------
The accuracy of FFTDF integrals are affected by the kinetic energy cutoff.  The
default kinetic energy cutoff is a conservative estimation based on the basis
set and the lattice parameter.  You can adjust the attribute :attr:`FFTDF.mesh`
(the numbers of grid points in each positive direction) to change the kinetic
energy cutoff.  If any values in :attr:`FFTDF.mesh` is too small to reach the
required accuracy :attr:`cell.precision`, :class:`FFTDF` may output a warning
message, e.g.::

  WARN: ke_cutoff/mesh (12.437 / [7, 9, 9]) is not enough for FFTDF to get integral accuracy 1e-08.
  Coulomb integral error is ~ 2.6 Eh.
  Recomended ke_cutoff/mesh are 538.542 / [40 40 40].

In this warning message, ``Coulomb integral error`` is a rough estimation for
the largest error of the matrix elements of the two-electron Coulomb integrals.
The overall computational error may be varied by 1 - 2 orders of magnitude.


AFTDF --- Analytic Fourier transform based density fitting
==========================================================

The AFTDF method implements the Fourier transform of the orbital pair density
analytically instead of numerically in the FFTDF case.

To enable AFTDF in the calculation, :class:`AFTDF` object can be initialized
and assigned to :attr:`with_df` object of mean-field object::

    >>> import numpy as np
    >>> from pyscf.pbc import gto, df, scf
    >>> cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    >>> aft = df.AFTDF(cell)
    >>> print(aft)
    <pyscf.pbc.df.aft.AFTDF object at 0x7ff8b1893d90>
    >>> mf = scf.RHF(cell)
    >>> mf.with_df = aft

Generally, AFTDF is slower than FFTDF method.

:class:`AFTDF` class offers the same methods as the :class:`FFTDF` class.
Nuclear and PP integrals, Hartree-Fock J and K matrices, electron repulsion
integrals and integral transformation can be computed with functions
:func:`AFTDF.get_nuc`, :func:`AFTDF.get_pp`, :func:`AFTDF.get_jk`,
:func:`AFTDF.get_eri` and :func:`AFTDF.ao2mo` using the same calling APIs as the
analogy functions in :ref:`fftdf`.


Kinetic energy cutoff
---------------------

:class:`AFTDF` also makes estimation on the kinetic energy cutoff.  When the
any values of :attr:`AFTDF.mesh` are too small for required accuracy
:attr:`cell.precision`, this class also outputs the
``Coulomb integral error`` warning message as the :class:`FFTDF` class.


.. _pbc_gdf:

GDF --- Gaussian density fitting
================================

GDF is an analogy of the conventional density fitting method with periodic
boundary condition.  The auxiliary fitting basis in PBC GDF is periodic Gaussian
function (To ensure the long range Coulomb integrals converging in the real
space lattice summation, the multipoles are removed from the auxiliary basis).
:class:`GDF` object can be initialized and enabled in the SCF calculation in two
ways::

    >>> import numpy as np
    >>> from pyscf.pbc import gto, df, scf
    >>> cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    >>> gdf = df.GDF(cell)
    >>> mf = scf.RHF(cell)
    >>> mf.with_df = gdf
    >>> mf.run()
    >>> # Using SCF.density_fit method
    >>> mf = scf.RHF(cell).density_fit().run()
    >>> print(mf.with_df)
    <pyscf.pbc.df.df.GDF object at 0x7fec7722aa10>

Similar to the molecular code, :func:`SCF.density_fit` method returns a
mean-field object with :class:`GDF` as the integral engine.

In the :class:`GDF` method, the DF-integral tensor is precomputed and stored
on disk.  :class:`GDF` method supports both the :math:`\Gamma`-point ERIs and
the ERIs of different k-points.  :attr:`GDF.kpts` should be specified before
initializing :class:`GDF` object.  :class:`GDF` class provides the same APIs as
the :class:`FFTDF` class to compute nuclear integrals and electron Coulomb
repulsion integrals::

    >>> import numpy as np
    >>> from pyscf.pbc import gto, df, scf
    >>> cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    >>> gdf = df.GDF(cell)
    >>> gdf.kpts = cell.make_kpts([2,2,2])
    >>> gdf.get_eri([kpts[0],kpts[0],kpts[1],kpts[1]])

In the mean-field calculation, assigning :attr:`kpts` attribute to mean-field
object updates the :attr:`kpts` attribute of the underlying DF method::

    >>> import numpy as np
    >>> from pyscf.pbc import gto, df, scf
    >>> cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    >>> mf = scf.KRHF(cell).density_fit()
    >>> kpts = cell.make_kpts([2,2,2])
    >>> mf.kpts = kpts
    >>> mf.with_df.get_eri([kpts[0],kpts[0],kpts[1],kpts[1]])

Once the GDF integral tensor was initialized, the :class:`GDF` can be only used
with certain k-points calculations.  An incorrect :attr:`kpts` argument can lead
to a runtime error::

    >>> import numpy as np
    >>> from pyscf.pbc import gto, df, scf
    >>> cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    >>> gdf = df.GDF(cell, kpts=cell.make_kpts([2,2,2]))
    >>> kpt = np.random.random(3)
    >>> gdf.get_eri([kpt,kpt,kpt,kpt])
    RuntimeError: j3c for kpts [[ 0.53135523  0.06389596  0.19441766]
     [ 0.53135523  0.06389596  0.19441766]] is not initialized.
    You need to update the attribute .kpts then call .build() to initialize j3c.

The GDF initialization is very expensive.  To reduce the initialization cost in
a series of calculations, it would be useful to cache the GDF integral tensor in
a file then load them into the calculation when needed.  The GDF integral tensor
can be saved and loaded the same way as we did for the molecular DF method (see
:ref:`sl_cderi`)::

    import numpy as np
    from pyscf.pbc import gto, df, scf
    cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    gdf = df.GDF(cell, kpts=cell.make_kpts([2,2,2]))
    gdf._cderi_to_save = 'df_ints.h5'  # To save the GDF integrals
    gdf.build()

    mf = scf.KRHF(cell, kpts=cell.make_kpts([2,2,2])).density_fit()
    mf.with_df._cderi = 'df_ints.h5'   # To load the GDF integrals
    mf.run()


Auxiliary Gaussian basis
------------------------

GDF method requires a set of Gaussian functions as the density fitting auxiliary basis.
See also :ref:`df_auxbasis` and :ref:`df_etb_auxbasis` for the choices of DF auxiliary
basis in PySCF GDF code.  There are not many optimized auxiliary basis sets available
for PBC AO basis.  You can use the even-tempered Gaussian functions as the
auxiliary basis in the PBC GDF method::

    import numpy as np
    from pyscf.pbc import gto, df, scf
    cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    gdf = df.GDF(cell, kpts=cell.make_kpts([2,2,2]))
    gdf.auxbasis = df.aug_etb(cell, beta=2.0)
    gdf.build()


Kinetic energy cutoff
---------------------

GDF method does not require the specification of kinetic energy cutoff.
:attr:`cell.ke_cutoff` and :attr:`cell.mesh` are ignored in the :class:`GDF`
class.  Internally, a small set of planewaves is used in the GDF method to
converge the long-range interactions of GDF integrals in the real space lattice
summation. An estimation of energy cutoff is made for the planewaves.
The estimated energy cutoff is converted to the planewave mesh and assigned to
the attribute :class:`GDF.mesh` in the :class:`GDF` class.  It is not
recommended to change this parameter.

The energy cutoff estimation is briefly documented below. In the GDF method, we
introduced a smooth Gaussian function :math:`g(\eta) = \frac{1}{N} r^l e^{-\eta r^2}`
to compensate the long range Coulomb potential of the auxiliary Gaussian basis.
The Coulomb interaction between the smooth auxiliary Gaussian and the rest other
auxiliary Gaussian basis or two-center Gaussian product is calculated in the
reciprocal space.

.. math::
    \sum_{\mathbf{G}} w_\mathbf{G} \frac{4\pi}{G^2}
    \mathrm{FT}(g(\eta)) \rho(\mathbf{G})

FT means Fourier transform.  Considering the leading term of the Fourier
transform :math:`\mathrm{FT}(g(\eta))`

.. math::
    \int e^{-i\mathbf{G}\cdot\mathbf{r}} \frac{1}{N} r^n e^{-\eta r^2} d\mathbf{r}
    = G^n e^{-\frac{k^2}{4\eta}} + \cdots

the Coulomb integral can be estimated

.. math::
    &w \sum_{\mathbf{G}} \frac{4\pi}{G^2} G^n e^{-\frac{G^2}{4\eta}}
    \rho(\mathbf{G})
    \approx \int_{-\infty}^\infty \frac{4\pi}{G^2}
    G^n e^{-\frac{G^2}{4\eta}}\rho(\mathbf{G}) d\mathbf{G} \\
    &\approx (4\pi)^2 \int_0^\infty G^n e^{-\frac{G^2}{4\eta}}\rho(G) dG \\
    &= (4\pi)^2 \int_0^{G_{max}} G^n e^{-\frac{G^2}{4\eta}}\rho(G) dG
    + \varepsilon(G_{max})

:math:`\varepsilon(G_{max})` is the error due to the energy cutoff
:math:`G_{max}`.  The largest error in this integral is the interaction
between :math:`g(\eta)` and a compact density distribution. For the regular
auxiliary DF basis or atomic orbital basis, the most compact function is s type
Gaussian function near nuclear core region. For the very compact function which
is close to the point charge distribution, the Fourier transform form is
approximately a constant :math:`\rho(\mathbf{G}) \sim 1`.

.. math::
    \varepsilon(G_{max})
    &=(4\pi)^2 \int_{G_{max}}^\infty G^{n} e^{-\frac{G^2}{4\eta}} dG
    \\
    &=(4\pi)^2\Big(2\eta G_{max}^{n-1} e^{-\frac{G_{max}^2}{4\eta}}
    + 2\eta(n-1) \int_{G_{max}}^{\infty} G^{n-2} e^{-\frac{G^2}{4\eta}} dG
    + \cdots\Big)

Assuming :math:`G_{max}^2 \gg 2\eta`, we can use the leading term to estimate
the error

.. math::
  \varepsilon(G_{max})
  < 32\pi^2 \eta G_{max}^{n-1} e^{-\frac{G_{max}^2}{4\eta}}

For certain precision requirement :math:`\epsilon`, the energy cutoff can be
evaluated

.. math::
  E = \frac{1}{2}G_{max}^2
  \geq 2\eta \Big((l_{max}-1)\log(G_{max}) - \log(\frac{\epsilon}{32\pi^2 \eta}) \Big)


.. _pbc_mdf:

MDF --- Mixed Gaussian and plane-wave density fitting
=====================================================

The MDF method combines the AFTDF and GDF methods in the same framework.
The MDF auxiliary basis is Gaussian and plane-wave mixed basis.
:class:`MDF` object can be created in two ways::

    >>> import numpy as np
    >>> from pyscf.pbc import gto, df, scf
    >>> cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g', ke_cutoff=10)
    >>> mdf = df.MDF(cell)
    >>> print(mdf)
    <pyscf.pbc.df.mdf.MDF object at 0x7f4025120a10>
    >>> mf = scf.RHF(cell).mix_density_fit().run()
    >>> print(mf.with_df)
    <pyscf.pbc.df.mdf.MDF object at 0x7f7963390a10>

The kinetic energy cutoff is specified in this example to constrain the number of
planewaves.  The number of planewaves can also be controlled through
attribute :attr:`MDF.mesh`.

In principle, the accuracy of MDF method can be increased by adding
more plane waves in the auxiliary basis.  In practice, the linear dependency
between plane waves and Gaussians may lead to numerical stability issue.
The optimal accuracy (with reasonable computational cost) requires a reasonable
size of plan wave basis with a reasonable linear dependency threshold.  A
threshold too large would remove many auxiliary functions while a threshold too
small would cause numerical instability.
.. In our preliminary test, ``ke_cutoff=10`` is able to produce 0.1 mEh accuracy in
.. total energy.
The default linear dependency threshold is 1e-10.  The threshold can be adjusted
through the attribute :attr:`MDF.linear_dep_threshold`.

Like the GDF method, it is also very demanding to initialize the 3-center
Gaussian integrals in the MDF method.  The 3-center Gaussian integral tensor can
be cached in a file and loaded to :class:`MDF` object at the runtime::

    import numpy as np
    from pyscf.pbc import gto, df, scf
    cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g')
    mdf = df.MDF(cell, kpts=cell.make_kpts([2,2,2]))
    mdf._cderi_to_save = 'df_ints.h5'  # To save the GDF integrals
    mdf.build()

    mf = scf.KRHF(cell, kpts=cell.make_kpts([2,2,2])).mix_density_fit()
    mf.with_df._cderi = 'df_ints.h5'   # To load the GDF integrals
    mf.run()


All-electron calculation
========================

All-electron calculations with FFTDF or AFTDF methods requires high energy cutoff
for most elements.  It is recommended to use GDF or MDF methods in the
all-electron calculations.  In fact, GDF and MDF can also be used in PP
calculations to reduce the number of planewave basis if steep functions are
existed in the AO basis.


Low-dimension system
====================

.. In 1.4 release, FFTDF module does not support low-dimension pbc system.

:class:`AFTDF` supports the systems with 0D (molecule), 1D and 2D periodic
boundary conditions.  When computing the integrals of low-dimension systems, an
infinite vacuum is placed on the free boundary.  You can set the
:attr:`cell.dimension`, to enable the integral algorithms for
low-dimension systems in :class:`AFTDF` class::

    import numpy as np
    from pyscf.pbc import gto, df, scf
    cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g', dimension=1)
    aft = df.AFTDF(cell)
    aft.get_eri()

:class:`GDF` and :class:`MDF` all support the integrals of low-dimension system.
Similar to the usage of AFTDF method, you need to set :attr:`cell.dimension` for
the low-dimension systems::

    import numpy as np
    from pyscf.pbc import gto, df, scf
    cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g', dimension=1)
    gdf = df.GDF(cell)
    gdf.get_eri()

See more examples in ``examples/pbc/31-low_dimensional_pbc.py``


Interface to molecular DF-post-HF methods
=========================================

PBC DF object is compatible to the molecular DF object.  The
:math:`\Gamma`-point PBC SCF object can be directly passed to molecular DF
post-HF methods for an electron correlation calculations in PBC::

    import numpy as np
    from pyscf.pbc import gto, df, scf
    from pyscf import cc as mol_cc
    cell = gto.M(atom='He 1 1 1', a=np.eye(3)*2, basis='3-21g', dimension=1)
    mf = scf.RHF(cell).density_fit()
    mol_cc.RCCSD(mf).run()


