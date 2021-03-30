.. _user_dft:

*******************************
Density functional theory (DFT)
*******************************

*Modules*: :mod:`dft`, :mod:`pbc.dft`

.. _user_dft_intro:

Introduction
============

Kohn-Sham density functional theory (KS-DFT) has been implemented through derived classes of the :class:`pyscf.scf.hf.SCF` parent class. As such, the methods and capabilities introduced in :numref:`user_scf` are also available to the :mod:`dft` module.

A minimal example of using the :mod:`dft` module reads:

  >>> from pyscf import gto, dft
  >>> mol = gto.M(atom = 'H 0 0 0; F 0 0 1.1', basis = 'ccpvdz', symmetry = True)
  >>> mf = dft.RKS(mol)
  >>> mf.xc = 'lda,vwn' # default
  >>> mf.kernel()

This will run a restricted, closed-shell Kohn-Sham DFT calculation with the default LDA functional.

.. _user_dft_theory:

Theory
======

In KS-DFT, as first proposed by Kohn and Sham :cite:`KohSha1965`, the electron density of a reference noninteracting system is used to represent the density of the true interacting system. As a result, the computational formulation of KS-DFT resembles that of Hartree-Fock (HF) theory, but with a different effective Fock potential. In KS-DFT, the total electronic energy is defined as follows:

.. math::

    E = T_s + E_{\rm ext} + E_J + E_{\rm xc} \ .

Here, :math:`T_s` is the noninteracting kinetic energy, :math:`E_{\rm ext}` is the energy due to the external potential, :math:`E_J` is the Coulomb energy, and
:math:`E_{\rm xc}` is the exchange-correlation (*xc*) energy. In practice, :math:`E_{\rm xc}` is approximated by a density functional approximation, which themselves may be divided into several classes along different rungs of Jacob's ladder :cite:`perdew_jacobs_ladder_aip_conf_proc_2001`:

- local density approximations (e.g. LDA; *xc* energy depends only on the electron density, :math:`\rho`), 
- generalized gradient approximations (GGA; *xc* energy also depends on the density gradient, :math:`|\nabla\rho|`), 
- meta-GGAs (*xc* energy also depends on the kinetic energy density and Laplacian, :math:`\sum_i |\nabla \psi_i|^2`, :math:`\nabla^2\rho`),
- non-local correlation functionals (*xc* energy involves a double integral)
- hybrid density functionals (a fraction of exact exchange is used), and
- long-range corrected density functionals (exact exchange is used with a modified interaction kernel)

Variationally minimizing the total energy with respect to the density yields the KS equations for obtaining the non-interacting reference orbitals, on par with HF theory, and these have the same general form as the Fock equations in :numref:`theory_scf`. However, the exact exchange, :math:`\hat{K}`, is replaced by the *xc* potential, :math:`\hat{v}_{\rm xc}=\delta E_{\rm xc}/\delta \rho`. For hybrid and meta-GGA calculations, PySCF uses the generalized KS formalism :cite:`GKS`, in which the so-called generalized KS equations minimize the total energy with respect to the orbitals themselves. 

.. _user_dft_predef_func:

Predefined *xc* Functionals and Functional Aliases
==================================================

The choice of *xc* functional is assigned via the attribute :attr:`DFT.xc`. This is a comma separated string (precise grammar discussed :ref:`below <user_dft_custom_func>`, e.g., ``xc = 'pbe,pbe'`` denotes PBE exchange plus PBE correlation. In common usage, a single name (alias) is often used to refer to the combination of a particular exchange and correlation approximation instead. To support this, PySCF will first examine a lookup table to see if :attr:`DFT.xc` corresponds to a common compound name, and if so, the implementation dispatches to the appropriate exchange and correlation forms, e.g., ``xc = 'pbe'`` directly translates to ``xc = 'pbe,pbe'``. However, if the name is not found in the compound functional table, and only a single string is given, it will be treated as an exchange functional only, e.g., ``xc = 'b86'`` leads to B86 exchange only (without correlation). Please note that earlier PySCF versions (1.5.0 or earlier)
did not support compound functional aliases, and both exchange and correlation always had to be explicitly assigned. 

PySCF supports two independent libraries of *xc* functional implementations, namely `Libxc <https://www.tddft.org/programs/libxc/>`_ and `XCFun <https://xcfun.readthedocs.io/en/latest/>`_. The former of these is the default, but the latter may be chosen upon by setting ``mf._numint.libxc = dft.xcfun``, cf. `dft/32-xcfun_as_default.py <https://github.com/pyscf/pyscf/blob/master/examples/dft/32-xcfun_as_default.py>`_. For complete lists of available functionals, the user is referred to `pyscf/dft/libxc.py <https://github.com/pyscf/pyscf/blob/master/pyscf/dft/libxc.py>`_ and `pyscf/dft/xcfun.py <https://github.com/pyscf/pyscf/blob/master/pyscf/dft/xcfun.py>`_, respectively.

.. _user_dft_custom_func:

Customizing XC functionals
==========================

XC functionals of DFT methods can be customized. The simplest way to customize
the XC functional is to assign a string expression to :attr:`mf.xc`::

    from pyscf import gto, dft
    mol = gto.M(atom='H  0  0  0; F  0.9  0  0', basis='6-31g')
    mf = dft.RKS(mol)
    mf.xc = 'HF*0.2 + .08*LDA + .72*B88, .81*LYP + .19*VWN'
    mf.kernel()
    mf.xc = 'HF*0.5 + .08*LDA + .42*B88, .81*LYP + .19*VWN'
    mf.kernel()
    mf.xc = 'HF*0.8 + .08*LDA + .12*B88, .81*LYP + .19*VWN'
    mf.kernel()
    mf.xc = 'HF'
    mf.kernel()

The XC functional string is parsed against the rules, as described below.

* The given functional description must be a one-line string.

* The functional description is case-insensitive.

* The functional description string has two parts, separated by ``,``.  The
  first part describes the exchange functional, the second part sets the
  correlation functional.

  - If "," does not appear in the string, the entire string is treated as the name of a
    compound functional (containing both the exchange and the correlation
    functional) which should be in the functional aliases list. See
    the list of predefined XC functionals in the section above.

    If the string is not found in the aliased functional list, it is treated as
    an X functional.

  - To input only an X functional (without a C functional), leave the second part
    blank. E.g. description='slater,' means a functional with the LDA contribution
    only.

  - To neglect the contribution of the X functional (i.e. to just use a C functional), leave
    the first part blank, e.g. description=',vwn' means a functional with VWN
    only.

  - If compound XC functional is specified, no matter whether it is in the X
    part (the string in front of the comma) or the C part (the string behind the comma),
    both X and C functionals of the compound XC functional will be used.

* The functional name can be placed in an arbitrary order.  Two names need to be
  separated by operators "+" or "-".  Blank spaces are ignored.  NOTE the parser
  only reads operators "+" "-" "*".  / is not supported.

* A functional name can have at most one factor.  If the factor is not given, it
  is set to 1.  Compound functionals can be scaled as a unit. For example
  '0.5*b3lyp' is equivalent to 'HF*0.1 + .04*LDA + .36*B88, .405*LYP + .095*VWN'

* String "HF" stands for exact exchange (HF K matrix).  "HF" can be put in the
  correlation functional part (after the comma). Putting "HF" in the correlation
  part is the same as putting "HF" in the exchange part.

* String "RSH" means range-separated operator. Its format is RSH(alpha; beta;
  omega).  Another way to input RSH is to use keywords SR_HF and LR_HF:
  "SR_HF(0.1) * alpha_plus_beta" and "LR_HF(0.1) * alpha" where the number in
  parenthesis is the value of omega.

* Be careful with the libxc convention of GGA functionals, in which the LDA
  contribution is included.


There is also another way to customize XC functionals which uses the :py:meth:`eval_xc`
method of the numerical integral class::

    mol = gto.M(atom='H 0 0 0; F 0.9 0 0', basis = '6-31g')
    mf = dft.RKS(mol)
    def eval_xc(xc_code, rho, spin=0, relativity=0, deriv=1, verbose=None):
        # A fictitious XC functional to demonstrate the usage
        rho0, dx, dy, dz = rho
        gamma = (dx**2 + dy**2 + dz**2)
        exc = .01 * rho0**2 + .02 * (gamma+.001)**.5
        vrho = .01 * 2 * rho0
        vgamma = .02 * .5 * (gamma+.001)**(-.5)
        vlapl = None
        vtau = None
        vxc = (vrho, vgamma, vlapl, vtau)
        fxc = None  # 2nd order functional derivative
        kxc = None  # 3rd order functional derivative
        return exc, vxc, fxc, kxc
    dft.libxc.define_xc_(mf._numint, eval_xc, xctype='GGA')
    mf.kernel()

By calling the :func:`dft.libxc.define_xc_` function, the customized :func:`eval_xc`
function is patched to the numerical integration class :attr:`mf._numint`
dynamically.

More examples of DFT XC functional customization can be found in
:source:`examples/dft/24-custom_xc_functional.py` and
:source:`examples/dft/24-define_xc_functional.py`.

Numerical integration grids
===========================
PySCF implements several numerical integration grids,
which can be tuned in DFT calculations following the examples in 
:source:`examples/dft/11-grid_scheme.py`.
In addition, these grids can be used for the general numerical evaluation of
basis functions, electron densities, and integrals.
Some examples can be found in 
:source:`examples/dft/30-ao_value_on_grid.py`, and
:source:`examples/dft/31-xc_value_on_grid.py`.
The following is an example that computes the kinetic energy from the 
nonnegative kinetic energy density

.. math::

    t_s(\mathbf{r}) = \frac{1}{2} \sum_{i\in occ} |\nabla\psi_i(\mathbf{r})|^2 \;,

.. math::

    T_s = \int d\mathbf{r} t_s(\mathbf{r}) \;.

.. code-block:: python

    from pyscf.dft import gen_grid, numint
    orbo = mf.mo_coeff[:,mf.mo_occ>0]
    grids = gen_grid.Grids(mol)
    grids.build(with_non0tab=True)
    weights = grids.weights
    ao1 = numint.eval_ao(mol, grids.coords, deriv=1, non0tab=grids.non0tab)
    ts = 0.5 * numpy.einsum('xgp,pi,xgq,qi->g', ao1[1:], orbo, ao1[1:], orbo)
    Ts = numpy.einsum('g,g->', weights, ts)

    Ts_ao = mol.intor("int1e_kin")
    Ts_anal = np.einsum("ui,uv,vi->", orbo, Ts_ao, orbo)
    print(asb(Ts - Ts_anal))

Dispersion corrections
======================
Grimme's "D3" dispersion correction :cite:`DFTD3` can be added with
an interface to the external library `libdftd3 <https://github.com/cuanto/libdftd3>`_.
See :mod:`dftd3`.

References
==========
.. bibliography:: ref_dft.bib
   :style: unsrt
