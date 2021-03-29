.. _user_dft:

*******************************
Density functional theory (DFT)
*******************************

*Modules*: :mod:`dft`, :mod:`pbc.dft`

Introduction
============
Kohn-Sham density functional theory (KS-DFT) is implemented through derived 
classes of the :class:`pyscf.scf.hf.SCF` class, thus, the methods and capabilities introduced in 
:numref:`user_scf` also apply to the :mod:`dft` module.

A minimal example of using the :mod:`dft` module is as follows. ::

    from pyscf import gto, scf
    mol = gto.M(
        atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
        basis = 'ccpvdz',
        symmetry = True,
    )
    mf = scf.RKS(mol)
    mf.xc = 'lda,vwn' #default functional
    mf.kernel()

This will run a restricted closed-shell Kohn-Sham DFT calculation with the default LDA functional.

Theory
======
KS-DFT, first proposed by Kohn and Sham :cite:`KohSha1965`,
uses the electron density of a reference noninteracting system
to represent the density of the true interacting system. 
As a result, the computational formulation of KS-DFT resembles that of Hartree-Fock (HF) theory.
KS-DFT defines the total electronic energy as follows

.. math::

    E = T_s + E_{\rm ext} + E_J + E_{\rm xc} \;,

where :math:`T_s` is the noninteracting kinetic energy,
:math:`E_{\rm ext}` is the energy due to the external potential,
:math:`E_J` is the Coulomb energy, and
:math:`E_{\rm xc}` is the exchange-correlation (XC) energy.
In practice :math:`E_{\rm xc}` is approximated by a density functional approximation. These
may be divided into several classes, such as

- local density approximations (e.g. LDA; XC energy depends only on the electron density :math:`\rho`), 
- generalized gradient approximations (GGA; XC energy also depends on the density gradient :math:`|\nabla\rho|`), 
- meta-GGAs (XC energy also depends on the kinetic energy density and Laplacian :math:`\sum_i |\nabla \psi_i|^2`, :math:`\nabla^2\rho`),
- non-local correlation functionals (XC energy involves a double integral)
- hybrid density functionals (a fraction of exact exchange is used), and
- long-range corrected density functionals (exact exchange is used with a modified interaction kernel)

Variationally minimizing the total energy with respect to the density leads to the KS equations
for the non-interacting reference orbitals. These have the same form as the Fock equations in :numref:`user_scf`,
but with exact exchange :math:`\hat{K}` replaced by the XC potential :math:`\hat{v}_{\rm xc}=\delta E_{\rm xc}/\delta \rho`.
For hybrid and meta-GGa calculations, PySCF uses the generalized KS formalism :cite:`GKS` where the
generalized KS equations minimize the total energy with respect to the orbitals themselves. 


Predefined XC functionals and functional aliases
================================================
The XC functional is assigned via the attribute :attr:`xc`. This
is a comma separated string (precise grammar discussed below [LINK]) e.g.
xc = 'pbe,pbe' denotes PBE exchange plus PBE correlation. 
In common usage, a single name is often used to refer to the combination of a particular
exchange and correlation approximation.
To support this, PySCF will first examine a lookup table
to see if :attr:`xc` corresponds to a common compound name. 
If so,  the implementation  dispatches to the
appropriate exchange and correlation forms, e.g.  xc = 'pbe' translates to 
xc = 'pbe,pbe'. If the name is not found in
the compound functional table and only a single string is found, it will be treated as an exchange
functional only, e.g. xc = 'b86' leads to B86 exchange only (without correlation).
Note that  earlier PySCF versions (1.5.0 or earlier)
did not support compound functional aliases, and both exchange and correlation always had to be
explicitly assigned.

See examples in 
:source:`examples/dft/00-simple_dft.py`,
:source:`examples/dft/32-xcfun_as_default.py`, and
:source:`examples/dft/33-nlc_functionals.py`.

libxc
-----

The hybrid functionals and predefined compound functionals for libxc (the
default XC library) are

=============  ========================================
Functional     Comments
-------------  ----------------------------------------
B3PW91         The original (ACM) hybrid of Becke
B3LYP          The (in)famous B3LYP
B3P86          Perdew 86 hybrid similar to B3PW91
O3LYP          hybrid using the optx functional
PBE0
PBE1PBE
PBEH           aka PBE0 or PBE1PBE
X3LYP          hybrid by Xu and Goddard
MPW3PW         mixture with the mPW functional
B1LYP          Becke 1-parameter mixture of B88 and LYP
B1PW91         Becke 1-parameter mixture of B88 and PW91
MPW1PW         Becke 1-parameter mixture of mPW91 and PW91
MPW3LYP        mixture of mPW and LYP
HSE03          the 2003 version of the screened hybrid HSE
HSE06          the 2006 version of the screened hybrid HSE
CAMB3LYP       CAM version of B3LYP
BHANDH         Becke half-and-half
BHANDHLYP      Becke half-and-half with B88 exchange
REVB3LYP       Revised B3LYP
CAMYBLYP       BLYP with yukawa screening
B3LYPS         B3LYP* functional
WB97           Chai and Head-Gordon
WB97X_V        Mardirossian and Head-Gordon
LC_VV10        Vydrov and Van Voorhis
CAMYB3LYP      B3LYP with Yukawa screening
WB97X_D        Chai and Head-Gordon
LRC_WPBE       Long-range corrected functional by Rorhdanz et al
B3LYP5         B3LYP with VWN functional 5 instead of RPA
LC_WPBE        Long-range corrected functional by Vydrov and Scuseria
HSE12          HSE12 by Moussa, Schultz and Chelikowsky
MPW1LYP        Becke 1-parameter mixture of mPW91 and LYP
MPW1PBE        Becke 1-parameter mixture of mPW91 and PBE
B88B95         Mixture of B88 with BC95 (B1B95)
B86B95         Mixture of B86 with BC95
PW86B95        Mixture of PW86 with BC95
M06            M06 functional from Minnesota
M06_2X         M06-2X functional from Minnesota
PW6B95         Mixture of PW91 with BC95 from Zhao and Truhlar
PWB6K          Mixture of PW91 with BC95 from Zhao and Truhlar for kinetics
TPSSH          TPSS hybrid
REVTPSSH       revTPSS hybrid
M11            M11 functional from Minnesota
WB97M_V        Mardirossian and Head-Gordon
B3LYP          aka B3LYP5
B3LYP5         .2*HF + .08*SLATER + .72*B88, .81*LYP + .19*VWN5
B3LYPG         .2*HF + .08*SLATER + .72*B88, .81*LYP + .19*VWN3
B3P86          aka B3P865
B3P865         .2*HF + .08*SLATER + .72*B88, .81*P86 + .19*VWN5
B3P86G         .2*HF + .08*SLATER + .72*B88, .81*P86 + .19*VWN3
B3PW91         aka B3PW915
B3PW915        .2*HF + .08*SLATER + .72*B88, .81*PW91 + .19*VWN5
B3PW91G        .2*HF + .08*SLATER + .72*B88, .81*PW91 + .19*VWN3
O3LYP          .1161*HF + 0.071006917*SLATER + .8133*OPTX, .81*LYP + .19*VWN5
MPW3PW         aka MPW3PW5
MPW3PW5        .2*HF + .08*SLATER + .72*MPW91, .81*PW91 + .19*VWN5
MPW3PWG        .2*HF + .08*SLATER + .72*MPW91, .81*PW91 + .19*VWN3
MPW3LYP        aka MPW3LYP5
MPW3LYP5       .218*HF + .073*SLATER + .709*MPW91, .871*LYP + .129*VWN5
MPW3LYPG       .218*HF + .073*SLATER + .709*MPW91, .871*LYP + .129*VWN3
REVB3LYP       aka REVB3LYP5
REVB3LYP5      .2*HF + .13*SLATER + .67*B88, .84*LYP + .16*VWN5
REVB3LYPG      .2*HF + .13*SLATER + .67*B88, .84*LYP + .16*VWN3
X3LYP          aka X3LYP5
X3LYP5         .218*HF + .073*SLATER + .542385*B88 + .166615*PW91, .871*LYP + .129*VWN5
X3LYPG         .218*HF + .073*SLATER + .542385*B88 + .166615*PW91, .871*LYP + .129*VWN3
B5050LYP       .5*HF + .08*SLATER + .42*B88, .81*LYP + .19*VWN
MPW1LYP        .25*HF + .75*MPW91, LYP
MPW1PBE        .25*HF + .75*MPW91, PBE
PBE50          .5*HF + .5*PBE, PBE
REVPBE0        .25*HF + .75*PBE_R, PBE
TPSS0          .25*HF + .75*TPSS, TPSS
OPTXCORR       0.7344536875999693*SLATER - 0.6984752285760186*OPTX
=============  ========================================

The list above is incomplete. Please refer to the libxc manual
(https://www.tddft.org/programs/libxc/functionals/) for a complete list of 
hybrid functionals.

The (aliased) pure functionals for libxc are

==================  ==================  ====================
Functional alias    exchange part       correlation part
------------------  ------------------  --------------------
BLYP                B88                 LYP
BP86                B88                 P86
PW91                PW91                PW91
PBE                 PBE                 PBE
REVPBE              PBE_R               PBE
PBESOL              PBE_SOL             PBE_SOL
PKZB                PKZB                PKZB
TPSS                TPSS                TPSS
REVTPSS             REVTPSS             REVTPSS
SCAN                SCAN                SCAN
SOGGA               SOGGA               PBE
BLOC                BLOC                TPSSLOC
OLYP                OPTX                LYP
OPBE                OPTX                PBE
RPBE                RPBE                PBE
BPBE                B88                 PBE
MPW91               MPW91               PW91
HFLYP               HF                  LYP
HFPW92              HF                  PW_MOD
SPW92               SLATER              PW_MOD
SVWN                SLATER              VWN
MS0                 MS0                 REGTPSS
MS1                 MS1                 REGTPSS
MS2                 MS2                 REGTPSS
MS2H                MS2H                REGTPSS
MVS                 MVS                 REGTPSS
MVSH                MVSH                REGTPSS
SOGGA11             SOGGA11             SOGGA11
SOGGA11-X           SOGGA11_X           SOGGA11_X
KT1                 KT1                 VWN
DLDF                DLDF                DLDF
GAM                 GAM                 GAM
M06-L               M06_L               M06_L
M11-L               M11_L               M11_L
MN12-L              MN12_L              MN12_L
MN15-L              MN15_L              MN15_L
N12                 N12                 N12
N12-SX              N12_SX              N12_SX
MN12-SX             MN12_SX             MN12_SX
MN15                MN15                MN15
MBEEF               MBEEF               PBE_SOL
SCAN0               SCAN0               SCAN
PBEOP               PBE                 OP_PBE
BOP                 B88                 OP_B88
REVSCAN             REVSCAN             REVSCAN
REVSCAN_VV10        REVSCAN             REVSCAN_VV10
SCAN_VV10           SCAN                SCAN_VV10
SCAN_RVV10          SCAN                SCAN_RVV10
==================  ==================  ====================

Libxc provides implementation of individual exchange and correlation
functionals, such as B86, P88, LYP, VWN, etc.  Please refer to libxc manual
(https://www.tddft.org/programs/libxc/functionals/) for the supported
functionals.

xcfun
-----

Another XC functional library that PySCF supports is xcfun
(http://dftlibs.org/xcfun/). The Xcfun library can evaluate arbitrary derivatives of
XC functionals. The predefined compound functionals in xcfun are

=============  ========================================
Functional     Comments
-------------  ----------------------------------------
PBE0           .25*HF + .75*PBEX + PBEC
PBE1PBE        aka PBE0
PBEH           aka PBE0
B3P86          .2*HF + .08*SLATER + .72*B88 + .81*P86C + .19*VWN5C
B3P86G         .2*HF + .08*SLATER + .72*B88 + .81*P86C + .19*VWN3C
B3PW91         .2*HF + .08*SLATER + .72*B88 + .81*PW91C + .19*VWN5C
B3PW91G        .2*HF + .08*SLATER + .72*B88 + .81*PW91C + .19*VWN3C
B3LYP          aka B3LYP5
B3LYP5         .2*HF + .08*SLATER + .72*B88 + .81*LYP + .19*VWN5C
B3LYPG         .2*HF + .08*SLATER + .72*B88 + .81*LYP + .19*VWN3C
O3LYP          .1161*HF + 0.071006917*SLATER + .8133*OPTX, .81*LYP + .19*VWN5
X3LYP          .218*HF + .073*SLATER + 0.542385*B88 + .166615*PW91X + .871*LYP + .129*VWN5C
X3LYPG         .218*HF + .073*SLATER + 0.542385*B88 + .166615*PW91X + .871*LYP + .129*VWN3C
CAMB3LYP       0.19*SR_HF(0.33) + 0.65*LR_HF(0.33) + BECKECAMX + VWN5C*0.19 + LYPC*0.81
B97XC          B97X + B97C + HF*0.1943
B97_1XC        B97_1X + B97_1C + HF*0.21
B97_2XC        B97_2X + B97_2C + HF*0.21
M05XC          .28*HF + .72*M05X + M05C
TPSSH          0.1*HF + 0.9*TPSSX + TPSSC
OLYP           2.4832*SLATER - 1.43169*OPTX + LYP
HFLYP          HF + LYP
KT1            Keal-Tozer 1, JCP, 119, 3015
               SLATERX - 0.006*KTX
KT2XC          Keal-Tozer 2, JCP, 119, 3015
               1.07173*SLATER - .006*KTX + 0.576727*VWN5
KT3XC          Keal-Tozer 3, JCP, 121, 5654
               SLATERX*1.092 + KTX*-0.004 + OPTXCORR*-0.925452 + LYPC*0.864409
=============  ========================================

The (aliased) pure functionals are

==================  ==================  ====================
Functional alias    exchange part       correlation part
------------------  ------------------  --------------------
BPW91               B88                 PW91C
BPW92               B88                 PW92C
BLYP                B88                 LYP
BP86                B88                 P86
PW91                PW91                PW91
PBE                 PBE                 PBE
REVPBE              REVPBE              PBE
PBESOL              PBESOL              PBESOL
TPSS                TPSS                TPSS
REVTPSS             REVTPSS             REVTPSS
SCAN                SCAN                SCAN
BLOC                BLOC                TPSSLOC
OLYP                OPTX                LYP
RPBE                RPBE                PBE
BPBE                B88                 PBE
SVWN                SLATER              VWN5
KT1                 KT1                 VWN
M06-L               M06L                M06L
==================  ==================  ====================

Individual exchange functionals (and kinetic functionals) in xcfun are

=============  ========================================
Functional     Comments
-------------  ----------------------------------------
SLATER         Slater LDA exchange
LDA            aka SLATER
PW86           PW86 exchange
PBE            PBE Exchange Functional
BECKE          Becke 88 exchange
BECKECORR      Becke 88 exchange correction
B88            aka BECKECORRX
BECKESR        Short range Becke 88 exchange
BECKECAM       CAM Becke 88 exchange
BR             Becke-Roussells exchange with jp dependence
LDAERF         Short-range spin-dependent LDA exchange functional
OPT            OPTX Handy & Cohen exchange
REVPBE         Revised PBE Exchange Functional
RPBE           RPBE Exchange Functional
KT             KT exchange GGA correction
PW91           Perdew-Wang 1991 GGA Exchange Functional
M05            M05 exchange
M052X          M05-2X exchange
M06            M06 exchange
M062X          M06-2X exchange
M06L           M06-L exchange
M06HF          M06-HF exchange
TPSS           TPSS original exchange functional
REVTPSS        Reviewed TPSS exchange functional
SCAN           SCAN exchange functional
RSCAN          rSCAN exchange functional
RPPSCAN        r++SCAN exchange functional
R2SCAN         r2SCAN exchange functional
R4SCAN         r4SCAN exchange functional
B97            B97 exchange
B97_1          B97-1 exchange
B97_2          B97-2 exchange
APBE           APBE Exchange Functional
BLOC           BLOC exchange functional
PBEINT         PBEint Exchange Functional
PBESOL         PBEsol Exchange Functional
TF             Thomas-Fermi Kinetic Energy Functional
BT             Borgoo-Tozer TS
VW             von Weizsaecker kinetic energy
TW             von Weizsacker Kinetic Energy Functional
=============  ========================================

Individual correlation functionals in xcfun are

=============  ========================================
Functional     Comments
-------------  ----------------------------------------
VWN3           VWN3 LDA Correlation functional
VWN5           VWN5 LDA Correlation functional
VWN            aka VWN5
PBE            PBE correlation functional
BR             Becke-Roussells correlation with jp dependence
LDAERF         Short-range spin-dependent LDA correlation functional
LDAERFC_JT     Short-range spin-unpolarized LDA correlation functional
LYP            LYP correlation
SPBE           sPBE correlation functional
VWN_PBE        PBE correlation functional using VWN LDA correlation.
PW91K          PW91 GGA Kinetic Energy Functional
PW92           PW92 LDA correlation
M052X          M05-2X Correlation
M05            M05 Correlation
M06            M06 Correlation
M06HF          M06-HF Correlation
M06L           M06-L Correlation
M062X          M06-2X Correlation
TPSS           TPSS original correlation functional
REVTPSS        Revised TPSS correlation functional
SCAN           SCAN correlation functional
RSCAN          rSCAN correlation functional
RPPSCAN        r++SCAN correlation functional
R2SCAN         r2SCAN correlation functional
R4SCAN         r4SCAN correlation functional
PZ81           PZ81 LDA correlation
P86            P86C GGA correlation
B97            B97 correlation
B97_1          B97-1 correlation
B97_2          B97-2 correlation
CS             Colle-Salvetti correlation functional
APBE           APBE correlation functional.
ZVPBESOL       zvPBEsol correlation Functional
PBEINT         PBEint correlation Functional
PBELOC         PBEloc correlation functional.
TPSSLOC        TPSSloc correlation functional
ZVPBEINT       zvPBEint correlation Functional
PW91           PW91 Correlation
=============  ========================================


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
