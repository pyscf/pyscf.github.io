
.. _user_qmmm:

QM/MM methods
*************

*Modules*: :mod:`qmmm`

Introduction
============

In the Quantum Mechanics - Molecular Mechanics (QM/MM) method,
the system is divided into a QM region and an MM region.
QM/MM models mainly differ in the treatment of the interaction
between QM and MM regions.

The :mod:`qmmm` module implements the electronic embedding model,
where classical point charges in the MM region are
added into the QM Hamiltonian analogously to the nuclear charges:

.. math::
    \mathbf{H}_{\mathrm{QMMM}} =&\ \mathbf{H}_{\mathrm{QM}}
        + \mathbf{V}_{\mathrm{nuc}-\mathrm{MM}} + \mathbf{V}_{\mathrm{elec}-\mathrm{MM}} \\
        =&\ \mathbf{H}_{\mathrm{QM}} + \sum_{A}^{N_{\mathrm{MM}}} \sum_{B}^{N_{\mathrm{QM}}}
            \frac{Q_{A}Q_{B}}{|\mathbf{R}_A - \mathbf{R}_B|}
        - \sum_{A}^{N_{\mathrm{MM}}} \sum_i^{N_{\mathrm{elec}}}
            \frac{Q_{A}}{|\mathbf{R}_A - \mathbf{r}_i|}

The total energy calculated using the above model includes
the regular QM energy, the interaction between
the nuclei in QM region and the MM charges, and the static Coulomb
interaction between the electron density and the MM charges. It does not
include the static Coulomb interactions of the MM point charges, the MM
energy, the vdw interaction or other bonding/non-bonding effects between
QM region and MM particles.

Another QM/MM model - polarizable embedding model - is provided in the
:mod:`solvent` module. See :ref:`user_solvent` for details.

SCF methods with MM charges
===========================

MM background charges support SCF methods
:func:`scf.RHF`, :func:`scf.UHF`, :func:`scf.ROHF`, 
:func:`scf.RKS`, :func:`scf.UKS` and :func:`scf.ROKS`, 
by decorating the underlying SCF objects with :func:`qmmm.mm_charge`.
A minimal example of using the :mod:`qmmm` module is as follows. ::

    >>> from pyscf import gto, scf, qmmm
    >>> mol = gto.M(atom='H 0 0 0; F 0 0 1', basis='ccpvdz')
    >>> coords = [(0.5,0.6,0.8)]
    >>> charges = [-0.3]
    >>> mf = qmmm.mm_charge(scf.RHF(mol), coords, charges)
    >>> mf.kernel()
    converged SCF energy = -100.045455504517

In the above example, the coordinates (in the same unit as ``mol.unit``) and
charges of the point charges in the MM region are given by ``coords`` and ``charges``,
respectively.

.. note::

    Currently, MM charges do not support :func:`scf.GHF` and :func:`scf.GKS`.

Analytical nuclear gradients are calculated with the background charges. ::

    >>> mf.nuc_grad_method().run()
    --------------- QMMM gradients ---------------
         x                y                z
    0 H    -0.0157686538    -0.0189223846    -0.1102601870
    1 F    -0.0830715173    -0.0996858207     0.1182587572
    ----------------------------------------------

.. note::

    The gradients obtained from ::

    >>> from pyscf import grad
    >>> grad.RHF(mf).run()

    for a :func:`qmmm.mm_charge` decorated ``mf`` object
    will be missing the contributions from the background charges.

If MM charges and X2C correction are used together, function
:func:`qmmm.mm_charge` needs to be applied after X2C decoration. ::

    >>> qmmm.mm_charge(scf.RHF(mol).x2c(), coords, charges).run()
    converged SCF energy = -100.126131355203
    >>> qmmm.mm_charge(scf.RHF(mol).x2c1e(), coords, charges).run()
    converged SCF energy = -100.126131355203
    >>> qmmm.mm_charge(scf.RHF(mol).sfx2c1e(), coords, charges).run()
    converged SCF energy = -100.126131355203

.. note::

    X2C gradients with MM charges are not supported.

MM charges can also be used together with second order scf and solvation models.

Post-SCF methods with MM charges
================================

Once function :func:`qmmm.mm_charge` is
applied on the SCF object, it affects all the
post-HF calculations, eg. MP2, CCSD, MCSCF, etc. ::

    >>> from pyscf import mcscf
    >>> mf = qmmm.mm_charge(scf.RHF(mol), coords, charges).run()
    >>> mc = mcscf.CASSCF(mf, 4, 4)
    >>> mc.run(conv_tol=1E-10)
    CASSCF energy = -100.101848457578
    CASCI E = -100.101848457578  E(CI) = -6.74400107375546  S^2 = 0.0000000
