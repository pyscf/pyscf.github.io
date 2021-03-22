.. _user_cc:

**********************
Coupled-cluster theory 
**********************

*Modules*: :mod:`cc`, :mod:`pbc.cc`

The MP2 and coupled-cluster functionalities of PySCF are similar.  See
also :numref:`user_mp2`.

Introduction
============
Coupled-cluster (CC) theory is a post-Hartree-Fock method capable of describing
electron correlation in the ground state.  It is size extensive but not
variational. 
PySCF has extensive support for CC calculations with single and double excitations (CCSD).
It can also include a perturbative treatment of triple excitations (CCSD(T)),
which is a very accurate method for single-reference quantum chemistry.
CC calculations can be performed with or without density fitting,
depending on the initial SCF calculation.
Correlated excited states are
accessible through the equation-of-motion (EOM) CCSD framework, described below.

A minimal example of a CCSD and CCSD(T) calculation is as follows::

    from pyscf import gto, scf, cc
    mol = gto.M(
        atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
        basis = 'ccpvdz',
        symmetry = True,
    )
    mf = scf.HF(mol).run()
    # Note that the line following these comments could be replaced by
    # mycc = cc.CCSD(mf)
    # mycc.kernel()
    mycc = cc.CCSD(mf).run()
    print('CCSD total energy', mycc.e_tot)
    et = mycc.ccsd_t()
    print('CCSD(T) total energy', mycc.e_tot + et)

Spin symmetry
=============
The CC module in PySCF supports a number of reference wave functions with
broken spin symmetry.  In particular, CC can be performed with a
spin-restricted, spin-unrestricted, and general (spin-mixed) Hartree-Fock
solution, leading to the RCCSD, UCCSD, and GCCSD methods.

The module-level ``cc.CCSD(mf)`` constructor can infer the correct method based
on the level of symmetry-breaking in the mean-field argument.  For more explicit
control or inspection, the respective classes and functions can be found in
``ccsd.py`` (restricted with real orbitals), ``rccsd.py`` (restricted with
potentially complex orbitals), ``uccsd.py`` (unrestricted), and ``gccsd.py``
(general).

For example, a spin-unrestricted calculation on triplet oxygen can be performed
as follows::

    from pyscf import gto, scf, cc
    mol = gto.M(
        atom = 'O 0 0 0; O 0 0 1.2',  # in Angstrom
        basis = 'ccpvdz',
        spin = 2
    )
    mf = scf.HF(mol).run() # this is UHF
    mycc = cc.CCSD(mf).run() # this is UCCSD
    print('UCCSD total energy = ', mycc.e_tot)


Properties
==========

A number of properties are available at the CCSD level.

Unrelaxed 1- and 2-electron reduced density matrices can be calculated. 
They are returned in the MO basis::

    dm1 = mycc.make_rdm1()
    dm2 = mycc.make_rdm2()

Analytical nuclear gradients can be calculated::

    mygrad = mycc.nuc_grad_method().run()

The CCSD Lambda equations can be solved::

    l1, l2 = mycc.solve_lambda()


Frozen orbitals
===============

By default, CCSD calculations in PySCF correlate all electrons in all available
orbitals. To freeze the lowest-energy core orbitals,
use the ``frozen`` keyword argument::

    mycc = cc.CCSD(mf, frozen=2).run()

To freeze occupied and/or unoccupied orbitals with finer control, a list of
0-based orbital indices can be provided as the ``frozen`` keyword argument::
    
    # freeze 2 core orbitals
    mycc = cc.CCSD(mf, frozen=[0,1]).run()
    # freeze 2 core orbitals and 3 unoccupied orbitals
    mycc = cc.CCSD(mf, frozen=[0,1,16,17,18]).run()


Equation-of-motion coupled-cluster theory 
=========================================

EOM-CCSD can be used to calculate neutral excitation energies (EE-EOM-CCSD),
spin-flip excitations (SF-EOM-CCSD),
or charged excitations, i.e. ionization potentials (IP-EOM-CCSD) or electron affinities
(EA-EOM-CCSD).  The EOM functions return the requested number of 
eigenvalues and right-hand eigenvectors. For example::
    
    mycc.kernel()
    e_ip, c_ip = mycc.ipccsd(nroots=1)
    e_ea, c_ea = mycc.eaccsd(nroots=1)
    e_ee, c_ee = mycc.eeccsd(nroots=1)
    e_sf, c_sf = mycc.eomsf_ccsd(nroots=1)

The ``eecsd()`` function returns neutral excitations with all possible spin
multiplicities.  For closed-shell calculations (RHF and RCCSD), singlet and triplet 
excitations can be requested explicitly::

    e_s, c_s = mycc.eomee_ccsd_singlet(nroots=1)
    e_t, c_t = mycc.eomee_ccsd_triplet(nroots=1)

By default, PySCF calculates the ``nroots`` eigenvalues with the lowest energy,
which may include states with dominant double-excitation character.  To only
calculate states with dominant single-excitation character, use the ``koopmans``
keyword argument::

    e, c = mycc.eeccsd(nroots=3, koopmans=True)

An initial guess wavefunction may be provided, in which case PySCF will try to
find the most similar EOM solution vector::

    from pyscf.cc.eom_rccsd import amplitudes_to_vector_ee
    r1 = np.zeros((nocc,nvir))
    r2 = np.zeros((nocc,nocc,nvir,nvir))
    r1[occ_index,vir_index] = 1.0
    myguess = amplitudes_to_vector_ee(r1,r2)
    e_s, c_s = mycc.eomee_ccsd_singlet(nroots=1, guess=myguess)



Job control
===========

Saving and restarting
---------------------

To allow for future restarts, the SCF information
and the CCSD DIIS information must be saved::

    mf = scf.HF(mol)
    mf.chkfile = 'hf.chk'
    mf.kernel()

    mycc = cc.CCSD(mf)
    mycc.diis_file = 'ccdiis.h5'
    mycc.kernel()

To restart a CCSD calculation, first the molecule and SCF information must
be restored::

    mol = lib.chkfile.load_mol('hf.chk')
    mf = scf.HF(mol)
    mf.__dict__.update(lib.chkfile.load('hf.chk', 'scf'))

Next, the CCSD calculation can be restarted by using the previous 
CCSD amplitudes as the initial guess::

    mycc = cc.CCSD(mf)
    mycc.restore_from_diis_('ccdiis.h5')
    mycc.kernel(mycc.t1, mycc.t2)

Modifying DIIS
--------------

The parameters of the DIIS algorithm can be tuned in cases where
convergence is difficult.  To increase the size of the DIIS space::

    mycc = cc.CCSD(mf)
    mycc.diis_space = 10
    mycc.kernel()

By default, DIIS is activated on the first CCSD iteration.  Sometimes
it can be helpful to postpone the use of DIIS::

    mycc = cc.CCSD(mf)
    mycc.diis_start_cycle = 4
    mycc.kernel()

Integral-direct CCSD 
--------------------

In order to avoid large memory requirements, the default behavior in CCSD calculations 
is to store most two-electron integral tensors on disk.  This leads to a
potential I/O bottleneck.  For medium-sized molecules, an integral-direct
AO-driven implementation can be more efficient.  The user must manually
request an integral-direct CCSD calculation::

    mycc = cc.CCSD(mf)
    mycc.direct = True
    e_corr, t1, t2 = mycc.kernel()


