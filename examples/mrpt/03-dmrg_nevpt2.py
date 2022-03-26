#!/usr/bin/env python
#
# Author: Sheng Guo <shengg@princeton.edu>
#         Qiming Sun <osirpt.sun@gmail.com>
#         Seunghoon Lee <seunghoonlee89@gmail.com>
#

'''
DMRG-CASCI then DMRG-NEVPT2 calculation.

There are two NEVPT2 implementations available for DMRG Block program.  The slow
version (default) strictly follows the formula presented in JCP, 117(2002), 9138
in which the 4-particle density matrix is explictly computed.  Typically 26
orbitals is the upper limit of the slow version due to the large requirements
on the memory usage.  The fast version employs the so called MPS-pertuber
technique.  It is able to handle much larger systems, up to about 30 orbitals.
'''

from pyscf import gto, scf, mcscf, mrpt, dmrgscf
mol = gto.M(atom = [['H', (0.,0.,i-3.5)] for i in range(8)],basis = 'sto-3g',symmetry='d2h')
m = scf.RHF(mol).run()

##############################################################################
#
# FCI-based CASCI + NEVPT2.  Two roots are computed.  mc.ci holds the two CI
# wave functions.  Root ID needs to be specified for the state-specific NEVPT2
# calculation.  By default the lowest root is computed.
#
##############################################################################
mc = mcscf.CASCI(m, 4, 4)
mc.fcisolver.nroots = 2
mc.fix_spin_(shift=.5, ss=0)
mc.kernel()

ci_nevpt_e1 = mrpt.NEVPT(mc, root=0).kernel()
ci_nevpt_e2 = mrpt.NEVPT(mc, root=1).kernel()

print('FCI NEVPT correlation E = %.15g %.15g' % (ci_nevpt_e1, ci_nevpt_e2))

##################################################
#
# DMRG-NEVPT2 fast version
# Use compressed MPS as perturber functions for SC-NEVPT2.
# 4-particle density matrix is not computed.
#
##################################################

mc = mcscf.CASCI(m, 4, 4)
dmrgscf.settings.MPIPREFIX =''
mc.fcisolver = dmrgscf.DMRGCI(mol, maxM=200)
mc.fcisolver.nroots = 2
mc.kernel()

mps_nevpt_e1 = mrpt.NEVPT(mc, root=0).compress_approx(maxM=100).kernel()
mps_nevpt_e2 = mrpt.NEVPT(mc, root=1).compress_approx(maxM=100).kernel()

print('MPS NEVPT correlation E = %.15g %.15g' % (mps_nevpt_e1, mps_nevpt_e2,))