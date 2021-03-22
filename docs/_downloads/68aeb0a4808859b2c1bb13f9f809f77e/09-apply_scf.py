#!/usr/bin/env python
#
# Author: Qiming Sun <osirpt.sun@gmail.com>
#

'''
Apply SCF and post-SCF methods
'''

from pyscf import gto

mol = gto.M(atom='H 0 0 0; F 0 0 1.1', basis='ccpvdz')

#
# After initializing the molecular object, call .apply method to run SCF or
# post-SCF calculations
#
mf = mol.apply('RHF').run()

ccobj = mol.apply('CCSD').run()

tdobj = mol.apply('TDHF').run(nstates=5)
