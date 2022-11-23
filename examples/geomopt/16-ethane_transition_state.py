from pyscf.qsdopt.qsd_optimizer import QSD
from pyscf import gto, scf

mol = gto.M(atom='''
C  -2.17177  -0.46068  -0.00000
C  -0.02194  -0.50455  0.00000
H  0.35133   0.30506   -0.59169
H  0.33677   -0.40796  1.00344
H  0.31586   -1.43259  -0.41175
H  -2.68285  -0.25342  0.93134
H  -2.40018  0.15510   -0.90287
H   -2.53568   -1.52338   -0.17790''', basis='minao', verbose=0)
mf = scf.RHF(mol)

optimizer = QSD(mf, stationary_point="TS")
optimizer.kernel(hess_update_freq=0, step=0.5, hmin=1e-3)
