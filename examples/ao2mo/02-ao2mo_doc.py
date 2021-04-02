from pyscf import gto, scf, ao2mo
import tempfile

mol = gto.Mole()
mol.build(
    atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
    basis = 'ccpvdz',
    symmetry = True,
)

mf = scf.RHF(mol)
mf.kernel()

orb = mf.mo_coeff
# get the two-electron integrals as a numpy array
eri = ao2mo.kernel(mol, orb)

ftmp = tempfile.NamedTemporaryFile()
# saves the two-electron integrals in the file ftmp.name
ao2mo.kernel(mol, orb, ftmp.name)

ftmp = tempfile.NamedTmpFile()
ao2mo.kernel(mol, orb, ftmp)
with ao2mo.load(ftmp) as eri:
    print(eri.shape)

import numpy
with ao2mo.load(ftmp) as eri:
    eri1 = ao2mo.restore(1, numpy.asarray(eri), orb.shape[1])
    eri4 = ao2mo.restore('4', numpy.asarray(eri), orb.shape[1)
    eri8 = ao2mo.restore('s8', numpy.asarray(eri), orb.shape[1])
    print(eri1.shape)
    print(eri4.shape)
    print(eri8.shape)

eri = ao2mo.kernel(mol, orb)
print(type(eri))
with ao2mo.load(erifile) as eri:
    print(eri.shape)
    
ao2mo.kernel(mol, orb, erifile = 'hf.h5', dataname = 'test')
# load 2e integrals by filename and dataname
with ao2mo.load('hf.h5', 'test') as eri:
    print(eri.shape)
    
import numpy
erirand = numpy.random.random((5,5,5,5))
    # load 2e integrals from numpy array
with ao2mo.load(erirand) as eri:
    print(eri.shape)
