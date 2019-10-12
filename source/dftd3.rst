dftd3 --- DFT plus Dispersion Correction
**********************************
An interface to libdftd3:
https://github.com/cuanto/libdftd3
A minimal example for this::

    from pyscf import gto,scf,dftd3
    mol = gto.Mole()
    mol.atom = ''' O    0.00000000    0.00000000   -0.11081188
                   H   -0.00000000   -0.84695236    0.59109389
                   H   -0.00000000    0.89830571    0.52404783 '''
    mol.basis = 'sto3g'
    mol.build()
    mf = dftd3.dftd3(scf.RHF(mol))
    mf.kernel()
Examples
========

Relevant examples
:file:`examples/dftd3/00-hf_with_dftd3.py`


Program reference
=================

.. automodule:: pyscf.dftd3

itrf
------

.. automodule:: pyscf.dftd3.itrf
   :members:
