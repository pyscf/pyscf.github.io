qmmm --- QM/MM interface
************************

To run HF with background charges, one could follow the following example::

    import numpy
    from pyscf import gto, scf, qmmm
    mol = gto.M(atom='''
    C       1.1879  -0.3829 0.0000
    C       0.0000  0.5526  0.0000
    O       -1.1867 -0.2472 0.0000
    H       -1.9237 0.3850  0.0000
    H       2.0985  0.2306  0.0000
    H       1.1184  -1.0093 0.8869
    H       1.1184  -1.0093 -0.8869
    H       -0.0227 1.1812  0.8852
    H       -0.0227 1.1812  -0.8852
                ''',
                basis='3-21g',
                verbose=4)
    numpy.random.seed(1)
    coords = numpy.random.random((5,3)) * 10
    charges = (numpy.arange(5) + 1.) * -.1
    mf = scf.UHF(mol)
    mf = qmmm.mm_charge(mf, coords, charges)
    mf.run()
    
Examples
========

Relevant examples
:file:`examples/qmmm/00-hf.py`
:file:`examples/qmmm/01-dft.py`
:file:`examples/qmmm/02-mcscf.py`
:file:`examples/qmmm/03-ccsd.py`
:file:`examples/qmmm/04-cisd.py`
:file:`examples/qmmm/05-mp2.py`
:file:`examples/qmmm/06-tddft.py`
:file:`examples/qmmm/10-x2c.py`
:file:`examples/qmmm/11-newton.py`
:file:`examples/qmmm/11-soscf.py`
:file:`examples/qmmm/20-grad.py`
:file:`examples/qmmm/21-geom_opt.py`

.. automodule:: pyscf.qmmm.itrf
   :members:

The :mod:`qmmm` module implements a generic interface for use with MM
programs.
