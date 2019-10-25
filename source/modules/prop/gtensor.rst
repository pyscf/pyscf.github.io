gtensor ---
**************

The :mod:`gtensor` module computes electronic g tensors for
mean-field wavefunctions.

To calculate the g-tensor of CN molecule, first run an unrestricted HF using
the module :mod:`scf`, then calculate the g-tensor using the module 
:mod:`gtensor`::

    from pyscf import gto, scf, dft
    from pyscf.prop import gtensor
    mol = gto.M(atom='''
                C 0 0 0
                N 0 0 1.1747
                ''',
                basis='ccpvdz', spin=1, charge=0, verbose=3)
    mf = scf.UHF(mol).run()
    gobj = gtensor.uhf.GTensor(mf).set(verbose=4)
    gobj.kernel()

Details about the examples can be find from here:

Examples
========

Relevant examples
:file:`examples/prop/02-g_tensor.py`


Program reference
=================

.. automodule:: pyscf.prop.gtensor

dhf
---

.. automodule:: pyscf.prop.gtensor.dhf
   :members:

uhf
---

.. automodule:: pyscf.prop.gtensor.uhf
   :members:

uks
---

.. automodule:: pyscf.prop.gtensor.uks
   :members:
