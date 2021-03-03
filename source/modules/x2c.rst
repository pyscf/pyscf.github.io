.. _x2c:

:mod:`x2c` --- exact-two-component approach
*******************************************

An example to apply scalar relativistic effects by decorating the scf object with module :mod:`x2c` is shown as following::

    from pyscf import gto
    from pyscf import scf
    mol = gto.M(
        verbose = 0,
        atom = '''8  0  0.     0
                  1  0  -0.757 0.587
                  1  0  0.757  0.587''',
        basis = 'ccpvdz',
    )
    mf = scf.RHF(mol).x2c().run()
    mol.spin = 1
    mol.charge = 1
    mol.build(0, 0)
    mf = scf.UKS(mol).x2c1e()
    energy = mf.kernel()

More examples can be find here:

Examples
========

* :source:`examples/x2c/01-spin_free_x2c.py`
* :source:`examples/x2c/02-basis_for_x.py`
* :source:`examples/x2c/10-picture_change.py`

Program reference
=================

.. automodule:: pyscf.x2c
 
X2C
---

.. automodule:: pyscf.x2c.x2c
   :members:

1e spin-free x2c
----------------

.. automodule:: pyscf.x2c.sfx2c1e
   :members:

1e spin-free x2c gradient
-------------------------

.. automodule:: pyscf.x2c.sfx2c1e_grad
   :members:

1e spin-free x2c hessian
------------------------

.. automodule:: pyscf.x2c.sfx2c1e_hess
   :members:

