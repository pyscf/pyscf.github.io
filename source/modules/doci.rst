.. _doci:

:mod:`doci` --- Doubly occupied configuration interaction
*********************************************************

.. module:: doci
The :mod:`doci` implements doubly-occupied configuration interaction. A minimal example on N2 calculation::

    from pyscf import gto, doci
    mol = gto.M(atom='N 0 0 0; N 0 0 2.', basis='6-31g')
    mf = mol.RHF().run()
    mc = doci.CASSCF(mf, 18, 14)
    mc.kernel()

Examples
========

* :source:`examples/doci/00-simple_doci_casscf.py`

Program reference
=================

.. automodule:: pyscf.doci

doci_mcscf
------

.. automodule:: pyscf.doci.doci_mcscf
   :members:

doci_slow
-------

.. automodule:: pyscf.doci.doci_slow
   :members:
