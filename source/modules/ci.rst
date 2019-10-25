.. _ci:

ci --- Configuration interaction
********************************

.. module:: ci
   :synopsis: CI energy and properties

The :mod:`ci` module implements the truncated CI model. A simple example for cisd::

    import pyscf
    mol = pyscf.M(atom = 'H 0 0 0; F 0 0 1.1',basis = 'ccpvdz')
    mf = mol.HF().run()
    mycc = mf.CISD().run()

Examples
========

:file:`00-simple_cisd.py`
:file:`01-density_matrix.py`
:file:`20-from_fci.py`
:file:`32-wfn_overlap.py`


Program reference
=================

ci.cisd
-------

.. automodule:: pyscf.ci.cisd
   :members:

ci.addons
---------

.. automodule:: pyscf.ci.addons
   :members:

ci.gcisd
--------

.. automodule:: pyscf.ci.gcisd
   :members:

ci.ucisd
--------

.. automodule:: pyscf.ci.ucisd
   :members:

