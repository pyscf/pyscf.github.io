fci --- Full configuration interaction
**************************************

.. automodule:: pyscf.fci

The :mod:`fci` module implements Full Configuration Interaction. Different
implementations are available for different Hamiltonian and spin symmetries.
One example for FCI calculation on HF molecule::

    from pyscf import gto, scf, fci
    mol = gto.Mole(atom = 'H 0 0 0; F 0 0 1.1', basis= '6-31g', symmetry=True)
    myhf = scf.RHF(mol)
    myhf.kernel()
    cisolver = fci.FCI(mol, myhf.mo_coeff)
    cisolver.kernel()

Examples
========

For a complete list of FCI examples, see
``pyscf/examples/fci``.

Program reference
=================

direct CI
---------

.. automodule:: pyscf.fci.direct_spin1
   :members:

-----

.. automodule:: pyscf.fci.direct_spin1_symm
   :members:

-----

.. automodule:: pyscf.fci.direct_spin0
   :members:

-----

.. automodule:: pyscf.fci.direct_spin0_symm
   :members:

-----

.. automodule:: pyscf.fci.direct_uhf
   :members:


cistring
--------

.. automodule:: pyscf.fci.cistring
   :members:


spin operator
-------------

.. automodule:: pyscf.fci.spin_op
   :members:


rdm
---

.. automodule:: pyscf.fci.rdm
   :members:


addons
------

.. automodule:: pyscf.fci.addons
   :members:
