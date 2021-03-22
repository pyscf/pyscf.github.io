.. _fci:

:mod:`fci` --- Full configuration interaction
*********************************************

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

* :source:`examples/fci/00-simple_fci.py`
* :source:`examples/fci/01-given_h1e_h2e.py`
* :source:`examples/fci/02-selected_ci.py`
* :source:`examples/fci/10-spin.py`
* :source:`examples/fci/11-large_ci.py`
* :source:`examples/fci/12-multiple_roots.py`
* :source:`examples/fci/13-wfn_symmetry.py`
* :source:`examples/fci/14-density_matrix.py`
* :source:`examples/fci/15-FCI_hamiltonian.py`
* :source:`examples/fci/30-h6_scan.py`
* :source:`examples/fci/31-apply_2nd_quantized_op.py`
* :source:`examples/fci/32-wfn_overlap.py`
* :source:`examples/fci/33-rotate_wfn.py`

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
