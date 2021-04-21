.. _hessian:

:mod:`hessian` --- Analytical nuclear Hessian
*********************************************

The :mod:`hessian` module implements the analytical nuclear Hessian for
mean-field methods. This module also provides functions for thermo-chemical
corrections (entropic corrections) using the harmonic model obtained
from the Hessian frequencies,

For example::

    from pyscf import gto
    mol = gto.M(
        atom = [
            ['O' , 0. , 0.     , 0],
            ['H' , 0. , -0.757 , 0.587],
            ['H' , 0. ,  0.757 , 0.587]],
        basis = '631g')
    
    mf = mol.RHF().run()
    h = mf.Hessian().kernel()

The resulting structure of the Hessian is ::h[Atom_1, Atom_2, Atom_1_XYZ, Atom_1_XYZ]:: in this case ::(3,3,3,3)::.

Examples
========

* :source:`examples/hessian/01-scf_hessian.py`
* :source:`examples/hessian/10-thermochemistry.py`

Program reference
=================

.. automodule:: pyscf.hessian
 
Spin-restricted Hartree-Fock
----------------------------

.. automodule:: pyscf.hessian.rhf
   :members:

Spin-unrestricted Hartree-Fock
------------------------------

.. automodule:: pyscf.hessian.uhf
   :members:

Spin-restricted DFT
-------------------

.. automodule:: pyscf.hessian.rks
   :members:

Spin-unrestricted DFT
---------------------

.. automodule:: pyscf.hessian.uks
   :members:

Thermo-chemistry analysis
-------------------------

.. automodule:: pyscf.hessian.thermo
   :members:
