.. _grad:

:mod:`grad` --- Analytical nuclear gradients
********************************************

.. module:: grad
	    
The :mod:`grad` module provides gradients for mean-field and correlated methods.
For example, the RHF gradient can be computed by::

    from pyscf import gto, scf
    mol = gto.M(
        atom = [
            ['O' , 0. , 0.     , 0],
            ['H' , 0. , -0.757 , 0.587],
            ['H' , 0. ,  0.757 , 0.587]],
            basis = '631g')
    mf = scf.RHF(mol)
    mf.kernel()
    g = mf.nuc_grad_method()
    g.kernel()

 
Examples
========

* :source:`examples/grad/01-scf_grad.py`
* :source:`examples/grad/02-dft_grad.py`
* :source:`examples/grad/03-mp2_grad.py`
* :source:`examples/grad/04-cisd_grad.py`
* :source:`examples/grad/05-ccsd_grad.py`
* :source:`examples/grad/06-tddft_gradients.py`
* :source:`examples/grad/10-excited_state_cisd_grad.py`
* :source:`examples/grad/11-excited_state_casci_grad.py`
* :source:`examples/grad/16-scan_force.py`

Program reference
=================

.. automodule:: pyscf.grad

casci
-----

.. automodule:: pyscf.grad.casci
   :members:


casscf
------

.. automodule:: pyscf.grad.casscf
   :members:
 
ccsd
----

.. automodule:: pyscf.grad.ccsd
   :members:

ccsd_slow
---------

.. automodule:: pyscf.grad.ccsd_slow
   :members:

ccsd_t
------

.. automodule:: pyscf.grad.ccsd_t
   :members:

cisd
----

.. automodule:: pyscf.grad.cisd
   :members:

dhf
---

.. automodule:: pyscf.grad.dhf
   :members:

mp2
---

.. automodule:: pyscf.grad.mp2
   :members:

rhf
---

.. automodule:: pyscf.grad.rhf
   :members:

rks
---

.. automodule:: pyscf.grad.rks
   :members:

rohf
----

.. automodule:: pyscf.grad.rohf
   :members:

roks
----

.. automodule:: pyscf.grad.roks
   :members:

tdrhf
-----

.. automodule:: pyscf.grad.tdrhf
   :members:

tdrhf_slow
----------

.. automodule:: pyscf.grad.tdrhf_slow
   :members:

tdrks
-----

.. automodule:: pyscf.grad.tdrks
   :members:

tduhf
-----

.. automodule:: pyscf.grad.tduhf
   :members:

tduks
-----

.. automodule:: pyscf.grad.tduks
   :members:

uccsd
-----

.. automodule:: pyscf.grad.uccsd
   :members:

uccsd_t
-------

.. automodule:: pyscf.grad.uccsd_t
   :members:

ucisd
-----

.. automodule:: pyscf.grad.ucisd
   :members:

uhf
---

.. automodule:: pyscf.grad.uhf
   :members:

uks
---

.. automodule:: pyscf.grad.
   :members:

ump2
----

.. automodule:: pyscf.grad.
   :members:
