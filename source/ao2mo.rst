ao2mo --- Integral transformations
**********************************

.. module:: ao2mo
   :synopsis: transformations from AO to MO integrals for
      various permutational and spin symmetries.

The :mod:`ao2mo` module implements transformations from AO to MO integrals for various permutational and spin symmetries. A simple example to perform integral transformation::
    mol = gto.Mole(atom='H 0 0 0; F 0 0 1.1', basis = 'ccpvdz')
    myhf = scf.RHF(mol)
    myhf.kernel()
    orb = myhf.mo_coeff
    eri_4fold = ao2mo.kernel(mol, orb)

Examples
========

Relevant examples
:file:`examples/ao2mo/00-mo_integrals.py`
:file:`examples/ao2mo/01-outcore.py`
:file:`examples/ao2mo/10-diff_orbs_for_ijkl.py`
:file:`examples/ao2mo/11-ump2.py`
:file:`examples/ao2mo/20-eri_grad_hess.py`
:file:`examples/ao2mo/21-spin_orbit_coupling.py`


Program reference
=================

.. automodule:: pyscf.ao2mo
 
incore
------

.. automodule:: pyscf.ao2mo.incore
   :members:


semi_incore
------

.. automodule:: pyscf.ao2mo.semi_incore
   :members:


outcore
-------

.. automodule:: pyscf.ao2mo.outcore
   :members:


addons
------

.. automodule:: pyscf.ao2mo.addons
   :members:

