.. _pbc_cc:

pbc.cc --- PBC coupled cluster
******************************

The module :mod:`pbc.cc` carries out PBC coupled cluster calculation, with optional usage of k-point symmetry.

Examples
========

Relevant examples
:file:`examples/pbc/12-gamma_point_post_hf.py`
:file:`examples/pbc/22-k_points_ccsd.py`
:file:`examples/pbc/24-k_points_vs_gamma.py`
:file:`examples/pbc/25-k_points_mpi_ccsd.py`
:file:`examples/pbc/29-eom_ccsd_Ta.py`
:file:`examples/pbc/36-ccsd_level_shift.py`



Program reference
=================

.. .. automodule:: pyscf.pbc.cc

ccsd
----
.. .. autoclass:: pyscf.pbc.cc.rccsd.RCCSD
.. .. autoclass:: pyscf.pbc.cc.rccsd.UCCSD
.. .. autoclass:: pyscf.pbc.cc.rccsd.GCCSD

.. .. automodule:: pyscf.pbc.cc.ccsd
      :members:

kccsd
-----
.. .. autoclass:: pyscf.pbc.cc.kccsd.GCCSD
.. .. automodule:: pyscf.pbc.cc.kccsd
      :members:


kccsd_rhf
---------

.. .. autoclass:: pyscf.pbc.cc.kccsd.RCCSD
.. .. automodule:: pyscf.pbc.cc.kccsd_rhf
      :members:

kccsd_uhf
---------

.. .. autoclass:: pyscf.pbc.cc.kccsd.KUCCSD
.. .. automodule:: pyscf.pbc.cc.kccsd_uhf
      :members:

kccsd_t
-------

.. .. automodule:: pyscf.pbc.cc.kccsd_t
      :members:

kccsd_t_rhf
-----------

.. .. automodule:: pyscf.pbc.cc.kccsd_t_rhf
      :members:

kccsd_t_rhf_slow
----------------

.. .. automodule:: pyscf.pbc.cc.kccsd_t_rhf_slow
      :members:

kintermediates
--------------

.. .. automodule:: pyscf.pbc.cc.kintermediates
      :members:

kintermediates_rhf
------------------

.. .. automodule:: pyscf.pbc.cc.kintermediates_rhf
      :members:

kintermediates_uhf
------------------

.. .. automodule:: pyscf.pbc.cc.kintermediates_uhf
      :members:

eom_kccsd_ghf
-------------

.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_ghf.EOMIP
.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_ghf.EOMEA
.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_ghf.EOMIP_Ta
.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_ghf.EOMEA_Ta
.. .. automodule:: pyscf.pbc.cc.eom_kccsd_ghf
      :members:

eom_kccsd_rhf
-------------

.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_rhf.EOMIP
.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_rhf.EOMEA
.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_rhf.EOMIP_Ta
.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_rhf.EOMEA_Ta
.. .. automodule:: pyscf.pbc.cc.eom_kccsd_rhf
      :members:

eom_kccsd_rhf_ip
----------------

.. .. automodule:: pyscf.pbc.cc.eom_kccsd_rhf_ip
      :members:

eom_kccsd_rhf_ea
----------------

.. .. automodule:: pyscf.pbc.cc.eom_kccsd_rhf_ea
      :members:

eom_kccsd_uhf

-------------
.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_uhf.EOMIP
.. .. autoclass:: pyscf.pbc.cc.eom_kccsd_uhf.EOMEA
.. .. automodule:: pyscf.pbc.cc.eom_kccsd_uhf
      :members:
