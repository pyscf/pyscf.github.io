.. _pbc_cc:

pbc.cc --- Coupled cluster with PBCs
************************************

The :mod:`pbc.cc` module implements the coupled cluster methods with periodic boundary conditions.
It is analogous to the :mod:`cc` module.


Examples
========

* :source:`examples/pbc/12-gamma_point_post_hf.py`
* :source:`examples/pbc/22-k_points_ccsd.py`
* :source:`examples/pbc/24-k_points_vs_gamma.py`
* :source:`examples/pbc/25-k_points_mpi_ccsd.py`
* :source:`examples/pbc/29-eom_ccsd_Ta.py`
* :source:`examples/pbc/36-ccsd_level_shift.py`

Program reference
=================

.. automodule:: pyscf.pbc.cc

.. automodule:: pyscf.pbc.cc.ccsd
   :members:

.. automodule:: pyscf.pbc.cc.kccsd
   :members:

.. automodule:: pyscf.pbc.cc.kccsd_rhf
   :members:

.. automodule:: pyscf.pbc.cc.kccsd_uhf
   :members:

.. automodule:: pyscf.pbc.cc.kccsd_t
   :members:

.. automodule:: pyscf.pbc.cc.kccsd_t_rhf
   :members:

.. automodule:: pyscf.pbc.cc.kccsd_t_rhf_slow
   :members:

.. automodule:: pyscf.pbc.cc.kintermediates
   :members:

.. automodule:: pyscf.pbc.cc.kintermediates_rhf
   :members:

.. automodule:: pyscf.pbc.cc.kintermediates_uhf
   :members:

.. automodule:: pyscf.pbc.cc.eom_kccsd_ghf
   :members:

.. automodule:: pyscf.pbc.cc.eom_kccsd_rhf
   :members:

.. automodule:: pyscf.pbc.cc.eom_kccsd_rhf_ip
   :members:

.. automodule:: pyscf.pbc.cc.eom_kccsd_rhf_ea
   :members:

.. automodule:: pyscf.pbc.cc.eom_kccsd_uhf
   :members:
