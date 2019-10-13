.. _dmrgscf:

dmrgscf
*******

An interface to DMRG and DMRG-SCF/CASSCF. The DMRG calculation
 must be carried out with an external solver, such as
 Block (https://sanshar.github.io/Block) or StackBlock (https://github.com/sanshar/StackBlock) or CheMPS2 (https://github.com/SebWouters/CheMPS2).
 See also :mod:`mcscf`.

To perform a DMRG CASSCF calculation, first modify the pyscf/dmrgscf/settings.py
and set the correct path for DMRG solver, and then run e.g.
`examples/dmrg/01-dmrg_casscf_with_block.py`
  
Examples
========
Relevant examples
:file:`examples/dmrg/01-dmrg_casscf_with_block.py`
:file:`examples/dmrg/01-dmrg_casscf_with_stackblock.py`
:file:`examples/dmrg/02-dmrg_nevpt2.py`
:file:`examples/dmrg/03-density_matrix.py`
:file:`examples/dmrg/10-state_average.py`
:file:`examples/dmrg/11-excited_states.py`
:file:`examples/dmrg/30-dmrg_casscf_nevpt2_for_Cr2.py`
:file:`examples/dmrg/31-dmrg_casscf_for_feporph.py`
:file:`examples/dmrg/32-dmrg_casscf_nevpt2_for_FeS.py`

Program Reference
=================

.. automodule:: pyscf.dmrgscf.chemps2
   :members:
.. automodule:: pyscf.dmrgscf.dmrg_sym
   :members:
.. automodule:: pyscf.dmrgscf.dmrgci
   :members:
.. automodule:: pyscf.dmrgscf.p_dmrg
   :members:
.. automodule:: pyscf.dmrgscf.nvevpt_mpi
   :members:
