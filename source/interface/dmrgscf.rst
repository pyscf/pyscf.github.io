.. _dmrgscf:

:mod:`dmrgscf` --- DMRG and DMRG-SCF/CASSCF
*******************************************

The :mod:`dmrgscf` module interfaces PySCF with external DMRG solvers, such as
`Block <https://sanshar.github.io/Block>`_, 
`StackBlock <https://github.com/sanshar/StackBlock>`_ and 
`CheMPS2 <https://github.com/SebWouters/CheMPS2>`_.
See also :mod:`mcscf`.

To perform a DMRG CASSCF calculation, first modify pyscf/dmrgscf/settings.py
and set the correct path for the DMRG solver, and then run *e.g.*, ::

    from pyscf import gto
    from pyscf import scf
    from pyscf import mcscf
    from pyscf import dmrgscf

    import os
    from pyscf.dmrgscf import settings
    if 'SLURMD_NODENAME' in os.environ:  # slurm system
        settings.MPIPREFIX = 'srun'
    elif 'PBS_NODEFILE' in os.environ:   # PBS system
        settings.MPIPREFIX = 'mpirun'
    else:  # MPI on single node
        settings.MPIPREFIX = 'mpirun -np 4'

    b = 1.2
    mol = gto.M(
        verbose = 4,
        atom = 'N 0 0 0; N 0 0 %f'%b,
        basis = 'cc-pvdz',
        symmetry = True,
    )
    mf = scf.RHF(mol)
    mf.kernel()

    mc = dmrgscf.DMRGSCF(mf, 8, 8)
    mc.state_average_([0.5, 0.5])
    mc.kernel()
    print(mc.e_tot)

  
Examples
========

* :source:`examples/dmrg/01-dmrg_casscf_with_block.py`
* :source:`examples/dmrg/01-dmrg_casscf_with_stackblock.py`
* :source:`examples/dmrg/02-dmrg_nevpt2.py`
* :source:`examples/dmrg/03-density_matrix.py`
* :source:`examples/dmrg/10-state_average.py`
* :source:`examples/dmrg/11-excited_states.py`
* :source:`examples/dmrg/30-dmrg_casscf_nevpt2_for_Cr2.py`
* :source:`examples/dmrg/31-dmrg_casscf_for_feporph.py`
* :source:`examples/dmrg/32-dmrg_casscf_nevpt2_for_FeS.py`

Program Reference
=================

.. automodule:: pyscf.dmrgscf

.. automodule:: pyscf.dmrgscf.chemps2
   :members:

.. automodule:: pyscf.dmrgscf.dmrgci
   :members:

.. automodule:: pyscf.dmrgscf.p_dmrg
   :members:

.. automodule:: pyscf.dmrgscf.nevpt_mpi
   :members:

.. automodule:: pyscf.dmrgscf.dmrg_sym
   :members:
