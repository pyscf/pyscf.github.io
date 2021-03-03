.. _lo:

:mod:`lo` --- Orbital localization and analysis tools 
*****************************************************

The :mod:`lo` module implements various orbital localizations, such
as intrinsic atomic orbitals and natural atomic orbitals. 

For example, to obtain the natural atomic orbital coefficients (in terms
of the original atomic orbitals)::

    import numpy
    from pyscf import gto, scf, lo
    
    x = .63
    mol = gto.M(atom=[['C', (0, 0, 0)],
                      ['H', (x ,  x,  x)],
                      ['H', (-x, -x,  x)],
                      ['H', (-x,  x, -x)],
                      ['H', ( x, -x, -x)]],
                basis='ccpvtz')
    mf = scf.RHF(mol).run()
    
    # C matrix stores the AO to localized orbital coefficients
    C = lo.orth_ao(mf, 'nao')
        

Examples
========

* :source:`examples/local_orb/01-pop_with_meta_lowdin.py`
* :source:`examples/local_orb/01-pop_with_nao.py`
* :source:`examples/local_orb/02-pop_with_iao.py`
* :source:`examples/local_orb/03-split_localization.py`
* :source:`examples/local_orb/04-ibo_benzene_cubegen.py`
* :source:`examples/local_orb/05-ibo_periodic_diamond_cubegen.py`
* :source:`examples/local_orb/06-vvo_livvo_water_cubegen.py`
* :source:`examples/local_orb/07-pipek_mezey.py`
* :source:`examples/local_orb/10-modify_valence_space_for_meta_lowdin.py`
* :source:`examples/local_orb/40-hubbard_model_PM_localization.py`
* :source:`examples/local_orb/nlocal.py`
* :source:`examples/local_orb/pmloc.py`
* :source:`examples/local_orb/ulocal.py`

Program reference
=================

.. automodule:: pyscf.lo


Foster-Boys, Edmiston-Ruedenberg, Pipek-Mezey localization
==========================================================

.. automodule:: pyscf.lo.boys
   :members:

.. automodule:: pyscf.lo.edmiston
   :members:

.. automodule:: pyscf.lo.pipek
   :members:


Meta-Lowdin
===========
.. automodule:: pyscf.lo.orth
   :members:


Natural atomic orbitals
=======================

.. automodule:: pyscf.lo.nao
   :members:

Intrinsic Atomic Orbitals
=========================

.. automodule:: pyscf.lo.iao
   :members:

Intrinsic Bond Orbitals
=========================

.. automodule:: pyscf.lo.ibo
   :members:
