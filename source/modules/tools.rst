.. _tools:

:mod:`tools` --- Useful tools
*****************************

.. module:: tools

The :mod:`tools` module contains useful tools to dump out integrals,
orbital densities, cube files, and other things to interface with
external utilities.

For example, one can write molecular orbitals in molden format as::
    
    from pyscf import gto, scf
    from pyscf import lo
    from pyscf.tools import molden
    mol = gto.M(
        atom = '''
      C  3.2883    3.3891    0.2345
      C  1.9047    3.5333    0.2237
      C  3.8560    2.1213    0.1612
      C  1.0888    2.4099    0.1396
      C  3.0401    0.9977    0.0771
      C  1.6565    1.1421    0.0663
      H  3.9303    4.2734    0.3007
      H  1.4582    4.5312    0.2815
      H  4.9448    2.0077    0.1699
      H  0.0000    2.5234    0.1311
      H  3.4870    0.0000    0.0197
      H  1.0145    0.2578    0.0000
               ''',
        basis = 'cc-pvdz',
        symmetry = 1)
    mf = scf.RHF(mol)
    mf.kernel()
    with open('C6H6mo.molden', 'w') as f1:
        molden.header(mol, f1)
        molden.orbital_coeff(mol, f1, mf.mo_coeff, ene=mf.mo_energy, occ=mf.mo_occ)

Cube files can be generated as::
    
    from pyscf import gto, scf
    from pyscf.tools import cubegen
    mol = gto.M(atom='''
                O 0.0000000, 0.000000, 0.00000000
                H 0.761561 , 0.478993, 0.00000000
                H -0.761561, 0.478993, 0.00000000''', basis='6-31g*')
    mf = scf.RHF(mol).run()
    cubegen.density(mol, 'h2o_den.cube', mf.make_rdm1())


Examples
========

* :source:`examples/tools/01-fcidump.py`
* :source:`examples/tools/02-molden.py`
* :source:`examples/tools/03-print_mo_and_dm.py`
* :source:`examples/tools/04-analyze_local_orbitals.py`
* :source:`examples/tools/05-cubegen.py`
* :source:`examples/tools/06-chgcar.py`
* :source:`examples/tools/11-davidson_eigh.py`
* :source:`examples/tools/12-einsum.py`

Program reference
=================

.. automodule:: pyscf.tools
 
FCIDUMP
=======

.. automodule:: pyscf.tools.fcidump
   :members:


Molden
======

.. automodule:: pyscf.tools.molden
   :members:


GAMESS WFN
==========

.. automodule:: pyscf.tools.wfn_format
   :members:


Cubegen
=======

.. automodule:: pyscf.tools.cubegen
   :members:


Print Matrix
============

.. automodule:: pyscf.tools.dump_mat
   :members:

VASP CHGCAR
===========

.. automodule:: pyscf.tools.chgcar
   :members:

Molpro to PySCF
===============

.. automodule:: pyscf.tools.Molpro2Pyscf
   :members:

chkfile
==========

.. automodule:: pyscf.tools.chkfile_util
   :members:


MO mapping
==========

.. automodule:: pyscf.tools.mo_mapping
   :members:

Augmented Hessian Newton-Raphson
================================

.. automodule:: pyscf.tools.rhf_newtonraphson
   :members:

ring
==========

.. automodule:: pyscf.tools.ring
   :members:

c60struct
=========

.. automodule:: pyscf.tools.c60struct
   :members:


