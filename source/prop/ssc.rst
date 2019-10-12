ssc --- Spin-spin coulping
*******

To compute the nuclear spin-spin coupling constants, one could follow the following example::

    from pyscf import gto, scf, dft
    from pyscf.prop import ssc
    mol = gto.M(atom='''
                O 0 0      0
                H 0 -0.757 0.587
                H 0  0.757 0.587''',
                basis='ccpvdz')
    
    mf = scf.UHF(mol).run()
    ssc.UHF(mf).kernel()
    
    mf = dft.UKS(mol).set(xc='b3lyp').run()
    ssc.UKS(mf).kernel()
    
Examples
========

Relevant examples
:file:`examples/prop/05-ssc.py`


Program reference
=================

.. automodule:: pyscf.prop.ssc

dhf
----

.. automodule:: pyscf.prop.ssc.dhf
   :members:    

rhf
----

.. automodule:: pyscf.prop.ssc.rhf
   :members:    

rks
----

.. automodule:: pyscf.prop.ssc.rks
   :members:    

uhf
----

.. automodule:: pyscf.prop.ssc.uhf
   :members:    


uks
----

.. automodule:: pyscf.prop.ssc.uks
   :members:    
