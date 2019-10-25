nmr --- NMR
*******

To compute NMR shielding constants, one could follow the following example::

    from pyscf import gto, dft
    from pyscf.prop import nmr
    mol = gto.M(atom='''
                C 0 0 0
                O 0 0 1.1747
                ''',
                basis='ccpvdz', verbose=3)
    mf = dft.RKS(mol)
    mf.xc = 'b3lyp'
    mf.run()
    nmr.RKS(mf).kernel()

Examples
========

Relevant examples
:file: `examples/prop/04-nmr.py`


Program reference
=================

.. automodule:: pyscf.prop.nmr

dhf
----

.. automodule:: pyscf.prop.nmr.dhf
   :members:    

rhf
----

.. automodule:: pyscf.prop.nmr.rhf
   :members:    

rks
----

.. automodule:: pyscf.prop.nmr.rks
   :members:    

uhf
----

.. automodule:: pyscf.prop.nmr.uhf
   :members:    


uks
----

.. automodule:: pyscf.prop.nmr.uks
   :members:    
