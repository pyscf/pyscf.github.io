hfc --- Hyperfine coupling tensor
*********************************

To compute the hyperfine coupling tensor, one could follow the following example::
    
    from pyscf import gto, scf, dft
    from pyscf.prop import hfc
    mol = gto.M(atom='''
                C 0 0 0
                N 0 0 1.1747
                ''',
                basis='ccpvdz', spin=1, charge=0, verbose=3)
    mf = scf.UHF(mol).run()
    gobj = hfc.uhf.HFC(mf).set(verbose=4)
    gobj.para_soc2e = 'SSO+SOO'
    gobj.so_eff_charge = False
    gobj.kernel()
    
Further examples can be found here:

Examples
========

Relevant examples
:file:`examples/prop/03-hfc.py`

Program reference
=================

.. automodule:: pyscf.prop.hfc

dhf
---

.. automodule:: pyscf.prop.hfc.dhf
   :members:    


uhf
---

.. automodule:: pyscf.prop.hfc.uhf
   :members:    


uks
---

.. automodule:: pyscf.prop.hfc.uks
   :members:    
