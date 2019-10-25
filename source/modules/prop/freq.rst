freq ---
*******

The :mod:`freq` module handles vibrational frequencies for
mean-field wavefunctions.

To compute the frequency of a water molecule under restricted spin symmetry, 
first perform a restricted Kohn-Sham calculation using density functional 
theory module :mod:`dft`, then calculate the frequency using the restricted
Kohn-Sham module :mod:`rks`::

    from pyscf import gto, dft
    from pyscf.prop.freq import rks
    mol = gto.M(atom='''
                O 0 0      0
                H 0 -0.757 0.587
                H 0  0.757 0.587''',
                basis='ccpvdz', verbose=4)
    mf = dft.RKS(mol).run()
    w, modes = rks.Freq(mf).kernel()

Examples
========

Relevant examples
:file:`examples/prop/01-freq.py`

.. automodule:: pyscf.prop.freq


Program reference
=================
       

rks
------

.. automodule:: pyscf.prop.freq.rks
         :members:

uks
------

.. automodule:: pyscf.prop.freq.uks
         :members:


rhf
------

.. automodule:: pyscf.prop.freq.rhf
         :members:


uhf
------

.. automodule:: pyscf.prop.freq.uhf
         :members:
