.. _sgx:

:mod:`sgx` --- Pseudo-spectral methods (COSX, PS, SN-K)
*******************************************************

.. module:: sgx

The :mod:`sgx` module implements pseudo-spectral methods for Hartree-Fock exchange.

Here is a simple example of how to use pseudo spectral integrals in SCF calculation::

    from pyscf import gto
    from pyscf import scf
    from pyscf import sgx
    mol = gto.M(
        atom='''O    0.   0.       0.
                H    0.   -0.757   0.587
                H    0.   0.757    0.587
        ''',
        basis = 'ccpvdz',
    )
    mf = sgx.sgx_fit(scf.RHF(mol))
    mf.kernel()
    # Using RI for Coulomb matrix while K-matrix is constructed with COS-X method
    mf.with_df.dfj = True
    mf.kernel()



Examples
========

* :source:`examples/sgx/00-simple_sgx.py`

Program reference
=================

.. automodule:: pyscf.sgx

Main class
----------

.. automodule:: pyscf.sgx.sgx
   :members:

Get JK
------

.. automodule:: pyscf.sgx.sgx_jk
   :members:

