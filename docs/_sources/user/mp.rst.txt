.. _user_mp2:

Second-Order Møller–Plesset Perturbation Theory (MP2)
*****************************************************

*Modules*: :mod:`mp`, :mod:`pbc.mp`

The MP2 and coupled-cluster functionalities of PySCF are similar.  See
also :numref:`user_cc`.

Introduction
============

Second-order Møller–Plesset perturbation theory (MP2) :cite:`Moller1934`
is a post-Hartree--Fock method.
MP2 calculations can be performed with or without density fitting,
depending on the initial SCF calculation.

A simple example (see :source:`examples/mp/00-simple_mp2.py`) of running
an MP2 calculation is

.. literalinclude:: ../../examples/mp/00-simple_mp2.py

which outputs

.. code::

  converged SCF energy = -99.9873974403487
  E(MP2) = -100.198764900659  E_corr = -0.211367460310054

namely, the Hartree--Fock energy, the MP2 energy and their difference, the
MP2 correlation energy.

.. note::

  The last line in the code example above could have been replaced by

  .. code ::
    
    pyscf.mp.MP2(mf).kernel()

  or

  .. code ::
    
    pyscf.mp.MP2(mf).run()
  
  for the same result.


Spin symmetry
=============

The MP2 module in PySCF supports a number of reference wavefunctions with
broken spin symmetry.  In particular, MP2 can be performed with a
spin-restricted, spin-unrestricted, and general (spin-mixed) Hartree-Fock
solution, leading to the RMP2, UMP2, and GMP2 methods.

The module-level ``mp.MP2(mf)`` constructor can infer the correct method based
on the level of symmetry-breaking in the mean-field argument.  For more explicit
control or inspection, the respective classes and functions can be found in
``mp2.py`` (restricted), ``ump2.py`` (unrestricted), and ``gmp2.py``
(general).

For example, a spin-unrestricted calculation on triplet oxygen can be performed
as follows::

    from pyscf import gto, scf, mp
    mol = gto.M(
        atom = 'O 0 0 0; O 0 0 1.2',  # in Angstrom
        basis = 'ccpvdz',
        spin = 2
    )
    mf = scf.HF(mol).run() # this is UHF
    mymp = mp.MP2(mf).run() # this is UMP2
    print('UMP2 total energy = ', mymp.e_tot)


Properties
==========

A number of properties are available at the MP2 level.

Unrelaxed 1- and 2-electron reduced density matrices can be calculated. 
They are returned in the MO basis::

    dm1 = mymp.make_rdm1()
    dm2 = mymp.make_rdm2()

Analytical nuclear gradients can be calculated :cite:`Pople1979,Handy1985,Yamaguchi2011m` ::

    mygrad = mymp.nuc_grad_method().run()


Frozen orbitals
===============

By default, MP2 calculations in PySCF correlate all electrons in all available
orbitals. To freeze the lowest-energy core orbitals, use
the ``frozen`` keyword argument::

    mymp = mp.MP2(mf, frozen=2).run()

To freeze occupied and/or unoccupied orbitals with finer control, a
list of 0-based orbital indices can be provided as the ``frozen``
keyword argument::
    
    # freeze 2 core orbitals
    mymp = mp.MP2(mf, frozen=[0,1]).run()
    # freeze 2 core orbitals and 3 unoccupied orbitals
    mymp = mp.MP2(mf, frozen=[0,1,16,17,18]).run()


Job control
===========

Avoid t2 storage
----------------
If the t2 amplitudes are not required after the MP2 calculation, they
don't need to be saved in memory::

    mymp = mp.MP2(mf)
    # by default, with_t2=True 
    mymp.kernel(with_t2=False)


References
==========

.. bibliography:: ref_mp.bib
   :style: unsrt
