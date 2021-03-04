.. _theory_mp2:

Second-Order Møller–Plesset Perturbation Theory
***********************************************

*Modules*: :mod:`mp`, :mod:`pbc.mp`

The interfaces to coupled cluster and MP2 calculations are similar
in PySCF.  Therefore, note the parallels between this section and the
:ref:`theory_cc` documentation.


Introduction
============

Second-Order Møller–Plesset Perturbation Theory (MP2) :cite:`Moller1934`
is a post-Hartree--Fock method.

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
  
  for the same result.


Spin symmetry
=============

The MP2 module in PySCF supports a number of reference wavefunctions with
broken spin symmetry wavefunctions.  In particular, MP2 can be performed with a
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

Gradients (see e.g. :cite:`Pople1979`, :cite:`Handy1985` for MP2 gradients)
can be calculated::

    from pyscf import grad
    mygrad = mymp.Gradients()
    grad = mygrad.kernel()


Frozen orbitals
===============

By default, MP2 calculations in PySCF include all electrons; that is, also
core orbitals are correlated. To freeze the lowest-energy core orbitals, use
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

Saving memory
-------------
If the t2 amplitudes are not required after the MP2 calculation, they
don't need to be saved.  This can be set by setting the ``with_t2``
keyword argument::

    mymp = mp.MP2(mf)
    # set with_t2 to False to not save t2 amplitudes (default is True)
    mymp.kernel(with_t2=False)


MP2 and MO integral transformation
==================================

For the triplet
state, we can write a function to compute the correlation energy as

.. math::

  E_{\text{corr}} = \frac{1}{4}\sum_{ijab}
           \frac{\langle ij||ab \rangle \langle ab||ij \rangle}
           {\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b}.

.. code:: python

  def myump2(mf):
      import numpy
      from pyscf import ao2mo
      # As UHF objects, mo_energy, mo_occ, mo_coeff are two-item lists
      # (the first item for alpha spin, the second for beta spin).
      mo_energy = mf.mo_energy
      mo_occ = mf.mo_occ
      mo_coeff = mf.mo_coeff
      o = numpy.hstack((mo_coeff[0][:,mo_occ[0]>0] ,mo_coeff[1][:,mo_occ[1]>0]))
      v = numpy.hstack((mo_coeff[0][:,mo_occ[0]==0],mo_coeff[1][:,mo_occ[1]==0]))
      eo = numpy.hstack((mo_energy[0][mo_occ[0]>0] ,mo_energy[1][mo_occ[1]>0]))
      ev = numpy.hstack((mo_energy[0][mo_occ[0]==0],mo_energy[1][mo_occ[1]==0]))
      no = o.shape[1]
      nv = v.shape[1]
      noa = sum(mo_occ[0]>0)
      nva = sum(mo_occ[0]==0)
      eri = ao2mo.general(mf.mol, (o,v,o,v)).reshape(no,nv,no,nv)
      eri[:noa,nva:] = eri[noa:,:nva] = eri[:,:,:noa,nva:] = eri[:,:,noa:,:nva] = 0
      g = eri - eri.transpose(0,3,2,1)
      eov = eo.reshape(-1,1) - ev.reshape(-1)
      de = 1/(eov.reshape(-1,1) + eov.reshape(-1)).reshape(g.shape)
      emp2 = .25 * numpy.einsum('iajb,iajb,iajb->', g, g, de)
      return emp2

.. code:: python

  >>> print('E(UMP2) = %.9g' % myump2(uhf3))
  -0.346926068

In this example, we concatenate :math:`\alpha` and :math:`\beta` orbitals to
mimic the spin-orbitals.  After integral transformation, we zeroed out the
integrals of different spin.  Here, the :mod:`ao2mo` module provides the general
2-electron MO integral transformation, using chemists' notation for the
integrals.  Using this module, you can perform
*arbitrary* integral transformation for *arbitrary* integrals. For example, the
following code gives the ``(ov|vv)`` type integrals::

  >>> from pyscf import ao2mo
  >>> import h5py
  >>> mocc = m.mo_coeff[:,m.mo_occ>0]
  >>> mvir = m.mo_coeff[:,m.mo_occ==0]
  >>> ao2mo.general(mol, (mocc,mvir,mvir,mvir), 'tmp.h5', compact=False)
  >>> feri = h5py.File('tmp.h5')
  >>> ovvv = numpy.array(feri['eri_mo'])
  >>> print(ovvv.shape)
  (160, 400)

We pass ``compact=False`` to :func:`ao2mo.general` to prevent the
function using the permutation symmetry between the virtual-virtual pair
of ``|vv)``.  So the shape of ``ovvv`` corresponds to 8 occupied
orbitals by 20 virtual orbitals for electron 1 ``(ov|`` and 20 by 20 for
electron 2 ``|vv)``.  In the following example, we transformed the
analytical gradients of 2-electron integrals

.. math::

  \langle (\frac{\partial}{\partial R} \varphi_i) \varphi_k | \varphi_j \varphi_l \rangle
  = \int \frac{\frac{\partial\varphi_i(r_1)}{\partial R}
  \varphi_j(r_1) \varphi_k(r_2)\varphi_l(r_2)}{|r_1-r_2|} dr_1 dr_2

.. code:: python

  >>> nocc = mol.nelectron // 2
  >>> co = mf.mo_coeff[:,:nocc]
  >>> cv = mf.mo_coeff[:,nocc:]
  >>> nvir = cv.shape[1]
  >>> eri = ao2mo.general(mol, (co,cv,co,cv), intor='int2e_ip1_sph', comp=3)
  >>> eri = eri.reshape(3, nocc, nvir, nocc, nvir)
  >>> print(eri.shape)
  (3, 8, 20, 8, 20)


References
==========

.. bibliography:: ref_mp.bib
   :style: unsrt