.. _theory_mp2:

Second-Order Møller–Plesset Perturbation Theory
***********************************************

*Modules*: :mod:`mp`, :mod:`pbc.mp`



MP2 and MO integral transformation
----------------------------------
For the triplet
state, we can write a function to compute the correlation energy

.. math::

  E_{corr} = \frac{1}{4}\sum_{ijab}
           \frac{\langle ij||ab \rangle \langle ab||ij \rangle}
           {\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b}

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
2-electron MO integral transformation.  Using this module, you are able to do
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
