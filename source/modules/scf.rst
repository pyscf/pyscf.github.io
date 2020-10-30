.. _scf:

:mod:`scf` --- Self-consistent field methods
********************************************

.. module:: scf
   :synopsis: restricted and unrestricted, closed shell and open shell Hartree-Fock methods

The :mod:`scf` module implements restricted and unrestricted, closed shell and open shell Hartree-Fock methods.


Examples
========

:download:`examples/scf/00-simple_hf.py </../examples/scf/00-simple_hf.py>`
:download:`examples/scf/01-h2o.py </../examples/scf/01-h2o.py>`
:download:`examples/scf/02-ghf.py </../examples/scf/02-ghf.py>`
:download:`examples/scf/02-rohf_uhf.py </../examples/scf/02-rohf_uhf.py>`
:download:`examples/scf/03-level_shift.py </../examples/scf/03-level_shift.py>`
:download:`examples/scf/04-dirac_hf.py </../examples/scf/04-dirac_hf.py>`
:download:`examples/scf/05-breit_gaunt.py </../examples/scf/05-breit_gaunt.py>`
:download:`examples/scf/10-glycine.py </../examples/scf/10-glycine.py>`
:download:`examples/scf/11-linear_dep.py </../examples/scf/11-linear_dep.py>`
:download:`examples/scf/12-fast_newton.py </../examples/scf/12-fast_newton.py>`
:download:`examples/scf/13-symmetry.py </../examples/scf/13-symmetry.py>`
:download:`examples/scf/14-restart.py </../examples/scf/14-restart.py>`
:download:`examples/scf/15-initial_guess.py </../examples/scf/15-initial_guess.py>`
:download:`examples/scf/16-h2_scan.py </../examples/scf/16-h2_scan.py>`
:download:`examples/scf/17-stability.py </../examples/scf/17-stability.py>`
:download:`examples/scf/20-density_fitting.py </../examples/scf/20-density_fitting.py>`
:download:`examples/scf/21-x2c.py </../examples/scf/21-x2c.py>`
:download:`examples/scf/22-newton.py </../examples/scf/22-newton.py>`
:download:`examples/scf/23-decorate_scf.py </../examples/scf/23-decorate_scf.py>`
:download:`examples/scf/24-callback.py </../examples/scf/24-callback.py>`
:download:`examples/scf/24-tune_diis.py </../examples/scf/24-tune_diis.py>`
:download:`examples/scf/30-scan_pes.py </../examples/scf/30-scan_pes.py>`
:download:`examples/scf/31-cr_atom_rohf_tune_init_guess.py </../examples/scf/31-cr_atom_rohf_tune_init_guess.py>`
:download:`examples/scf/31-v_atom_rohf.py </../examples/scf/31-v_atom_rohf.py>`
:download:`examples/scf/32-break_spin_symm.py </../examples/scf/32-break_spin_symm.py>`
:download:`examples/scf/40-apply_electric_field.py </../examples/scf/40-apply_electric_field.py>`
:download:`examples/scf/40-customizing_hamiltonian.py </../examples/scf/40-customizing_hamiltonian.py>`
:download:`examples/scf/41-hf_with_given_densityfit_ints.py </../examples/scf/41-hf_with_given_densityfit_ints.py>`
:download:`examples/scf/42-remove_linear_dep.py </../examples/scf/42-remove_linear_dep.py>`
:download:`examples/scf/43-custom_get_jk.py </../examples/scf/43-custom_get_jk.py>`
:download:`examples/scf/44-soc_ecp.py </../examples/scf/44-soc_ecp.py>`
:download:`examples/scf/50-mom-deltaSCF.py </../examples/scf/50-mom-deltaSCF.py>`
:download:`examples/scf/51-elecoup_mom.py </../examples/scf/51-elecoup_mom.py>`
:download:`examples/scf/52-dynamically_control_level_shift.py </../examples/scf/52-dynamically_control_level_shift.py>`
:download:`examples/scf/53-dynamic_sz.py </../examples/scf/53-dynamic_sz.py>`
:download:`examples/scf/54-fractional_occupancy.py </../examples/scf/54-fractional_occupancy.py>`
:download:`examples/scf/55-overload_convergence_criteria.py </../examples/scf/55-overload_convergence_criteria.py>`
:download:`examples/scf/56-h2_symm_breaking.py </../examples/scf/56-h2_symm_breaking.py>`
:download:`examples/scf/70-background_and_hot_tuning.py </../examples/scf/70-background_and_hot_tuning.py>`
:download:`examples/scf/71-patch_scf_kernel.py </../examples/scf/71-patch_scf_kernel.py>`


Addons
======

Special treatments may be required to the SCF methods in some situations.  These
special treatments cannot be universally applied for all SCF models.  They were
defined in the :mod:`scf.addons` module. For example, in an UHF calculation, we
may want the :math:`S_z` value to be changed (the numbers of alpha and beta
electrons not conserved) during SCF iteration while conserving the total number
of electrons.  :func:`scf.addons.dynamic_sz_` can provide this functionality::

    from pyscf import gto, scf
    mol = gto.M(atom='O 0 0 0; O 0 0 1')
    mf = scf.UHF(mol)
    mf.verbose=4
    mf = scf.addons.dynamic_sz_(mf)
    mf.kernel()
    print('S^2 = %s, 2S+1 = %s' % mf.spin_square())

This function automatically converges the ground sate of oxygen molecule to
triplet state although we didn't specify spin state in the :attr:`mol` object.

.. note:: Function :func:`scf.addons.dynamic_sz_` has side effects.  It changes
  the underlying mean-field object.

The `addons` mechanism increases the flexibility of PySCf program.  You can
define various addons to customize the default behaviour of pyscf program.  For
example, if you'd like to track the changes of the density (the diagonal
term of density matrix) of certain basis during the SCF iteration, you can write
the following addon to output the required density::

    def output_density(mf, basis_label):
        ao_labels = mf.mol.ao_labels()
        old_make_rdm1 = mf.make_rdm1
        def make_rdm1(mo_coeff, mo_occ):
            dm = old_make_rdm1(mo_coeff, mo_occ)
            print('AO         alpha             beta')
            for i,s in enumerate(ao_labels):
                if basis_label in s:
                    print(s, dm[0][i,i], dm[1][i,i])
            return dm
        mf.make_rdm1 = make_rdm1
        return mf
    from pyscf import gto, scf
    mol = gto.M(atom='O 0 0 0; O 0 0 1')
    mf = scf.UHF(mol)
    mf.verbose=4
    mf = scf.addons.dynamic_sz_(mf)
    mf = output_density(mf, 'O 2p')
    mf.kernel()


Caching two-electron integrals
==============================

When memory is enough (specified by the :attr:`max_memory` of SCF object), the
SCF object generates all two-electron integrals in memory and cache them in
:attr:`_eri` array.  The default :attr:`max_memory` (defined in
:data:`lib.parameters.MAX_MEMORY`, see :ref:`max_mem`) is 4 GB.  It roughly
corresponds to two-electron real integrals for 250 orbitals.  For small systems,
the cached integrals usually provide the best performance.  If you have enough
main memory in your computer, you can increase the :attr:`max_memory` of SCF
object to cache the integrals in memory.

The cached integrals :attr:`_eri` are treated as a dense tensor.  When system
becomes larger and the two-electron integral tensor becomes sparse, caching
integrals may lose performance advantage.  This is mainly due to the fact that
the implementation of J/K build for the cached integrals did not utilize the
sparsity of the integral tensor.  Also, the data locality was not considered in
the implementation which sometimes leads to bad OpenMP multi-threading speed up.
For large system, the AO-driven direct SCF method is more favorable.


.. _customize_h:

Customizing Hamiltonian
=======================

This integral object :attr:`_eri` is not merely used by the mean-field
calculation.  Along with the :meth:`get_hcore` method, this two-electron
integral object will be treated as the Hamiltonian in the post-SCF code whenever
possible.  This mechanism provides a way to model arbitrary fermion system in
PySCF.  You can customize a system by changing the 1-electron Hamiltonian and
the mean-field :attr:`_eri` attribute.  For example, the following code solves a
model system::

    import numpy
    from pyscf import gto, scf, ao2mo, ccsd
    mol = gto.M()
    n = 10
    mol.nelectron = 10
    mf = scf.RHF(mol)
    
    t = numpy.zeros((n,n))
    for i in range(n-1):
        t[i,i+1] = t[i+1,i] = -1.0
    t[n-1,0] = t[0,n-1] = 1.0  # anti-PBC
    eri = numpy.zeros((n,n,n,n))
    for i in range(n):
        eri[i,i,i,i] = 4.0
    
    mf.get_hcore = lambda *args: t
    mf.get_ovlp = lambda *args: numpy.eye(n)
    # ao2mo.restore(8, eri, n) to get 8-fold symmetry of the integrals
    # ._eri only supports the 2-electron integrals in 4-fold or 8-fold symmetry.
    mf._eri = ao2mo.restore(8, eri, n)
    
    mf.kernel()

    mycc = ccsd.RCCSD(mf).run()
    e,v = mycc.ipccsd(nroots=3)
    print('IP = ', e)
    e,v = mycc.eaccsd(nroots=3)
    print('EA = ', e)

Some post-SCF methods require the 4-index MO integrals.  Depending the available
memory (affected by the value of :attr:`max_memory` in each class), these
methods may not use the "AO integrals" cached in :attr:`_eri`.  To ensure the
post mean-field methods to use the :attr:`_eri` integrals no matter whether the
actual memory usage is over the :attr:`max_memory` limite, you can set the flag
:attr:`incore_anyway` in :class:`Mole` class to ``True`` before calling the
:meth:`kernel` function of the post-SCF methods.  In the following example,
without setting ``incore_anyway=True``, the CCSD calculations may crash::

    import numpy
    from pyscf import gto, scf, ao2mo, ccsd
    mol = gto.M()
    n = 10
    mol.nelectron = n
    mol.max_memory = 0
    mf = scf.RHF(mol)
    
    t = numpy.zeros((n,n))
    for i in range(n-1):
        t[i,i+1] = t[i+1,i] = 1.0
    t[n-1,0] = t[0,n-1] = -1.0
    eri = numpy.zeros((n,n,n,n))
    for i in range(n):
        eri[i,i,i,i] = 4.0
    
    mf.get_hcore = lambda *args: t
    mf.get_ovlp = lambda *args: numpy.eye(n)
    mf._eri = ao2mo.restore(8, eri, n)
    mf.kernel()

    mol.incore_anyway = True
    mycc = ccsd.RCCSD(mf).run()
    e,v = mycc.ipccsd(nroots=3)
    print('IP = ', e)
    e,v = mycc.eaccsd(nroots=3)
    print('EA = ', e)

Holding the entire two-particle interactions matrix elements in memory often
leads to high memory usage.  In the SCF calculation, the memory usage can be
optimized if :attr:`_eri` is sparse.  The SCF iterations requires only the Fock
matrix which in turn calls the J/K build function :meth:`SCF.get_jk` to compute
the Coulomb and HF-exchange matrix.  Overwriting the :meth:`SCF.get_jk` function
can reduce the memory footprint of the SCF part in the above example::

    import numpy
    from pyscf import gto, scf, ao2mo, ccsd
    mol = gto.M()
    n = 10
    mol.nelectron = n
    mol.max_memory = 0
    mf = scf.RHF(mol)
    
    t = numpy.zeros((n,n))
    for i in range(n-1):
        t[i,i+1] = t[i+1,i] = 1.0
    t[n-1,0] = t[0,n-1] = -1.0
    
    mf.get_hcore = lambda *args: t
    mf.get_ovlp = lambda *args: numpy.eye(n)
    def get_jk(mol, dm, *args):
        j = numpy.diag(dm.diagonal()) * 4.
        k = numpy.diag(dm.diagonal()) * 4.
        return j, k
    mf.get_jk = get_jk
    mf.kernel()

Another way to handle the two-particle interactions of large model system is to
use the density fitting/Cholesky decomposed integrals.  See also
:ref:`sl_cderi`.


Program reference
=================

.. automodule:: pyscf.scf
   :members: 

Non-relativistic Hartree-Fock
-----------------------------

.. automodule:: pyscf.scf.hf
   :members:

.. automodule:: pyscf.scf.hf_symm
   :members:

.. automodule:: pyscf.scf.rohf
   :members:

.. automodule:: pyscf.scf.uhf
   :members:

.. automodule:: pyscf.scf.uhf_symm
   :members:

.. automodule:: pyscf.scf.atom_hf
   :members:

.. automodule:: pyscf.scf.ghf
   :members:

.. automodule:: pyscf.scf.ghf_symm
   :members:

.. automodule:: pyscf.scf.cphf
   :members:

.. automodule:: pyscf.scf.ucphf
   :members:

Relativistic Hartree-Fock
-------------------------

.. automodule:: pyscf.scf.dhf
   :members:


Addons
------

.. automodule:: pyscf.scf.addons
   :members:

.. automodule:: pyscf.scf.chkfile
   :members:

.. automodule:: pyscf.scf.diis
   :members:

.. automodule:: pyscf.scf.jk
   :members:

Stability analysis
------------------

.. automodule:: pyscf.scf.stability
   :members:

.. automodule:: pyscf.scf.stability_slow
   :members:
