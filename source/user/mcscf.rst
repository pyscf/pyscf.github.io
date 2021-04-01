.. _user_mcscf:

Multi-configuration self-consistent field (MCSCF)
*************************************************

*Modules*: :mod:`mcscf`


Introduction
------------

Multiconfigurational self-consistent field (MCSCF) methods extend Hartree-Fock (HF) theory by describing the wave function as a linear combination of determinants. 
While there are a large variety of MCSCF methods, PySCF focuses on the complete active space (CAS) family of methods.
Unlike full configuration interaction (FCI), described in :numref:`theory_ci`, CAS methods perform the FCI procedure on a subset of the molecular orbitals, referred to as the "active space." 
These methods are crucial for systems that exhibit strong electron correlation, such as transition metal complexes.
For a detailed discussion of MCSCF methods, we direct the reader to References :cite:`Helgaker2013` and :cite:`esqc`.

The MCSCF module has two main flavors: CASCI and CASSCF. 
CASCI expands the wave function as a linear combination of slater determinants and variationally solves for the coefficients of this expansion.
CASSCF calculations optimize the orbitals for the set of determinants and the coefficients of their expansion.
This means performing multiple (and possibly many) CASCI step.

Minimal examples for each type of calculations are shown below.

CASCI
"""""
.. code-block:: python

  from pyscf import gto, mp, mcscf
  mol = gto.M(
      atom = 'O 0 0 0; O 0 0 1.2',
      basis = 'ccpvdz',
      spin = 2)
  myhf = mol.RHF().run()
  mymp = mp.UMP2(myhf).run()

  noons, natorbs = mcscf.addons.make_natural_orbitals(mymp)
  ncas, nelecas = (6,8)
  mycas = mcscf.CASCI(myhf, ncas, nelecas)
  mycas.kernel(natorbs)


CASSCF
""""""
.. code-block:: python

  import pyscf
  mol = pyscf.M(
      atom = 'O 0 0 0; O 0 0 1.2',
      basis = 'ccpvdz',
      spin = 2)
  myhf = mol.RHF().run()
  ncas, nelecas = (6,(5,3))
  mycas = myhf.CASSCF(ncas, nelecas).run()

Like many other modules in PySCF, :mod:`mcscf` works with density-fitting (:numref:`user_df`) and x2c (:numref:`user_x2c`).
It can also be used in its unrestricted formulation, please see :source:`examples/mcscf/60-uhf_based_ucasscf.py` for an example.
An important feature of :mod:`mcscf` is that it can interface with external CI solver such as DMRG, FCIQMC, or selected CI methods, see the external projects for more details.
You can even use MP2 or CC methods as solvers!
Check out :source:`/examples/cc/42-as_casci_fcisolver.py` for a working example.


Picking an Active Space
-----------------------
In general, selecting an active space can be cumbersome and PySCF offers severals ways to facilitate this process.
There are several strategies and they all contain two main components:

* Specify the number of electrons and orbitals in the active space.
* (Optional) Specify the which orbitals are "active" 

.. 
  warning::
  The total set orbitals (core, active, and virtual) used in active space methods can be specified or selected in a variety of ways, giving users substantial flexibility for CAS-type calculations.
  But users should note, "with great power comes great responsibility."
  Active space calculations are notoriously difficult and just because a calculation completes without error does not guarantee that the results will be chemically/physically meaningful so we urge users to select their active space orbitals with thought and care.

.. note::
  We always advise users to visualize their chosen active orbitals before starting large/expensive calculations.
  This involves dumping the MO coefficients to a ``molden`` file (see example :source:`examples/tools/02-molden.py`) and visualizing with your chosen program.
  While there are many great softwares available to visualize orbitals, `JMol <http://jmol.sourceforge.net/>`_ is one of the easiest to use and is recommended for less experienced users.


Below is a list of several general strategies one could employ to pick active space orbitals:

1) (Default) Specifying no additional information.
  This is the most minimal strategy for selecting an active space and chooses orbitals (and electrons) around the Fermi level that match the number of orbitals and electrons specified by the user.
  In most circumstances, this is not an ideal strategy and will lead to poor convergence or none at all.

  For example:
.. code-block:: python

  ncas, nelecas = (6,8)
  mycas = myhf.CASSCF(ncas, nelecas)


2) Specifying the molecular orbital (MO) index of the active space orbitals you want. 
  This is often useful after selecting (and typically visualizing) localized orbitals.
  The user can "manually" select the MO orbital indices (in a 1-based indexing scheme) and pass them to the ``sort_mo`` function.
  See :source:`examples/mcscf/10-define_cas_space.py` and :source:`examples/mcscf/34-init_guess_localization.py` for more details.

.. code-block:: python

  mycas = mcscf.CASSCF(myhf, 4, 4)
  # Note sort_mo by default take the 1-based orbital indices.
  mo = mycas.sort_mo([5,6,8,9])
  mycas.kernel(mo)


3) Specifying the number of orbitals in each symmetry group. 
This strategy can occasionally be helpful when the initial guess orbitals are not easily identifiable.

.. code-block:: python

  mycas = mcscf.CASSCF(mf, 12, 12)
  ncore = {'A1g':5, 'A1u':5}
  ncas = {'A1g':2, 'A1u':2,'E1ux':1, 'E1uy':1, 'E1gx':1, 'E1gy':1,
              'E2ux':1, 'E2uy':1, 'E2gx':1, 'E2gy':1}
  mo = mcscf.sort_mo_by_irrep(mycas, mf.mo_coeff, ncas, ncore)
  mycas.kernel(mo)

A similar approach where we specify the electron occupations by irreducible representation is also possible by setting ``mycas.fcisolver.irrep_nelec``.

.. code-block:: python

  mycas = mcscf.CASSCF(myhf, 8, 8)
  mycas.fcisolver.irrep_nelec = {"A1g": (2, 1), "A1u": (1, 1), "E1ux": (1, 1), "E1uy": (1, 0)}

4) Use automated strategies (``avas`` and ``dmet_cas``) to pick an active space based on AO orbitals you're targeting.
For more details, see :source:`examples/mcscf/43-avas.py` and :source:`examples/mcscf/43-dmet_cas.py`.

.. code-block:: python

  from pyscf.mcscf import avas
  ao_labels = ['Fe 3d', 'Fe 4d', 'C 2pz']
  ncas, nelecas, orbs = avas.avas(mf, ao_labels)
  mycas = mcscf.CASSCF(mf, ncas, nelecas)

  

.. code-block:: python

  from pyscf.mcscf import dmet_cas
  ao_labels = ['Fe 3d', 'Fe 4d', 'C 2pz']
  ncas, nelecas, mo = dmet_cas.guess_cas(mf, mf.make_rdm1(), ao_labels)
  mycas = mcscf.CASSCF(mf, ncas, nelecas)
  mycas.kernel(mo)


Frozen Core MCSCF
-----------------

To reduce to computational expense of CASSCF calculations, users can "freeze" orbitals thereby excluding them from optimization.

Users can specify the number of lowest orbitals to freeze:

.. code-block:: python

  mycas = mcscf.CASSCF(myhf, 6, 8)
  mycas.frozen = 2
  mycas.kernel()


Users can also specify a list of orbital indices (0-based).
These may be occupied, virtual, or active orbitals.

.. code-block:: python
  mycas = mcscf.CASSCF(myhf, 6, 8)
  mycas.frozen = [0,1,26,27]
  mycas.kernel()

See :source:`examples/mcscf/19-frozen_core.py` for a complete example.


State-Averaged Calculations
---------------------------

When dealing with states that are close in energy, it can be helpful to perform state average calculations where the orbitals are optimized for multiple states.
The ``state_average_`` function (note the hanging underscore) is a member function of ``CASCI``/``CASSCF`` objects and takes the weights of the states as input.
The weights can be any normalized and non-negative array of values, but typically they are all the same.
See Section 12.7.2 in Ref. :cite:`Helgaker2013` for more details.

.. code-block:: python

  n_states = 5
  weights = np.ones(n_states)/n_states
  mycas = mcscf.CASSCF(mf, 4, 4).state_average_(weights)

See :source:`examples/mcscf/15-state_average.py` for a complete example.


By default, only a single spin and/or point group symmetry is targeted, but it is possible to target a mixture of both:

.. code-block:: python

  weights = [.5, .5]
  solver1 = fci.direct_spin1_symm.FCI(mol)
  solver1.wfnsym= 'A1'
  solver1.spin = 0
  solver2 = fci.direct_spin1_symm.FCI(mol)
  solver2.wfnsym= 'A2'
  solver2.spin = 2

  mycas = mcscf.CASSCF(mf, 4, 4)
  mcscf.state_average_mix_(mycas, [solver1, solver2], weights)
  mycas.kernel()

See :source:`examples/mcscf/41-state_average.py` for a complete example.


Job Control
-----------

Optimization Settings
"""""""""""""""""""""

For CASSCF calculations, users may want to modify several of the convergence thresholds such as the energy (``conv_tol``), the orbital gradient (``conv_tol_grad``), and the maximum number of MCSCF iterations (``max_cycle_macro``).

.. code-block:: python

  mycas = mcscf.CASSCF(mf, 6, 6)
  mycas.conv_tol = 1e-12
  mycas.conv_tol_grad = 1e-6
  mycas.max_cycle_macro = 25
  mycas.kernel()


Initial Guess
"""""""""""""

Initial guess orbitals for the CASSCF calculation (starting orbitals) may be passed to the ``kernel`` member function of an MCSCF object.

.. code-block:: python

  mycas = mcscf.CASSCF(myhf, 8, 8)
  mycas.kernel(my_custom_mos)


CI coefficients from a previous calculation can also be passed as an initial guess to expedite the calculation:

.. code-block:: python
  mycas = mcscf.CASSCF(myhf, 8, 8)
  mycas.kernel(my_custom_mos, my_custom_ci)

Examples:

* :source:`examples/mcscf/14-project_init_guess.py`
* :source:`examples/mcscf/31-cr2_scan/cr2-scan.py`
* :source:`examples/mcscf/34-init_guess_localization.py`
* :source:`examples/mcscf/43-avas.py`
* :source:`examples/mcscf/43-dmet_cas.py`


Restarting
""""""""""

.. warning::
  When running large calculations, it's always recommended that you specify a checkpoint file for your calculation.

.. code-block:: python

  mycas.chkfile = "casscf.chk"

Much like :mod:`scf`, if a job is interrupted, users can restart the MCSCF calculations using checkpoint files from a crashed calculation.

.. code-block:: python

  from pyscf.lib import chkfile
  old_chk_file = "old_casscf.chk"
  mycas = mcscf.CASSCF(scf.RHF(mol), 6, 6)
  mycas.chkfile = "restarted_casscf.chk"
  mo = chkfile.load(old_chk_file, 'mcscf/mo_coeff')
  mycas.kernel(mo)


See :source:`examples/mcscf/13-restart.py` for a complete example.

Restarting calculations can be also be useful when using results from a smaller active space to speed up calculations on a larger one.


Observables and Properties
--------------------------

Wave Function Analysis
""""""""""""""""""""""

The ``analyze`` member functions of MCSCF objects prints many useful properties to ``stdout`` when the verbosity is >=4.

1) Natural orbital occupancies
2) Natural orbital AO expansions
3) Overlap between canonical MCSCF orbitals and the initial guess orbitals.
4) Analysis of the CI coefficients, i.e. the leading configurations and their weights
5) AO populations
6) Atomic populations
7) AO spin densities (if applicable)
8) Atomic spin densities (if applicable)

.. code-block:: python

  mycas = myhf.CASCI(6, 8).run()
  mycas.verbose = 4
  mycas.analyze()


Natural Orbitals
""""""""""""""""

Energy of CAS state is invariant under orbital rotation within inactive, active, and virtual sectors.
Inactive and virtual orbitals are by default canonicalized, i.e. they transformed such that Fock matrices within virtual and inactive sectors are diagonal, and orbital energy can be assigned to these orbitals.
By default active orbitals are kept untouched after orbital optimization and strictly speaking no energy or electron occupation can be assigned to them.
Users can request that active orbitals be transformed to the so-called natural representation, such that the one-body density matrix is diagonal and electron occupation can be assigned to them.

.. warning::
  When ``mycas.natorb`` is set, the natural orbitals may NOT be sorted by the active space occupancy.

.. code-block:: python

  mycas = mcscf.CASSCF(myhf, 6, 8)
  mycas.natorb = True
  mycas.kernel()



References
----------

.. bibliography:: ref_mcscf.bib
   :style: unsrt
