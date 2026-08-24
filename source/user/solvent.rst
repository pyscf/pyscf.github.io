
.. _user_solvent:

Solvation models
****************

*Modules*: :py:mod:`pyscf.solvent`

.. module:: solvent
   :synopsis: Solvation models and solvent effects
.. sectionauthor:: Qiming Sun <osirpt.sun@gmail.com> and Xiaojie Wu <wxj6000@gmail.com>.

Introduction
============

Solvation model allows the quantum chemistry calculations to include the
interactions between solvents and the quantum solute. Solvents can be treated
implicitly, known as continuum solvents, and explicitly. For continuum solvents,
we implemented the PCM (polarizable continuum model), ddCOSMO (domain-decomposition COSMO solvation model),
ddPCM (domain-decomposition polarizable continuum model),
and SMD (Solvation Model Density). For the explicit solvent environment,
we provided the interface to call the polarizable embedding library CPPE.

PCM model
=========
PySCF support four types of PCM solvent models, i.e. C-PCM, IEF-PCM, SS(V)PE, and COSMO
(See https://manual.q-chem.com/5.2/Ch12.S2.SS2.html for more detailed descriptions of these methods).
The analytical gradient and semi-analytical Hessian are also supported.

PCM solvent models can be applied on to an SCF object::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RKS(xc='b3lyp').PCM()
  mf.with_solvent.method = 'IEF-PCM' # C-PCM, SS(V)PE, COSMO
  mf.with_solvent.eps = 78.3553 # for water
  mf.run()

In regular MCSCF (CASCI or CASSCF), and post-SCF calculations, the setup for
self-consistent solvent is similar::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mc = mol.CASCI(4, 4).PCM().run()

  mp2_model = mol.MP2().PCM().run()

Solvent for excited states
--------------------------
When combining to TDDFT or other methods of excited states, solvent can be
modelled in the manner of self-consistency (fast solvent) or single shot (slow
solvent). Below we use TDDFT to demonstrate the treatments of fast solvent
and slow solvent.

In vertical excitations, the solvent almost does not respond to the change of
electronic structure. It should be viewed as slow solvent. The calculation
can be started with an SCF with fully relaxed solvent and followed by a regular
TDDFT method::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RHF().PCM().run()
  td = mf.TDA().run()

In the diabatic excitations, we would like to let the solvent rapidly responds
to the electronic structure of excited states. The entire solvation system
should converge to an equilibrium state between solvent and the excited state of
the solute. In this scenario, solvent model should be applied on to the excited
state methods::

  mf = mol.RHF().PCM().run()
  td = mf.TDA().PCM()
  td.with_solvent.equilibrium_solvation = True
  td.kernel()

Please note that the flag ``equilibrium_solvation`` needs to be set to ``True``
in this case. PySCF by default assumes the slow solvent model for TDDFT.

The slow solvent uses the optical dielectric constant ``eps_optical`` rather
than the static one. It is taken from the solvent database when the solvent is
specified by name, so prefer ``PCM('acetonitrile')`` over setting ``eps``
directly here, see `Solvent parameters`_.

In the complicated procedure which involves for example electronic states from
different states (typically in the MCSCF calculations with state-average or
state-specific approximations), PySCF PCM implementation allows to input a
density matrix and freeze the solvent equilibrated against the input density
matrix::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mc = mol.CASCI(4, 4).PCM()
  mc.fcisolver.nstates = 5
  mc.with_solvent.state_id = 1  # Slow solvent wrt the first excited state
  mc.run()

The slow solvent does not have to be corresponding to a particular state. It can
be even the solvent from a different geometry or an artificial quantum state of
solute::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  scf_dm = mol.RHF().run().make_rdm1()

  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    1.035   -1.082
       H  0.   -1.035   -1.082''',
                basis='6-31g*', verbose=4)
  mc = mol.CASCI(4, 4).PCM(dm=scf_dm).run()


Solvent parameters
------------------
The default solvent in the PCM module is water. Other solvents can be selected
by name::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RHF().PCM('methanol')
  mf.run()

The name can be given to the ``PCM`` constructor, to the ``.PCM()`` method of a
mean-field object, or assigned to the solvent object at any point::

  from pyscf.solvent import pcm
  cm = pcm.PCM(mol, 'methanol')
  cm.solvent = 'DMSO'

Names are matched case-insensitively and ignoring spaces and hyphens, so
``N,N-dimethylformamide``, ``N,N-DiMethylFormamide`` and ``dmf`` all select the
same solvent. The abbreviations ``h2o``, ``dmso``, ``dmf``, ``dma``, ``thf``,
``dcm``, ``ccl4``, ``chcl3``, ``mecn``, ``acn``, ``meoh``, ``etoh``, ``et2o``
and ``ether`` are recognized as well. The full list of solvents is the set of
keys of ``pyscf.solvent.smd.solvent_db``.

Assigning ``solvent`` sets two parameters:

* ``eps``, the static dielectric constant, which determines the ground state
  solvation;
* ``eps_optical``, the optical (high-frequency) dielectric constant, which is
  the square of the refractive index of the solvent. It is used only by the
  non-equilibrium solvation of excited states, see `Solvent for excited
  states`_.

Without a solvent name, ``eps_optical`` falls back to the value for water and a
warning is issued whenever ``eps`` indicates a different solvent. Specifying the
solvent by name is therefore the recommended way to set up a non-equilibrium
excited state calculation::

  mf = mol.RHF().PCM('tetrahydrofuran').run()
  td = mf.TDA()
  td.kernel()

Either parameter can still be assigned by hand, which overrides the value taken
from the database::

  mf = mol.RHF().PCM()
  mf.with_solvent.eps = 32.613   # methanol
  mf.run()

The solvent parameters are taken from the Minnesota Solvent Descriptor Database
(https://comp.chem.umn.edu/solvation/mnsddb.pdf), except for the dielectric
constant of water, which is the more accurate value 78.3553 used by
https://gaussian.com/scrf/. The SMD model reads its solvent names and its
descriptors from the same database.

ddCOSMO
=======

Self-consistent solvents for ground state
-----------------------------------------

Solvent model can be applied on to an SCF object::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RKS(xc='b3lyp').DDCOSMO().run()


SMD model
=========
SMD model is recommended for computing solvation free energy. The implementation of SMD model
in PySCF is based on IEF-PCM. Other SMx models are not supported yet
(See https://manual.q-chem.com/5.2/Ch12.S2.SS8.html). The source code for CDS contribution
is taken from NWChem (https://github.com/nwchemgit/nwchem/blob/master/src/solvation/mnsol.F)
SMD solvent solvent models can be applied on to an SCF object::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RKS(xc='b3lyp').SMD()
  mf.with_solvent.solvent = 'water'
  mf.run()

The solvent descriptor can be assigned to the ``with_solvent`` directly::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RKS(xc='b3lyp').SMD()
  mf.with_solvent.solvent_descriptors = [1.3843, 1.3766, 0.0, 0.45, 35.06, 13.45, 0.0, 0.0]
  mf.run()

The format of solvant names are the same as Minnesota Solvent Descriptor Database
(https://comp.chem.umn.edu/solvation/mnsddb.pdf). One can assign solvent descriptors in the format
``mf.with_solvent.solvent_descriptors = [n, n25, alpha, beta, gamma, epsilon, phi, psi]``

Polarizable embedding
=====================
To use polarizable embedding model for mean-field calculations, one would need
to first generate potential data for the input of CPPE library. The best way to
generate potential files is with `PyFraME <https://gitlab.com/FraME-projects/PyFraME>`_.
You can directly throw in a pdb file, select the QM region and how to parametrize
different parts of the environment (with either pre-defined potentials, or with LoProp).
Some guidance is also provided in the Tutorial Review paper about PE, section 4:
https://onlinelibrary.wiley.com/doi/full/10.1002/qua.25717
Therein, the format of the potential file is also explained (it’s the same
format as used in the original Dalton pelib implementation).

With the generated potential file, one can carry out the polarizable embedding
calculations::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RKS(xc='b3lyp')
  mf = pyscf.solvent.PE(mf, 'potfile')
  mf.run()
