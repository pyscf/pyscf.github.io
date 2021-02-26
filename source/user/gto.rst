.. _theory_gto:

Molecules and crystals
**********************

*Modules*: :mod:`gto`, :mod:`pbc.gto`

Initializing a molecule
=======================

There are three ways to define and initialize a molecule.  The first is to use
the keyword arguments of :func:`Mole.build` to initialize a molecule::

  >>> from pyscf import gto
  >>> mol = gto.Mole()
  >>> mol.build(
  ...     atom = '''O 0 0 0; H  0 1 0; H 0 0 1''',
  ...     basis = 'sto-3g')

The second way is to assign the geometry, basis etc. to :class:`Mole`
object, then call :meth:`~Mole.build` function to initialize the
molecule::

  >>> mol = gto.Mole()
  >>> mol.atom = '''O 0 0 0; H  0 1 0; H 0 0 1'''
  >>> mol.basis = 'sto-3g'
  >>> mol.build()

The third way is to use the shortcut function :func:`Mole.M`.  This
function pass all arguments to :func:`Mole.build`::

  >>> from pyscf import gto
  >>> mol = gto.M(
  ...     atom = '''O 0 0 0; H  0 1 0; H 0 0 1''',
  ...     basis = 'sto-3g')

Either way, you may have noticed two keywords ``atom`` and ``basis``.
They are used to hold the molecular geometry and basis sets.

For more examples, see :ref:`gto`. 

You can also impose symmetry::

  >>> o2_sym = mol.copy()
  >>> o2_sym.spin = 2
  >>> o2_sym.symmetry = 1
  >>> o2_sym.build(0, 0)
  >>> rhf3_sym = scf.RHF(o2_sym)
  >>> print(rhf3_sym.kernel())
  -149.609461122

Here we rebuild the molecule because we need to initialize the point group
symmetry information, symmetry adapted orbitals.  We can check the occupancy for
each irreducible representations::

  >>> import numpy
  >>> from pyscf import symm
  >>> def myocc(mf):
  ...     mol = mf.mol
  ...     irrep_id = mol.irrep_id
  ...     so = mol.symm_orb
  ...     orbsym = symm.label_orb_symm(mol, irrep_id, so, mf.mo_coeff)
  ...     doccsym = numpy.array(orbsym)[mf.mo_occ==2]
  ...     soccsym = numpy.array(orbsym)[mf.mo_occ==1]
  ...     for ir,irname in enumerate(mol.irrep_name):
  ...         print('%s, double-occ = %d, single-occ = %d' %
  ...               (irname, sum(doccsym==ir), sum(soccsym==ir)))
  >>> myocc(rhf3_sym)
  Ag, double-occ = 3, single-occ = 0
  B1g, double-occ = 0, single-occ = 0
  B2g, double-occ = 0, single-occ = 1
  B3g, double-occ = 0, single-occ = 1
  Au, double-occ = 0, single-occ = 0
  B1u, double-occ = 2, single-occ = 0
  B2u, double-occ = 1, single-occ = 0
  B3u, double-occ = 1, single-occ = 0

To label the irreducible representation of given orbitals,
:func:`symm.label_orb_symm` needs the information of the point group
symmetry which are initialized in ``mol`` object, including the `id` of
irreducible representations :attr:`Mole.irrep_id` and the symmetry
adapted basis :attr:`Mole.symm_orb`.  For each :attr:`~Mole.irrep_id`,
:attr:`Mole.irrep_name` gives the associated irrep symbol (A1, B1 ...).
In the SCF calculation, you can control the symmetry of the wave
function by assigning the number of alpha electrons and beta electrons
`(alpha,beta)` for some irreps::

  >>> rhf3_sym.irrep_nelec = {'B2g': (1,1), 'B3g': (1,1), 'B2u': (1,0), 'B3u': (1,0)}
  >>> rhf3_sym.kernel()
  >>> print(rhf3_sym.kernel())
  -148.983117701
  >>> rhf3_sym.get_irrep_nelec()
  {'Ag' : (3, 3), 'B1g': (0, 0), 'B2g': (1, 1), 'B3g': (1, 1), 'Au' : (0, 0), 'B1u': (1, 0), 'B2u': (0, 1), 'B3u': (1, 0)}

More informations of the calculation can be found in the output file ``o2.log``.

Geometry
--------

The molecular geometry can be input in Cartesian format::

  >>> mol = gto.Mole()
  >>> mol.atom = '''O 0, 0, 0
  ... H   0  1  0; H 0, 0, 1'''

The atoms in the molecule are represented by an element symbol plus
three numbers for coordinates.  Different atoms should be separated by
``;`` or a line break. In the same atom, ``,`` can be used to separate
different items.  The input parser also supports the Z-matrix input
format::

  >>> mol = gto.Mole()
  >>> mol.atom = '''O
  ... H, 1, 1.2;  H   1 1.2   2 105'''

Similarly, different atoms need to be separated by ``;`` or a line
break.  If you need to label an atom to distinguish it from the
others, you can prefix or suffix the atom symbol with a number
``1234567890`` or a special character ``~!@#$%^&*()_+.?:<>[]{}|`` (not
``,`` and ``;``). With this decoration, you can specify different
basis sets, masses, or nuclear models on different atoms::

  >>> mol = gto.Mole()
  >>> mol.atom = '''8 0 0 0; h:1 0 1 0; H@2 0 0'''
  >>> mol.basis = {'O': 'sto-3g', 'H': 'cc-pvdz', 'H@2': '6-31G'}
  >>> mol.build()
  >>> print(mol._atom)
  [['O', [0.0, 0.0, 0.0]], ['H:1', [0.0, 1.0, 0.0]], ['H@2', [0.0, 0.0]]]

Basis set
---------

The simplest way to define the basis set is to assign the name of the
basis as a string to :attr:`mol.basis`::

  mol.basis = 'sto3g'

This input will apply the specified basis set to all atoms. The name
of the basis set in the string is case insensitive.  White spaces,
dashes and underscores in the name are all ignored.  If different
basis sets are required for different elements, a Python ``dict`` can
be used::

  mol.basis = {'O': 'sto3g', 'H': '6-31g'}

You can find more examples in section :ref:`input_basis` and in the file
:file:`examples/gto/04-input_basis.py`.

Other parameters
----------------

You can assign more information to the molecular object::

  mol.symmetry = 1
  mol.charge = 1
  mol.spin = 1
  mol.nucmod = {'O1': 1} 
  mol.mass = {'O1': 18, 'H': 2} 

.. note:: :attr:`Mole.spin` is *2S*, the number of unpaired electrons
  i.e. the difference between the number of alpha and beta electrons.

:class:`Mole` also defines some global parameters.  You can control the
print level globally with :attr:`~Mole.verbose`::

  mol.verbose = 4

The print level can be 0 (quiet, no output) to 9 (very noisy).  The
most useful messages are printed at level 4 (info) and 5 (debug).  You
can also specify a place where to write the output messages::

  mol.output = 'path/to/my_log.txt'

If this variable is not assigned, messages will be dumped to
:attr:`sys.stdout`.

The maximum memory usage can be controlled globally::

  mol.max_memory = 1000 # MB
  
The default size can also be defined with the shell environment
variable `PYSCF_MAX_MEMORY`

:attr:`~Mole.output` and :attr:`~Mole.max_memory` can be assigned from command
line::

  $ python example.py -o /path/to/my_log.txt -m 1000


Initializing a crystal
======================

Initialization a crystal unit cell is very similar to the initialization
molecular object.  Here, :class:`pyscf.pbc.gto.Cell` class should be used
instead of the :class:`pyscf.gto.Mole` class::

  >>> from pyscf.pbc import gto
  >>> cell = gto.Cell()
  >>> cell.atom = '''H  0 0 0; H 1 1 1'''
  >>> cell.basis = 'gth-dzvp'
  >>> cell.pseudo = 'gth-pade'
  >>> cell.a = numpy.eye(3) * 2
  >>> cell.build()

The crystal initialization requires an extra parameter :attr:`cell.a` which
represents the lattice vectors. In the above example, we specified
:attr:`cell.pseudo` for the pseudo-potential of the system which is an optional
parameter.  The input format of basis set is the same to that of :class:`Mole`
object.  The other attributes of :class:`Mole` object such as :attr:`verbose`,
:attr:`max_memory`, :attr:`spin` can also be used in the crystal systems.
More details of the crystal :class:`Cell` object and the relevant input
parameters are documented in :ref:`pbc_gto`.

1D and 2D systems
-----------------

PySCF PBC module supports the low-dimensional PBC systems.  You can initialize
the attribute :attr:`cell.dimension` to specify the dimension of the system::

  >>> from pyscf.pbc import gto
  >>> cell = gto.Cell()
  >>> cell.atom = '''H  0 0 0; H 1 1 0'''
  >>> cell.basis = 'sto3g'
  >>> cell.dimension = 2
  >>> cell.a = numpy.eye(3) * 2
  >>> cell.build()

When :attr:`cell.dimension` is specified, a vacuum of infinite size will be
applied on certain dimension(s).  More specifically, when :attr:`cell.dimension`
is 2, the z-direction will be treated as infinite large and the xy-plane
constitutes the periodic surface. When :attr:`cell.dimension` is 1, y and z axes
are treated as vacuum thus wire is placed on the x axis.  When
:attr:`cell.dimension` is 0, all three directions are vacuum.  The PBC system is
actually the same to the molecular system.



gto.Mole
========

Molecular Structure
-------------------

Basis (GTO)
-----------

ECP
---

Spin and charge
---------------

Symmetry
--------


pbc.gto.Cell
============

Unit cell
---------

Basis (crystalline GTO)
-----------------------
PySCF uses the crystalline Gaussian-type atomic orbitals (AOs) as the one-particle basis
for solid calculations with periodic boundary conditions (PBCs).
These AOs are Bloch waves with the following form:

.. math::

   \phi_{\mu\mathbf{k}}(\mathbf{r}) = \sum_{\mathbf{T}} e^{i\mathbf{k}\cdot \mathbf{T}} \chi_{\mu}(\mathbf{r}-\mathbf{T}) \;,

where :math:`\mathbf{T}` is a lattice translation vector,
:math:`\mathbf{k}` is a crystal momentum vector in the first Brillouin zone,
and :math:`\chi_{\mu}` is a Gaussian type orbital (GTO).

The basis set used in a PBC calculation can be specified exactly the same way
as that in the molecular case.
The following shows a simple example:
:source:`examples/pbc/04-input_basis.py`.

Pseudopotential
---------------



Access AO integrals
===================

molecular integrals
-------------------

PySCF uses `Libcint <https://github.com/sunqm/libcint>`_ library as the AO
integral engine.  It provides simple interface function :func:`getints_by_shell`
to evaluate integrals.  The following example evaluates 3-center 2-electron
integrals with this function::

  import numpy
  from pyscf import gto, scf, df
  mol = gto.M(atom='O 0 0 0; h 0 -0.757 0.587; h 0 0.757 0.587', basis='cc-pvdz')
  auxmol = gto.M(atom='O 0 0 0; h 0 -0.757 0.587; h 0 0.757 0.587', basis='weigend')
  pmol = mol + auxmol
  nao = mol.nao_nr()
  naux = auxmol.nao_nr()
  eri3c = numpy.empty((nao,nao,naux))
  pi = 0
  for i in range(mol.nbas):
      pj = 0
      for j in range(mol.nbas):
          pk = 0
          for k in range(mol.nbas, mol.nbas+auxmol.nbas):
              shls = (i, j, k)
              buf = pmol.intor_by_shell('int3c2e_sph', shls)
              di, dj, dk = buf.shape
              eri3c[pi:pi+di,pj:pj+dj,pk:pk+dk] = buf
              pk += dk
          pj += dj
      pi += di

Here we load the Weigend density fitting basis to ``auxmol`` and append the
basis to normal orbital basis which was initialized in ``mol``.  In the result
``pmol`` object, the first ``mol.nbas`` shells are the orbital basis and
the next ``auxmol.nbas`` are auxiliary basis.  The three nested loops run over
all integrals for the three index integral `(ij|K)`.  Similarly, we can compute
the two center Coulomb integrals::

  eri2c = numpy.empty((naux,naux))
  pk = 0
  for k in range(mol.nbas, mol.nbas+auxmol.nbas):
      pl = 0
      for l in range(mol.nbas, mol.nbas+auxmol.nbas):
          shls = (k, l)
          buf = pmol.intor_by_shell('int2c2e_sph', shls)
          dk, dl = buf.shape
          eri2c[pk:pk+dk,pl:pl+dl] = buf
          pl += dl
      pk += dk

Now we can use the two-center integrals and three-center integrals to implement
the density fitting Hartree-Fock code.

.. code:: python

  def get_vhf(mol, dm, *args, **kwargs):
      naux = eri2c.shape[0]
      nao = mol.nao_nr()
      rho = numpy.einsum('ijp,ij->p', eri3c, dm)
      rho = numpy.linalg.solve(eri2c, rho)
      jmat = numpy.einsum('p,ijp->ij', rho, eri3c)
      kpj = numpy.einsum('ijp,jk->ikp', eri3c, dm)
      pik = numpy.linalg.solve(eri2c, kpj.reshape(-1,naux).T)
      kmat = numpy.einsum('pik,kjp->ij', pik.reshape(naux,nao,nao), eri3c)
      return jmat - kmat * .5
      
  mf = scf.RHF(mol)
  mf.verbose = 0
  mf.get_veff = get_vhf
  print('E(DF-HF) = %.12f, ref = %.12f' % (mf.kernel(), scf.density_fit(mf).kernel()))

Your screen should output

  | E(DF-HF) = -76.025936299702, ref = -76.025936299702


Evaluating the integrals with nested loops and :func:`mol.intor_by_shell` method is
inefficient.  It is preferred to load integrals in bulk and this can be done
with :func:`mol.intor` method::

  eri2c = auxmol.intor('int2c2e_sph')
  eri3c = pmol.intor('int3c2e_sph', shls_slice=(0,mol.nbas,0,mol.nbas,mol.nbas,mol.nbas+auxmol.nbas))
  eri3c = eri3c.reshape(mol.nao_nr(), mol.nao_nr(), -1)

:func:`mol.intor` method can be used to evaluate one-electron integrals,
two-electron integrals::

  hcore = mol.intor('int1e_nuc_sph') + mol.intor('int1e_kin_sph')
  overlap = mol.intor('int1e_ovlp_sph')
  eri = mol.intor('int2e_sph')

There is a long list of supported AO integrals.  See :ref:`gto_moleintor`.


PBC AO integrals
----------------

:func:`mol.intor` can only be used to evaluate the integrals with open boundary
conditions.  When the periodic boundary conditions of crystal systems are
studied, you need to use :func:`pbc.Cell.pbc_intor` function to evaluate the
integrals of short-range operators, such as the overlap, kinetic matrix::

  from pyscf.pbc import gto
  cell = gto.Cell()
  cell.atom = 'H 0 0 0; H 1 1 1'
  cell.a = numpy.eye(3) * 2.
  cell.build()
  overlap = cell.pbc_intor('int1e_ovlp_sph')

By default, :func:`pbc.Cell.pbc_intor` function returns the :math:`\Gamma`-point
integrals.  If k-points are specified, function :func:`pbc.Cell.pbc_intor` can
also evaluate the k-point integrals::

  kpts = cell.make_kpts([2,2,2])  # 8 k-points
  overlap = cell.pbc_intor('int1e_ovlp_sph', kpts=kpts)

.. note:: :func:`pbc.Cell.pbc_intor` can only be used to evaluate the short-range
  integrals.  PBC density fitting method has to be used to compute the
  long-range operator such as nuclear attraction integrals, Coulomb integrals.

The two-electron Coulomb integrals can be evaluated with PBC density fitting
methods::

    from pyscf.pbc import df
    eri = df.DF(cell).get_eri()

See also :ref:`pbc_df` for more details of the PBC density fitting module.


Other features
==============
Density fitting
---------------

.. literalinclude:: /../examples/scf/20-density_fitting.py

Customizing Hamiltonian
-----------------------

.. literalinclude:: /../examples/scf/40-customizing_hamiltonian.py
