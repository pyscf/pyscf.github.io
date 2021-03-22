.. _user_gto:

Molecular structure
*******************

*Modules*: :mod:`gto`

Initializing a molecule
=======================

There are three ways to define and initialize a molecule. The first is to use
the keyword arguments of the :meth:`Mole.build` method to initialize a molecule::

  >>> from pyscf import gto
  >>> mol = gto.Mole()
  >>> mol.build(
  ...     atom = '''O 0 0 0; H  0 1 0; H 0 0 1''',
  ...     basis = 'sto-3g')

The second way is to assign the geometry, basis etc., to the :class:`Mole`
object, followed by calling the :meth:`~Mole.build` method::

  >>> from pyscf import gto
  >>> mol = gto.Mole()
  >>> mol.atom = '''O 0 0 0; H  0 1 0; H 0 0 1'''
  >>> mol.basis = 'sto-3g'
  >>> mol.build()

The third way is to use the shortcut functions :func:`pyscf.M` or :func:`Mole.M`.
These functions pass all the arguments to the :func:`~Mole.build` method::

  >>> import pyscf
  >>> mol = pyscf.M(
  ...     atom = '''O 0 0 0; H  0 1 0; H 0 0 1''',
  ...     basis = 'sto-3g')

  >>> from pyscf import gto
  >>> mol = gto.M(
  ...     atom = '''O 0 0 0; H  0 1 0; H 0 0 1''',
  ...     basis = 'sto-3g')

In any of these, you may have noticed two keywords ``atom`` and ``basis``.
They are used to hold the molecular `geometry`_ and `basis sets`_, which 
can be defined along with other input options as follows.

.. _geometry:

Geometry
--------

The molecular geometry can be input in Cartesian format 
with the default unit being Angstrom 
(one can specify the unit by setting the attribute :attr:`unit`
to either ``'Angstrom'`` or ``'Bohr'``)::

    mol = gto.Mole()
    mol.atom = '''
        O   0. 0. 0.
        H   0. 1. 0.
        H   0. 0. 1.
    '''
    mol.unit = 'B' # case insensitive, any string not starts by 'B' or 'AU' is treated as 'Angstrom'

The atoms in the molecule are represented by an element symbol plus
three numbers for coordinates.  Different atoms should be separated by
``;`` or a line break. In the same atom, ``,`` can be used to separate
different items. Blank lines or lines started with ``#`` will be
ignored::

  >>> mol = pyscf.M(
  ... atom = '''
  ... #O 0 0 0
  ... H 0 1 0
  ...
  ... H 0 0 1
  ... ''')
  >>> mol.natm
  2

The input parser also supports the Z-matrix input
format::

    mol = gto.Mole()
    mol.atom = '''
        O
        H  1  1.2  
        H  1  1.2  2 105
    '''

Similarly, different atoms need to be separated by ``;`` or a line
break.

The geometry string is case-insensitive. It also supports to input the nuclear
charges of elements::

  >>> mol = gto.Mole()
  >>> mol.atom = '''8 0. 0. 0.; h 0. 1. 0; H 0. 0. 1.'''

If you need to label an atom to distinguish it from the
others, you can prefix or suffix the atom symbol with a number
``1234567890`` or a special character ``~!@#$%^&*()_+.?:<>[]{}|`` (not
``,`` and ``;``). With this decoration, you can specify different
basis sets, masses, or nuclear models on different atoms::

  >>> mol = gto.Mole()
  >>> mol.atom = '''8 0 0 0; h:1 0 1 0; H@2 0 0 1'''
  >>> mol.unit = 'B'
  >>> mol.basis = {'O': 'sto-3g', 'H': 'cc-pvdz', 'H@2': '6-31G'}
  >>> mol.build()
  >>> print(mol._atom)
  [['O', [0.0, 0.0, 0.0]], ['H:1', [0.0, 1.0, 0.0]], ['H@2', [0.0, 0.0, 1.0]]]

You can also input the geometry in the internal format of
:attr:`Mole.atom`::

  atom = [[atom1, (x, y, z)],
          [atom2, (x, y, z)],
          ...
          [atomN, (x, y, z)]]

This is convenient as you can use Python script to construct the geometry::

  >>> mol = gto.Mole()
  >>> mol.atom = [['O',(0, 0, 0)], ['H',(0, 1, 0)], ['H',(0, 0, 1)]]
  >>> mol.atom.extend([['H', (i, i, i)] for i in range(1,5)])

Besides Python list, tuple and numpy.ndarray are all supported by the internal
format::

  >>> mol.atom = (('O',numpy.zeros(3)), ['H', 0, 1, 0], ['H',[0, 0, 1]])

No matter which format or symbols were used in the input, :func:`Mole.build`
will convert :attr:`Mole.atom` to the internal format::

  >>> mol.atom = '''
      O        0,   0, 0             ; 1 0.0 1 0

          H@2,0 0 1
      '''
  >>> mol.unit = 'B'
  >>> mol.build()
  >>> print(mol._atom)
  [('O', [0.0, 0.0, 0.0]), ('H', [0.0, 1.0, 0.0]), ('H@2', [0.0, 0.0, 1.0])]

which is stored as the attribute :attr:`Mole._atom`.

Once the :class:`Mole` object is built, the molecular geometry can be
accessed through the :meth:`Mole.atom_coords` function.  
This function returns a (N,3) array for the coordinates of each atom::

  >>> print(mol.atom_coords(unit='Bohr')) # unit can be "ANG" or "Bohr"
  [[0. 0. 0.]
   [0. 1. 0.]
   [0. 0. 1.]]

.. _basis sets:

Ghost atoms can also be specified when inputing the geometry.
See :source:`examples/gto/03-ghost_atom.py` for examples.

Basis set
---------

The simplest way to define the basis set is to assign the name of the
basis as a string to :attr:`Mole.basis`::

  mol.basis = 'sto-3g'

This input will apply the specified basis set to all atoms. The name
of the basis set in the string is case insensitive. White spaces,
dashes and underscores in the name are all ignored. If different
basis sets are required for different elements, a Python ``dict`` can
be used::

  mol.basis = {'O': 'sto-3g', 'H': '6-31g'}

One can also input custom basis sets with the helper functions.
The function :func:`gto.basis.parse` can parse a basis string in the NWChem format
(https://bse.pnl.gov/bse/portal)::

  mol.basis = {'O': gto.basis.parse('''
  C    S
       71.6168370              0.15432897
       13.0450960              0.53532814
        3.5305122              0.44463454
  C    SP
        2.9412494             -0.09996723             0.15591627
        0.6834831              0.39951283             0.60768372
        0.2222899              0.70011547             0.39195739
  ''')}

The functions :func:`gto.basis.load` can load arbitrary basis sets from the database, even
if the basis set does not match the element::

  mol.basis = {'H': gto.basis.load('sto3g', 'C')}

Both :func:`gto.basis.parse` and :func:`gto.basis.load` return the basis set in the
internal format (see :ref:`gto_basis`).

The basis parser also supports ghost atoms::

  mol.basis = {'GHOST': gto.basis.load('cc-pvdz', 'O'), 'H': 'sto3g'}

More examples of ghost atoms can be found in
:source:`examples/gto/03-ghost_atom.py`.

Like the requirements of geometry input, you can use atomic symbols
(case-insensitive) or atomic nuclear charges as the keyword of the
:attr:`~Mole.basis` dictionary. Prefix and suffix of numbers and special
characters are allowed. If the decorated atomic symbol is appeared in
:attr:`~Mole.atom` but not in :attr:`~Mole.basis`, the basis parser will
remove all decorations and seek the pure atomic symbol in the
:attr:`~Mole.basis` dictionary. In the following example, the ``6-31G`` basis set
will be assigned to the atom ``H1``,
but the ``STO-3G`` basis will be used for the atom ``H2``:

  mol.atom = '8 0 0 0; h1 0 1 0; H2 0 0 1'
  mol.basis = {'O': 'sto-3g', 'H': 'sto3g', 'H1': '6-31G'}

See :source:`examples/gto/04-input_basis.py` for more examples.

ECP
---
Effective core potentials (ECPs) can be specified with the attribute :attr:`Mole.ecp`. 
Scalar type ECPs are available for all molecular and crystal methods. 
The built-in scalar ECP datasets include

============ ========================
Keyword      Comment
------------ ------------------------
bfd
cc-pvdz-pp
cc-pvtz-pp   same to cc-pvdz-pp
cc-pvqz-pp   same to cc-pvdz-pp
cc-pv5z-pp   same to cc-pvdz-pp
crenbl
crenbs
def2-svp
def2-svpd    same to def2-svp
def2-tzvp    same to def2-svp
def2-tzvpd   same to def2-svp
def2-tzvpp   same to def2-svp
def2-tzvppd  same to def2-svp
def2-qzvp    same to def2-svp
def2-qzvpd   same to def2-svp
def2-qzvpp   same to def2-svp
def2-qzvppd  same to def2-svp
lanl2dz
lanl2tz
lanl08
sbkjc
stuttgart
============ ========================

ECP parameters can be specified directly in input script using NWChem format.
Examples of ECP input can be found in :source:`examples/gto/05-input_ecp.py`.

Spin-orbit (SO) ECP integrals can be evaluated using PySCF's integral driver. 
However, SO-ECPs are not automatically applied to any methods in the current implementation. 
They need to be added to the core Hamiltonian as shown in the examples
:source:`examples/gto/20-soc_ecp.py` and
:source:`examples/scf/44-soc_ecp.py`.
PySCF provides the following SO-ECPs

============ ========================
Keyword      Comment
------------ ------------------------
crenbl
crenbs
============ ========================

.. note::
  Be careful with the SO-ECP conventions when inputing them directly in the
  input script. SO-ECP parameters may take different conventions in different
  packages. More particularly, the treatment of the factor 2/(2l+1).  PySCF
  assumes that this factor was multiplied in the SOC parameters. See also
  relevant discussions in `Dirac doc
  <http://www.diracprogram.org/doc/master/molecule_and_basis/molecule_with_ecp.html>`_
  and `NWChem doc <https://nwchemgit.github.io/ECP.html>`_.

Point group symmetry
--------------------
You can also invoke point group symmetry for molecular calculations
by setting the attribute :attr:`Mole.symmetry` to ``True``::

    >>> mol = pyscf.M(
    ...     atom = 'B 0 0 0; H 0 1 1; H 1 0 1; H 1 1 0',
    ...     symmetry = True
    ... )

The point group symmetry information is held in the :class:`Mole` object.  
The symmetry module (:mod:`symm`) of PySCF can detect arbitrary point groups. 
The detected point group is saved in :attr:`Mole.topgroup`,
and the supported subgroup is saved in :attr:`Mole.groupname`::

    >>> print(mol.topgroup)
    C3v
    >>> print(mol.groupname)
    Cs

Currently, PySCF supports the linear molecular symmetries
:math:`D_{\infty h}` (labelled as ``Dooh`` in the program) and :math:`C_{\infty v}`
(labelled as ``Coov``), the :math:`D_{2h}` group and its subgroups.
Sometimes it is necessary to use a lower symmetry instead of the detected
symmetry group. The subgroup symmetry can be specified in
by :attr:`Mole.symmetry_subgroup` and the program will first detect the highest
possible symmetry group and then lower the point group symmetry to the specified
subgroup::

    >>> mol = gto.Mole()
    >>> mol.atom = 'N 0 0 0; N 0 0 1'
    >>> mol.symmetry = True
    >>> mol.symmetry_subgroup = C2
    >>> mol.build()
    >>> print(mol.topgroup)
    Dooh
    >>> print(mol.groupname)
    C2

When a particular symmetry is assigned to :attr:`Mole.symmetry`,
the initialization function :meth:`Mole.build` will test whether the molecule
geometry is subject to the required symmetry. If not, initialization will stop
and an error message will be issued::

    >>> mol = gto.Mole()
    >>> mol.atom = 'O 0 0 0; C 0 0 1'
    >>> mol.symmetry = 'Dooh'
    >>> mol.build()
    RuntimeWarning: Unable to identify input symmetry Dooh.
    Try symmetry="Coov" with geometry (unit="Bohr")
    ('O', [0.0, 0.0, -0.809882624813598])
    ('C', [0.0, 0.0, 1.0798434997514639])

.. note::
  :attr:`Mole.symmetry_subgroup` has no effects
  when specific symmetry group is assigned to :attr:`Mole.symmetry`.

When symmetry is enabled in the :class:`Mole` object, the point group symmetry
information will be used to construct the symmetry adapted orbital basis (see
also :mod:`symm`).  The symmetry adapted orbitals are held in
:attr:`Mole.symm_orb` as a list of 2D arrays.  Each element of the list
is an AO (atomic orbital) to SO (symmetry-adapted orbital) transformation matrix
of an irreducible representation. The name of the irreducible representations
are stored in :attr:`Mole.irrep_name` and their internal IDs (see more details
in :mod:`symm`) are stored in :attr:`Mole.irrep_id`::

  >>> mol = gto.Mole()
  >>> mol.atom = 'O 0 0 0; O 0 0 1.2'
  >>> mol.spin = 2
  >>> mol.symmetry = "D2h"
  >>> mol.build()
  >>> for s,i,c in zip(mol.irrep_name, mol.irrep_id, mol.symm_orb):
  ...     print(s, i, c.shape)
  Ag 0 (10, 3)
  B2g 2 (10, 1)
  B3g 3 (10, 1)
  B1u 5 (10, 3)
  B2u 6 (10, 1)
  B3u 7 (10, 1)

These symmetry-adapted orbitals are used as basis functions for the 
following SCF or post-SCF calculations::

  >>> mf = scf.RHF(mol)
  >>> mf.kernel()
  converged SCF energy = -147.631655286561

and we can check the occupancy of the MOs in each irreducible representation::

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
  >>> myocc(mf)
  Ag, double-occ = 3, single-occ = 0
  B2g, double-occ = 0, single-occ = 0
  B3g, double-occ = 0, single-occ = 1
  B1u, double-occ = 0, single-occ = 1
  B2u, double-occ = 0, single-occ = 0
  B3u, double-occ = 2, single-occ = 0

To label the irreducible representation of given orbitals,
:func:`symm.label_orb_symm` needs the information of the point group
symmetry which are initialized in the :class:`Mole` object, including the IDs of
irreducible representations (:attr:`Mole.irrep_id`) and the symmetry
adapted basis :attr:`Mole.symm_orb`. For each :attr:`~Mole.irrep_id`,
:attr:`Mole.irrep_name` gives the associated irrep symbol (A1, B1 ...).
In the SCF calculation, you can control the symmetry of the wave
function by assigning the number of alpha electrons and beta electrons
`(alpha,beta)` for the irreps::

  >>> mf.irrep_nelec = {'B2g': (1,1), 'B3g': (1,1), 'B2u': (1,0), 'B3u': (1,0)}
  >>> mf.kernel()
  converged SCF energy = -146.97333043702
  >>> mf.get_irrep_nelec()
  {'Ag': (3, 3), 'B2g': (1, 1), 'B3g': (1, 1), 'B1u': (2, 2), 'B2u': (1, 0), 'B3u': (1, 0)}

Spin and charge
---------------

Charge and spin multiplicity can be assigned to :class:`Mole` object::

  mol.charge = 1
  mol.spin = 1

.. note:: 
  :attr:`Mole.spin` is the number of unpaired electrons *2S*,
  i.e. the difference between the number of alpha and beta electrons.

These two attributes do not affect any other parameters 
in the :attr:`Mole.build` initialization function.
They can be set or modified after the
:class:`Mole` object is initialized::

  >>> mol = gto.Mole()
  >>> mol.atom = 'O 0 0 0; h 0 1 0; h 0 0 1'
  >>> mol.basis = 'sto-6g'
  >>> mol.spin = 2
  >>> mol.build()
  >>> print(mol.nelec)
  (6, 4)
  >>> mol.spin = 0
  >>> print(mol.nelec)
  (5, 5)

The attribute :attr:`Mole.charge` is the parameter to define the total number of electrons in the
system. For custom systems such as the Hubbard lattice model, the total number
of electrons needs to be specified directly by setting the attribute :attr:`Mole.nelectron`::

  >>> mol = gto.Mole()
  >>> mol.nelectron = 10

Other parameters
----------------

You can assign more information to the molecular object::

  mol.nucmod = {'O1': 1} # nuclear charge model: 0-point charge, 1-Gaussian distribution
  mol.mass = {'O1': 18, 'H': 2}  # atomic mass

See :source:`examples/gto/07-nucmod.py` for more examples of nuclear charge models.

The :class:`Mole` class also defines some global parameters.  You can control the
print level globally with :attr:`~Mole.verbose`::

  mol.verbose = 4

The print level can be 0 (quiet, no output) to 9 (very noisy).  The
most useful messages are printed at level 4 (info) and 5 (debug).  You
can also specify a place where to write the output messages::

  mol.output = 'path/to/log.txt'

If this variable is not assigned, messages will be dumped to
:attr:`sys.stdout`.

The maximum memory usage can be controlled globally::

  mol.max_memory = 1000 # MB
  
The default size can also be defined with the shell environment
variable `PYSCF_MAX_MEMORY`.

The attributes :attr:`~Mole.output` and :attr:`~Mole.max_memory` can also be 
assigned from command line::

  $ python input.py -o /path/to/my_log.txt -m 1000

By default, command line has the highest priority, which means the
settings in the script will be overwritten by the command line arguments.  
To make the input parser ignore the command line arguments, you can call the
:func:`Mole.build` with::

  mol.build(0, 0)

The first ``0`` prevent :func:`~Mole.build` dumping the input file. 
The second ``0`` prevent :func:`~Mole.build` parsing the command line arguments.


Access AO integrals
===================

PySCF uses `libcint <https://github.com/sunqm/libcint>`_ library as the AO
integral engine. A simple interface function :func:`Mole.intor` is provided 
to obtain the one- and two-electron AO integrals::

  kin = mol.intor('int1e_kin')
  vnuc = mol.intor('int1e_nuc')
  overlap = mol.intor('int1e_ovlp')
  eri = mol.intor('int2e')

For a full list of supported AO integrals, see :ref:`gto_moleintor`.
