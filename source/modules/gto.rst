.. _gto:

:mod:`gto` --- Molecular structure and GTO basis
************************************************

The :mod:`gto` module provides the functions to parse the command line options,
the molecular geometry and format the basic functions for ``libcint``
integral library.  In :mod:`gto.mole`, a basic class :class:`Mole` is
defined to hold the global parameters, which will be used throughout the
package.

Examples
========

:download:`examples/gto/00-input_mole.py </../examples/gto/00-input_mole.py>`
:download:`examples/gto/01-input_geometry.py </../examples/gto/01-input_geometry.py>`
:download:`examples/gto/02-dump_input.py </../examples/gto/02-dump_input.py>`
:download:`examples/gto/03-ghost_atom.py </../examples/gto/03-ghost_atom.py>`
:download:`examples/gto/04-input_basis.py </../examples/gto/04-input_basis.py>`
:download:`examples/gto/05-input_ecp.py </../examples/gto/05-input_ecp.py>`
:download:`examples/gto/06-load_mol_from_chkfile.py </../examples/gto/06-load_mol_from_chkfile.py>`
:download:`examples/gto/07-nucmod.py </../examples/gto/07-nucmod.py>`
:download:`examples/gto/09-apply_scf.py </../examples/gto/09-apply_scf.py>`
:download:`examples/gto/10-atom_info.py </../examples/gto/10-atom_info.py>`
:download:`examples/gto/11-basis_info.py </../examples/gto/11-basis_info.py>`
:download:`examples/gto/12-serialization.py </../examples/gto/12-serialization.py>`
:download:`examples/gto/13-symmetry.py </../examples/gto/13-symmetry.py>`
:download:`examples/gto/20-ao_integrals.py </../examples/gto/20-ao_integrals.py>`
:download:`examples/gto/20-ao_integrals_sph_to_spinor.py </../examples/gto/20-ao_integrals_sph_to_spinor.py>`
:download:`examples/gto/20-soc_ao_integrals.py </../examples/gto/20-soc_ao_integrals.py>`
:download:`examples/gto/20-soc_ecp.py </../examples/gto/20-soc_ecp.py>`
:download:`examples/gto/21-concatenate_molecules.py </../examples/gto/21-concatenate_molecules.py>`
:download:`examples/gto/22-range_separated_coulomb.py </../examples/gto/22-range_separated_coulomb.py>`
:download:`examples/gto/23-orbitals_cart2sph.py </../examples/gto/23-orbitals_cart2sph.py>`
:download:`examples/gto/24-ao_value_on_grid.py </../examples/gto/24-ao_value_on_grid.py>`
:download:`examples/gto/25-multipole-integrals.py </../examples/gto/25-multipole-integrals.py>`
:download:`examples/gto/30-read_molpro_orb.py </../examples/gto/30-read_molpro_orb.py>`

Input
=====

Geometry
--------

There are multiple ways to input molecular geometry. The internal format of
:attr:`Mole.atom` is a python list::

  atom = [[atom1, (x, y, z)],
          [atom2, (x, y, z)],
          ...
          [atomN, (x, y, z)]]

You can input the geometry in this format.  You can use Python script to
construct the geometry::

  >>> mol = gto.Mole()
  >>> mol.atom = [['O',(0, 0, 0)], ['H',(0, 1, 0)], ['H',(0, 0, 1)]]
  >>> mol.atom.extend([['H', (i, i, i)] for i in range(1,5)])

Besides Python list, tuple and numpy.ndarray are all supported by the internal
format::

  >>> mol.atom = (('O',numpy.zeros(3)), ['H', 0, 1, 0], ['H',[0, 0, 1]])

Also, :attr:`~Mole.atom` can be a string of Cartesian format or Z-matrix format::

  >>> mol = gto.Mole()
  >>> mol.atom = '''
  ... O 0 0 0
  ... H 0 1 0
  ... H 0 0 1;
  ... '''

There are a few requirements for the string format.  The string input
takes ``;`` or ``\n`` to partition atoms. White space and ``,`` are used to
split items for each atom.  Blank lines or lines started with ``#`` will be
ignored::

  >>> mol = gto.M(
  ... mol.atom = '''
  ... #O 0 0 0
  ... H 0 1 0
  ...
  ... H 0 0 1;
  ... ''')
  >>> mol.natm
  2

The geometry string is case-insensitive.  It also supports to input the nuclear
charges of elements::

  >>> mol = gto.Mole()
  >>> mol.atom = [[8,(0, 0, 0)], ['h',(0, 1, 0)], ['H',(0, 0, 1)]]

If you need to label an atom to distinguish it from the rest, you can prefix
or suffix number or special characters ``1234567890~!@#$%^&*()_+.?:<>[]{}|``
(except ``,`` and ``;``) to an atomic symbol.  With this decoration, you can
specify different basis sets, or masses, or nuclear models for different atoms::

  >>> mol = gto.Mole()
  >>> mol.atom = '''8 0 0 0; h:1 0 1 0; H@2 0 0'''
  >>> mol.basis = {'O': 'sto-3g', 'H': 'cc-pvdz', 'H@2': '6-31G'}
  >>> mol.build()
  >>> print(mol.atom)
  [['O', [0.0, 0.0, 0.0]], ['H:1', [0.0, 1.0, 0.0]], ['H@2', [0.0, 0.0]]]

No matter which format or symbols were used in the input, :func:`Mole.build`
will convert :attr:`Mole.atom` to the internal format::

  >>> mol.atom = '''
      O        0,   0, 0             ; 1 0.0 1 0

          H@2,0 0 1
      '''
  >>> mol.build()
  >>> print(mol.atom)
  [['O', [0.0, 0.0, 0.0]], ['H', [0.0, 1.0, 0.0]], ['H@2', [0.0, 0.0, 1.0]]]

In the program, the molecular geometry is accessed with :meth:`Mole.atom_coords`
function.  This function returns a (N,3) array for the coordinates (in Bohr) of
each atom::

  >>> print(mol.atom_coords())
  [[ 0.          0.          0.        ]
   [ 0.          1.88972612  0.        ]
   [ 0.          0.          1.88972612]]


.. _input_basis:

Input Basis
-----------
There are various ways to input basis sets.  Besides the input of universal
basis string and basis ``dict``::

  mol.basis = 'sto3g'
  mol.basis = {'O': 'sto3g', 'H': '6-31g'}

basis can be input with helper functions.
Function :func:`basis.parse` can parse a basis string of NWChem format
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

Functions :func:`basis.load` can be load arbitrary basis from the database, even
the basis which does not match the element.

  mol.basis = {'H': gto.basis.load('sto3g', 'C')}

Both :func:`basis.parse` and :func:`basis.load` return the basis set in the
internal format (See the :ref:`gto_basis`).

Basis parser supports "Ghost" atom::

  mol.basis = {'GHOST': gto.basis.load('cc-pvdz', 'O'), 'H': 'sto3g'}

More examples of ghost atoms in input can be found in
:download:`examples/gto/03-ghost_atom.py </../examples/gto/03-ghost_atom.py>`

Like the requirements of geometry input, you can use atomic symbol
(case-insensitive) or the atomic nuclear charge, as the keyword of the
:attr:`~Mole.basis` dict.  Prefix and suffix of numbers and special
characters are allowed.  If the decorated atomic symbol is appeared in
:attr:`~Mole.atom` but not :attr:`~Mole.basis`, the basis parser will
remove all decorations then seek the pure atomic symbol in
:attr:`~Mole.basis` dict.  In the following example, ``6-31G`` basis
will be assigned to the second H atom, but ``STO-3G`` will be used for
the third atom::

  mol.atom = '8 0 0 0; h1 0 1 0; H2 0 0 1'
  mol.basis = {'O': 'sto-3g', 'H': 'sto3g', 'H1': '6-31G'}


ECP
---
Effective core potential (ECP) can be specified with the attribute :attr:`ecp`
of :attr:`Mole`. Scalar type ECP is available for all molecular and crystal
methods. The built-in scalar ECP datasets include

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
Examples of ECP input can be found in :download:`examples/gto/05-input_ecp.py`.

Spin-orbit ECP integrals can be evaluated using PySCF integral driver. However,
SO-ECP is not automatically applied to any methods in current version. It needs
to be added to core Hamiltonian as shown in example
:download:`examples/gto/20-soc_ecp.py` and
:download:`examples/scf/44-soc_ecp.py`.
PySCF provides the following SOC-ECPs

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



Command line
------------

Some of the input variables can be passed from command line::

  $ python example.py -o /path/to/my_log.txt -m 1000

This command line specifies the output file and the maximum of memory for the
calculation.  By default, command line has the highest priority, which means our
settings in the script will be overwritten by the command line arguments.  To
make the input parser ignore the command line arguments, you can call the
:func:`Mole.build` with::

  mol.build(0, 0)

The first 0 prevent :func:`~Mole.build` dumping the input file.  The
second 0 prevent :func:`~Mole.build` parsing command line.


Spin and charge
---------------

Charge and spin multiplicity can be assigned to :class:`Mole` object::

  mol.charge = 1
  mol.spin = 1

Note :attr:`Mole.spin` is the number of unpaired electrons (difference between
the numbers of alpha electrons and beta electrons).  These two attributes do not
affect any other parameters in the :attr:`Mole.build` initialization function.
They can be set or modified (although not recommended to do so) after the
molecular object is initialized::

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

:attr:`Mole.charge` is the parameter to define the total number electrons of the
system.  In a custom system such as the Hubbard lattice model, the total number
of electrons needs to be defined directly::

  >>> mol = gto.Mole()
  >>> mol.nelectron = 10

Symmetry
--------

Point group symmetry information is held in :class:`Mole` object.  The symmetry
module (:ref:`symm`) of PySCF program can detect arbitrary point groups. In the
SCF calculations, PySCF program supports linear molecular symmetry
:math:`D_{\infty h}` (labelled as Dooh in the program), :math:`C_{\infty v}`
(labelled as Coov) plus :math:`D_{2h}` and its subgroups.

If the attribute :attr:`Mole.symmetry` is set, :meth:`Mole.build` function will
detect the top point group and the supported subgroups.  The detected point
groups are saved in :attr:`Mole.topgroup` and :attr:`Mole.groupname`::

  >>> mol = gto.Mole()
  >>> mol.atom = 'B 0 0 0; H 0 1 1; H 1 0 1; H 1 1 0'
  >>> mol.symmetry = True
  >>> mol.build()
  >>> print(mol.topgroup)
  C3v
  >>> print(mol.groupname)
  Cs

Sometimes it is necessary to use a lower symmetry instead of the detected
symmetry group.  The subgroup symmetry can be specified in
:attr:`Mole.symmetry_subgroup` and the program will first detect the highest
possible symmetry group and lower the point group symmetry to the specified
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

When a particular symmetry was assigned to :attr:`Mole.symmetry`,
the initialization function :meth:`Mole.build` will test whether the molecule
geometry is subject to the required symmetry.  If not, initialization will stop
and an error message will be issued::

  >>> mol = gto.Mole()
  >>> mol.atom = 'N 0 0 0; N 0 0 1'
  >>> mol.symmetry = 'Dooh'
  >>> mol.build()
  RuntimeWarning: Unable to identify input symmetry Dooh.
  Try symmetry="Dooh" with geometry (unit="Bohr")
  ('N', [0.0, 0.0, -0.9448630622825309])
  ('N', [0.0, 0.0, 0.9448630622825309])

.. note::
  :attr:`Mole.symmetry_subgroup` does not have effects
  when specific symmetry group is assigned to :attr:`Mole.symmetry`.

When symmetry is enabled in the :class:`Mole` object, the point group symmetry
information will be used to construct the symmetry adapted orbital basis (see
also :ref:`symm`).  The symmetry adapted orbitals are held in
:attr:`Mole.symm_orb` as a list of 2D arrays.  Each element of the list
is an AO (atomic orbital) to SO (symmetry-adapted orbital) transformation matrix
of an irreducible representation.  The name of the irreducible representations
are stored in :attr:`Mole.irrep_name` and their internal IDs (see more details
in :ref:`symm`) are stored in :attr:`Mole.irrep_id`::

  >>> mol = gto.Mole()
  >>> mol.atom = 'N 0 0 0; N 0 0 1'
  >>> mol.symmetry = True
  >>> mol.build()
  >>> for s,i,c in zip(mol.irrep_name, mol.irrep_id, mol.symm_orb):
  ...     print(s, i, c.shape)
  A1g 0 (10, 3)
  E1gx 2 (10, 1)
  E1gy 3 (10, 1)
  A1u 5 (10, 3)
  E1uy 6 (10, 1)
  E1ux 7 (10, 1)


Program reference
=================

mole
----

:class:`Mole` class handles three layers: input, internal format, libcint arguments.
The relationship of the three layers are::

  .atom (input)  <=>  ._atom (for python) <=> ._atm (for libcint)
  .basis (input) <=> ._basis (for python) <=> ._bas (for libcint)

input layer does not talk to libcint directly.  Data are held in python
internal fomrat layer.  Most of methods defined in this class only operates
on the internal format.  Exceptions are make_env, make_atm_env, make_bas_env,
:func:`set_common_orig_`, :func:`set_rinv_orig_` which are used to
manipulate the libcint arguments.


.. automodule:: pyscf.gto.mole
   :members:

.. autoclass:: Mole
   :members:


.. _gto_moleintor:

moleintor
---------

.. automodule:: pyscf.gto.moleintor
   :members:


.. _gto_basis:

basis
-----

Internal format
^^^^^^^^^^^^^^^

This module loads basis set and ECP data from basis database and parse the basis
(mostly in NWChem format) and finally convert to internal format.  The internal
format of basis set is::

  basis = {atom_type1:[[angular_momentum
                        (GTO-exp1, contract-coeff11, contract-coeff12),
                        (GTO-exp2, contract-coeff21, contract-coeff22),
                        (GTO-exp3, contract-coeff31, contract-coeff32),
                        ...],
                       [angular_momentum
                        (GTO-exp1, contract-coeff11, contract-coeff12),
                        ...],
                       ...],
           atom_type2:[[angular_momentum, (...),],
                       ...],

For example::

  mol.basis = {'H': [[0,
                      (19.2406000, 0.0328280),
                      (2.8992000, 0.2312080),
                      (0.6534000, 0.8172380),],
                     [0,
                      (0.1776000, 1.0000000),],
                     [1,
                      (1.0000000, 1.0000000),]],
              }

Some basis sets, *e.g.*, :file:`pyscf/gto/basis/dzp_dunning.py`, are saved in the
internal format.

.. automodule:: pyscf.gto.basis
   :members:


