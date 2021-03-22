
:mod:`gto` --- Molecular structure and GTO basis
************************************************
.. module:: gto
   :synopsis: Molecular structure and GTO basis

The :mod:`gto` module provides the functions to parse the command line options,
the molecular geometry and format the basic functions for ``libcint``
integral library.  In :mod:`gto.mole`, a basic class :class:`Mole` is
defined to hold the global parameters, which will be used throughout the
package.

Examples
========

* :source:`examples/gto/00-input_mole.py`
* :source:`examples/gto/01-input_geometry.py`
* :source:`examples/gto/02-dump_input.py`
* :source:`examples/gto/03-ghost_atom.py`
* :source:`examples/gto/04-input_basis.py`
* :source:`examples/gto/05-input_ecp.py`
* :source:`examples/gto/06-load_mol_from_chkfile.py`
* :source:`examples/gto/07-nucmod.py`
* :source:`examples/gto/09-apply_scf.py`
* :source:`examples/gto/10-atom_info.py`
* :source:`examples/gto/11-basis_info.py`
* :source:`examples/gto/12-serialization.py`
* :source:`examples/gto/13-symmetry.py`
* :source:`examples/gto/20-ao_integrals.py`
* :source:`examples/gto/20-ao_integrals_sph_to_spinor.py`
* :source:`examples/gto/20-soc_ao_integrals.py`
* :source:`examples/gto/20-soc_ecp.py`
* :source:`examples/gto/21-concatenate_molecules.py`
* :source:`examples/gto/22-range_separated_coulomb.py`
* :source:`examples/gto/23-orbitals_cart2sph.py`
* :source:`examples/gto/24-ao_value_on_grid.py`
* :source:`examples/gto/25-multipole-integrals.py`
* :source:`examples/gto/30-read_molpro_orb.py`

Program reference
=================

Mole class
----------

The :class:`Mole` class handles three layers: input, internal format, libcint arguments.
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

AO integration
--------------

.. automodule:: pyscf.gto.moleintor
   :members:


.. _gto_basis:

Basis set
---------

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

Some basis sets, *e.g.*, :source:`pyscf/gto/basis/dzp_dunning.py`, are saved in the
internal format.

.. automodule:: pyscf.gto.basis
   :members:

