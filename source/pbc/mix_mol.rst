.. _mix_mol:

Mixing PBC and molecular modules
********************************
 
Post-HF methods, as standalone numerical solvers, do not require
knowledge of the boundary conditions.  Calculations on finite-sized systems
and extended systems are distinguished by the boundary conditions of the integrals (and
basis).  The same post-HF solver can thus be used for both finite-size problems
and the periodic boundary problems if they have a compatible Hamiltonian
structure.

In PySCF, many molecular post-HF solvers have two implementations: an incore and
outcore version.  These differ by the treatment of the 2-electron
integrals.  The incore solver takes the :attr:`_eri` (or :attr:`with_df`, see
:ref:`mol_df`) from the underlying mean-field object as the two-electron
interaction part of the Hamiltonian while the outcore solver generates the
2-electron integrals (with free boundary conditions) on the fly.
To use the molecular post-HF solvers in PBC code, we need to ensure that the incore
version solver is called.

Generating :attr:`_eri` in a mean-field object is the straightforward way to
trigger the incore post-HF solver.  If the allowed memory is large enough to
hold the entire 2-electron integral array, the Gamma point HF solver always
generates and holds this array.  A second choice is to set :attr:`incore_anyway`
in ``cell`` which forces the program to generate and hold :attr:`_eri` in
the mean-field object.

.. note::

  If the problem is big, :attr:`incore_anyway` may overflow the available
  physical memory.

Holding the full integral array :attr:`_eri` in memory limits the problem size
one can treat.  Using the density fitting object :attr:`with_df` to hold the
integrals can overcome this problem.  This architecture has been bound to PBC
and molecular mean-field modules. Not all post-HF methods are available with density fitting.

Aside from the 2-electron integrals, there are some attributes and methods
required by the post-HF solver.  They are :meth:`get_hcore`, and
:meth:`get_ovlp` for 1-electron integrals, :attr:`_numint`, :attr:`grids` for
the numerical integration of DFT exchange-correlation functionals.  These are all
overloaded in the PBC mean-field object to produce the PBC integrals. 

Examples
--------

.. literalinclude:: ../../../examples/pbc/12-gamma_point_post_hf.py
