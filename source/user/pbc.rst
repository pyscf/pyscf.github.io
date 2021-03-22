.. _user_pbc:

****************************
Periodic boundary conditions
****************************

PySCF supports electronic structure calculations of
extended systems with periodic boundary conditions (PBCs).  
PBC-specific functionality must be imported from :mod:`pyscf.pbc`,
which has a directory structure that mirrors that of the molecular :mod:`pyscf`
module, e.g. ``from pyscf.pbc import gto, scf``.
Details of how to specify the system details, including the unit cell,
basis sets, and pseudopotentials, are provided in
:ref:`user_pbc_gto`.

All electronic structure methods can be applied with periodic
boundary conditions at the Gamma point of the Brillouin zone, :math:`k=(0,0,0)`,
as described in :ref:`mix_mol`.  Converging to the thermodynamic limit
requires the use of larger and larger supercells, which might quickly become
prohibitive.  

More affordable convergence to the thermodynamic limit can be obtained with the
use of k-point sampling.  A calculation performed with a unit cell and multiple
uniformly distributed k-points is equivalent to one performed with a supercell
and a single k-point, but is significantly more affordable due to the explicitly
enforced translational symmetry (leading to crystal momentum conservation).
Only a subset of all methods in PySCF support k-point sampling and they are 
prefixed with "K", e.g., KHF, KCCSD, etc.

.. toctree::
   :maxdepth: 1

   pbc/gto.rst
   pbc/scf.rst
   pbc/df.rst
   pbc/mix_mol.rst
