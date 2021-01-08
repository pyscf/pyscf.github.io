.. _theory_pbc_gto:

Crystalline Gaussian-type atomic orbitals
*****************************************

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
