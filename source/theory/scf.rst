.. _theory_scf:

***********************************
Self-consistent field (SCF) methods
***********************************

*Modules*: :ref:`scf <scf>`, :ref:`pbc.scf <pbc_scf>`, :ref:`soscf <soscf>` 

Introduction
============
In self-consistent field (SCF) methods, the electron interactions are treated in a mean-field way.
Here, the SCF methods include both Hartree-Fock (HF) theory and Kohn-Sham (KS) density functional theory (DFT).
This Chapter summarizes the general SCF capabilities of PySCF. 
For more details specific to DFT, see :numref:`theory_dft`.

A minimal example of using the SCF module is as follows::

    from pyscf import gto, scf
    mol = gto.M(
        atom = 'H 0 0 0; F 0 0 1.1',  # in Angstrom
        basis = 'ccpvdz',
        symmetry = True,
    )
    mf = scf.HF(mol)
    mf.kernel()

This will run a HF calculation for the hydrogen fluoride molecule using the default SCF settings.


Theory
======
In HF and KS-DFT methods, the ground-state wavefunction is approximated by a single Slater determinant,

.. math::

   |\Phi_0\rangle = \frac{1}{\sqrt{N!}}
   \begin{vmatrix} 
   \psi_1(\mathbf{r}_1) &\psi_2(\mathbf{r}_1) &\dots  &\psi_N(\mathbf{r}_1)\\
   \psi_1(\mathbf{r}_2) &\psi_2(\mathbf{r}_2) &\dots  &\psi_N(\mathbf{r}_2)\\
   \vdots               &\vdots               &\ddots &\vdots\\
   \psi_1(\mathbf{r}_N) &\psi_2(\mathbf{r}_N) &\dots  &\psi_N(\mathbf{r}_N)
   \end{vmatrix} \;,

where the molecular orbitals (MO) :math:`|\psi_i\rangle` are obtained by solving the Hartree-Fock equtions:

.. math::

   \hat{F}|\psi_i\rangle = \varepsilon_i |\psi_i\rangle \;.

The Fock operator :math:`\hat{F}`, when represented within the atomic orbital (AO) basis, has the following form 
(in a spin unrestricted formalism):

.. math::

   F_{\mu\nu}^{\alpha} = h_{\mu\nu} + J_{\mu\nu} - K_{\mu\nu}^{\alpha} \;,

where 

.. math::

   h_{\mu\nu} = \left( \mu | \hat{T} |\nu \right) + \left( \mu | \hat{V}_{nuc} |\nu \right) \;, 

usually called the core Hamiltonian, is the sum of kinetic energy operator and nuclear attraction operator, and

.. math::

   J_{\mu\nu} = \sum_{\lambda\sigma} \left(\mu\nu|\lambda\sigma\right) \left(P_{\lambda\sigma}^{\alpha}+P_{\lambda\sigma}^{\beta}\right)

and 

.. math::

   K_{\mu\nu}^{\alpha} = \sum_{\lambda\sigma} \left(\mu\lambda|\nu\sigma\right) P_{\lambda\sigma}^{\alpha}  

define the Coulomb and exchange operators, respectively.
In the equations above, :math:`\mathbf{P}` labels the density matrix:

.. math::

   P_{\mu\nu}^{\alpha} = \sum_{i\in occ} C_{\mu i}^{\alpha} C_{i\nu}^{\alpha\dagger} \;,

where :math:`\mathbf{C}` represents the MO coefficients:

.. math::

   |\psi^{\alpha}\rangle  = \sum_{\mu} |\mu\rangle C_{\mu i}^{\alpha} \;.

Within the AO representation, the Hartree-Fock equations reduce to the Roothaan-Hall or Pople-Nesbet :cite:`pople1954self` equations:

.. math::

   \mathbf{F}^{\alpha} \mathbf{C}^{\alpha} = \bf{\varepsilon}^{\alpha} \mathbf{S} \mathbf{C}^{\alpha} \;,

where 

.. math::

   S_{\mu\nu} = \left( \mu | \nu \right) 

is the AO overlap matrix. 
Note that as the Coulomb and exchange operators depend on the electron density, 
solving the Hartree-Fock equation is a nonlinear procedure,
which requires iteratively updating the MOs and the Fock operator until reaching self-consistency.
Finally, with the converged MOs or density matrix, one may compute various ground-state properties.
For example, the ground-state electronic energy is expressed as

.. math::

   E_{\rm HF} = \frac{1}{2} \left\{{\rm Tr}[\mathbf{h}(\mathbf{P}^{\alpha}+\mathbf{P}^{\beta})]
              + {\rm Tr}(\mathbf{F}^{\alpha}\mathbf{P}^{\alpha}) + {\rm Tr}(\mathbf{F}^{\beta}\mathbf{P}^{\beta}) \right\} \;. 


Periodic boundary conditions
----------------------------


Initial guess
=============
As the Roothaan-Hall and Pople-Nesbet equations are solved iteratively, 
an initial guess for the MOs or the density matrices must be supplied.
Poor initial guess may cause slow convergence or even divergence of the procedure. 
Furthermore, when treating magnetic or open-shell systems, 
the initial guess must be carefully chosen in order to get the correct state.

There are several options available in PySCF for selecting the initial guess to solve the 
SCF problem. One can set the attribute :attr:`mf.init_guess`
to the following values to generate the initial guess in different ways:

* ``'minao'`` (default)

    The initial guess density matrix is first generated based on the atomic natural orbital (ANO) basis 
    :cite:`widmark1990density,roos2004relativistic,roos2004main,roos2005new,roos2005new_a,roos2008new`,
    then projected onto the basis set used for the SCF calculation.

* ``'hcore'``

    The core Hamiltonian is diagonalized to get the initial MOs. 

* ``'atom'``

    The initial guess density matrix is from the superposition of atomic HF
    density matrix. Commonly know as the 'SAD' method.

* ``'chk'``

    Read the existing SCF results from the checkpoint file, then the density matrix is projected onto the
    basis set used for the new SCF calculation.

Alternatively, the user could manually set the initial guess density matrix for an SCF calculation 
by using the ``'dm0'`` argument. 
For example, the followings script first computes the HF density matrix for :math:`\rm Cr^{6+}` cation,  
which is then used as the initial guess for the HF calculation of :math:`\rm Cr` atom. ::

    #
    # use cation to produce initial guess
    #
    mol = gto.Mole()
    mol.build(
        symmetry = 'D2h',
        atom = [['Cr',(0, 0, 0)], ],
        basis = 'cc-pvdz',
        charge = 6,
        spin = 0,
    )

    mf = scf.RHF(mol)
    mf.kernel()
    dm1 = mf.make_rdm1()

    mol.charge = 0
    mol.spin = 6
    mol.build(False,False)

    mf = scf.RHF(mol)
    mf.kernel(dm0=dm1)

More examples can be found in 
:download:`examples/scf/15-initial_guess.py </../examples/scf/15-initial_guess.py>` 

Accelerating SCF convergence
============================

Direct Inversion in the Iterative Subspace (DIIS)
-------------------------------------------------
At convergence of an SCF calcuation, one should expect the density matrix commute with 
the Fock matrix:

.. math::

   \mathbf{SPF} - \mathbf{FPS} = \mathbf{0} \;.

Prior to convergence, it is possible to define an error vector as

.. math::

   \mathbf{e}_i \equiv \mathbf{S}\mathbf{P}_i\mathbf{F}_i - \mathbf{F}_i\mathbf{P}_i\mathbf{S} \;,

where :math:`\mathbf{F}_i` is a linear combination of the Fock matrices in the previous SCF cycles:

.. math::

   \mathbf{F}_i = \sum_{k=i-L}^{i-1} c_k \mathbf{F}_k \;,

:math:`\mathbf{P}_i` is obtained by diagonalizing :math:`\mathbf{F}_i`, and
:math:`L` is the size of the DIIS subspace, which can be modified by setting the :attr:`mf.diis_space` attribute 
(the default size is 8).
The DIIS method :cite:`pulay1980convergence,pulay1982improved` 
minimizes the square of the error vector 
with respect to the DIIS coefficients :math:`c_k`
under the constraint that :math:`\sum_k c_k = 1`.
The Euler–Lagrange equation of such a constrained minimization problem reads:

.. math::

   \left( 
   \begin{array}{cccc} 
   \mathbf{e}_1\cdot\mathbf{e}_1  &\dots  &\mathbf{e}_1\cdot\mathbf{e}_L  &1      \\ 
   \vdots                         &\ddots &\vdots                         &\vdots \\
   \mathbf{e}_L\cdot\mathbf{e}_1  &\dots  &\mathbf{e}_L\cdot\mathbf{e}_L  &1      \\
   1                              &\dots  &1                              &0
   \end{array}
   \right) \left( 
   \begin{array}{c}
   c_1    \\
   \vdots \\
   c_{L}  \\
   \lambda
   \end{array} 
   \right) = \left(
   \begin{array}{c}
   0       \\
   \vdots  \\
   0       \\
   1
   \end{array}
   \right) 

PySCF also implements two other similar DIIS algorithms, 
namely, EDIIS :cite:`kudin2002black` and ADIIS :cite:`hu2010accelerating`. 
Interested readers should refer to the reference.
An example of selecting different DIIS schemes can be found in 
:download:`examples/scf/24-tune_diis.py </../examples/scf/24-tune_diis.py>`

Co-iterative augmented hessian (CIAH) second order SCF solver :cite:`sun2016co`
-------------------------------------------------------------------------------

References
==========
.. bibliography:: ref_scf.bib
   :style: unsrt
