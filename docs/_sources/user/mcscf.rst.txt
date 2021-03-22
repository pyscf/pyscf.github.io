.. _theory_mcscf:

Multi-configuration self-consistent field (MCSCF)
*************************************************

*Modules*: :mod:`mcscf`



Symmetry in CASSCF
------------------

.. literalinclude:: /../examples/mcscf/21-nosymhf_then_symcasscf.py

To restart a CASSCF calculation, you need prepare either CASSCF orbitals
or CI coefficients (not that useful unless doing a DMRG-CASSCF calculation) or
both.  For example:

.. literalinclude:: /../examples/mcscf/13-restart.py



 The third argument for CASCI/CASSCF is the size of CAS space; the
fourth argument is the number of electrons.  By default, the CAS
solver determines the alpha-electron number and beta-electron number
based on the attribute :attr:`Mole.spin`.  In the above example, the
number of alpha electrons is equal to the number of beta electrons,
since the ``mol`` object is initialized with ``spin=0``.  The spin
multiplicity of the CASSCF/CASCI solver can be changed by the fourth
argument::

  >>> mc = mcscf.CASSCF(m, 4, (4,2))
  >>> print('E(CASSCF) = %.9g' % mc.kernel()[0])
  E(CASSCF) = -149.609461
  >>> print('S^2 = %.7f, 2S+1 = %.7f' % mcscf.spin_square(mc))
  S^2 = 2.0000000, 2S+1 = 3.0000000

The two integers in the tuple represent the number of alpha and beta electrons.
Although it is a triplet state, the solution might not be correct since the
CASSCF is based on the incorrect singlet HF ground state.  Starting from the
ROHF ground state, we have::

  >>> mc = mcscf.CASSCF(rhf3, 4, 6)
  >>> print('E(CASSCF) = %.9g' % mc.kernel()[0])
  E(CASSCF) = -149.646746

The energy is lower than the RHF initial guess.
.. We can also use the UHF ground
.. state to start a CASSCF calculation::
.. 
..   >>> mc = mcscf.CASSCF(uhf3, 4, 6)
..   >>> print('E(CASSCF) = %.9g' % mc.kernel()[0])
..   E(CASSCF) = -149.661324
..   >>> print('S^2 = %.7f, 2S+1 = %.7f' % mcscf.spin_square(mc))
..   S^2 = 3.9713105, 2S+1 = 4.1091656
.. 
.. Woo, the total energy is even lower.  But the spin is contaminated.
