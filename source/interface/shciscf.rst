.. _shciscf:

:mod:`shciscf` --- Semistochastic heat bath configuration interaction (SHCI) 
****************************************************************************

The :mod:`shciscf` module interfaces the `Dice <https://sanshar.github.io/Dice/index.html>`_ program to PySCF.
Dice is an efficient implementation for the SHCI algorithm.
It can be used with the CASCI and CASSCF module to solve large active space problems.  
The path of Dice program and other configurations should be initialized in the configuration file
``pyscf/shciscf/settings.py`` before using the SHCI method.

Examples
========
Examples can be found in the `Dice documentation <https://sanshar.github.io/Dice/usingincasscf.html?highlight=shciscf>`_.

Program reference
=================

.. automodule:: pyscf.shciscf
 
.. automodule:: pyscf.shciscf.shci
   :members:

