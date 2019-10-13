.. PySCF documentation master file, created by
   sphinx-quickstart on Thu Jan 15 01:55:04 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PySCF documentation!
===============================
 
PySCF is a collection of electronic structure programs powered by Python.
The package aims to provide a simple, light-weight, and efficient platform for
quantum chemistry calculations and code development.  The program is developed
with the following principles:

* Easy to install, to use, to extend and to be embedded;

* Minimal requirements on libraries (no Boost or MPI) and computing
  resources (perhaps sacrificing efficiency to reduce I/O);

* 90/10 Python/C (only computational hot spots are written in C);

* 90/10 functional/OOP (unless performance critical, functions are pure).

In addition to the core libraries, PySCF supports a rich
ecosystem of plugins and external modules that, for example, provide
MPI versions of some routines, additional quantum chemistry methods and analysis, interface with quantum computing toolkits etc. See :ref:'.. _installing_plugin:'.

  
Contents
--------

.. toctree::
   :maxdepth: 2
   :numbered:

   overview.rst
   version.rst
   install.rst
   code-rule.rst
   benchmark.rst
   tutorial.rst
   gto.rst
   ao2mo.rst
   df.rst
   symm.rst
   scf.rst
   dft.rst
   dftd3.rst
   grad.rst
   hessian.rst
   soscf.rst
   sgx.rst
   tddft.rst
   tdscf.rst
   mp.rst
   ci.rst
   doci.rst
   cc.rst
   gw.rst
   fci.rst
   hci.rst
   cornell_shci.rst
   mcscf.rst
   dmrgscf.rst
   fciqmcscf.rst
   mrpt.rst
   icmpspt.rst
   x2c.rst
   pbc.rst
   qmmm.rst
   solvent.rst
   semiempirical.rst
   geomopt.rst
   lo.rst
   tools.rst
   advanced.rst
   ccn.rst
   nao.rst
   lib.rst
   prop.rst
   data.rst
   shciscf.rst

You can also download the `PDF version
<http://www.sunqm.net/pyscf/files/pdf/PySCF-1.7.pdf>`_ of this manual.


.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


