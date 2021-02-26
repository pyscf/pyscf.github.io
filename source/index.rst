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
MPI versions of some routines, additional quantum chemistry methods and analysis, interface with quantum computing toolkits *etc*.
See :ref:`installing_plugin`.

.. toctree::
   :maxdepth: 1
   :numbered:

   overview.rst
   install.rst
   tutorial.rst
   user.rst
   develop.rst
   modules.rst
   interface.rst

You can also download the `PDF version
<http://www.sunqm.net/pyscf/files/pdf/PySCF-1.7.pdf>`_ of this manual.


.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
