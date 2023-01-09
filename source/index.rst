.. PySCF

Home
****

.. image:: ../logo/pyscf-logo.png
   :width: 400

The Python-based Simulations of Chemistry Framework (PySCF) is an open-source
collection of electronic structure modules powered by Python. The package
provides a simple, lightweight, and efficient platform for quantum chemistry
calculations and methodology development.
PySCF can be used to simulate the properties of molecules, crystals, and
custom Hamiltonians using mean-field and post-mean-field methods.
To ensure ease of extensibility, almost all of the features in PySCF are
implemented in Python, while computationally critical parts are
implemented and optimized in C. Using this combined Python/C
implementation, the package is as efficient as the best existing C or Fortran
based quantum chemistry programs.
In addition to its core libraries, PySCF supports a rich
ecosystem of :ref:`installing_extproj`.

.. MPI versions of some routines, additional quantum chemistry methods and analysis, interface with quantum computing toolkits *etc*. See :ref:`installing_extproj`.

.. For installation instructions, see the :ref:`install` page.

.. For a guide to performing a variety of calculations with PySCF, see the :ref:`quickstart` guide.


.. COMMENTED OUT FOR NOW
.. Install <install.rst>
.. User Guide <user.rst>
.. Developer Guide <develop.rst>
.. API Docs <api_docs/pyscf.rst>
.. Blog <blog_wrapper.rst>
.. Modules <modules.rst>
.. about.md

.. toctree::
   :hidden:
   :maxdepth: 1

   Install <install.rst>
   Quickstart <quickstart.rst>
   User Guide <user.rst>
   Developer Guide <develop.rst>
   API <pyscf_api_docs/pyscf.rst>
   About <about.rst>

.. You can also download the `PDF version <http://www.sunqm.net/pyscf/files/pdf/PySCF-1.7.pdf>`_ of this manual.


.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`
