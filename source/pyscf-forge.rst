.. _pyscf_forge:

PySCF-Forge
***********

PySCF-Forge is a platform designed to host innovative methodologies that may
eventually be integrated into the main PySCF core branch.

Installation
============
PySCF-forge can be installed via PyPI::

    $ pip install pyscf-forge

Alternatively, you can build it from the source code::

    $ pip install git+https://github.com/pyscf/pyscf-forge

Using modules in PySCF-Forge
============================
Once you have installed both PySCF and PySCF-Forge, you can access the features
of PySCF-forge just like the regular features in PySCF, using `from pyscf import
...` statement. Although these source code of the modules are hosted in a
separate repository, they function as part of the PySCF package.

Please note that the APIs and usage of modules within PySCF-Forge may be subject to
change. Particularly, when a module is transferred from PySCF-Forge to the core
branch, modifications to module names or function APIs may occur.

Contributing to PySCF-Forge
===========================
Contributing to the PySCF-Forge repository provides an opportunity for your
methodologies to be reviewed and potentially integrated into the PySCF core
branch. By placing your code into PySCF-forge,
PySCF maintainers will provide feedback and help evaluate the integration
potential of your work.

It is important to note that not all code submitted to PySCF-Forge will be moved
to the core branch. Such decisions require approval from the board. For a
detailed guide on contributing, please refer to the
[Guidelines](https://github.com/pyscf/pyscf-forge/blob/master/CONTRIBUTING.md).

How to transfer a feature from forge to core?
=============================================
To ensure a smooth transition of features from PySCF-Forge to the PySCF core
branch, the following steps should be followed:

* Duration Requirement: A feature must have been part of PySCF-Forge for at
  least six months before it can be considered for transfer to the core branch.
  This period allows for adequate testing and refinement.

* Initiating Transfer: Developers who wish to transfer a feature should
  open an issue on the PySCF-Forge GitHub repository to request the migration.
  This proposal will initiate a discussion in the PySCF board to consider the transfer.

* Maintenance Commitment: Identify and propose a volunteer who commits to
  being responsible for the long-term maintenance of the feature.

Co-authorship
=============
When a feature is added to the PySCF-forge, it is crucial to recognize the
contributions of individuals. Here are the guidelines for co-authorship in the
context of the PySCF project:

* Initial Contribution: The original contributor of a feature to the
  PySCF-forge is eligible for co-authorship.

* Long-term Maintenance: Individuals who commit to the long-term
  maintenance of the feature also qualify for co-authorship. Their ongoing
  efforts to update, optimize, and support the feature are critical.

* Documentation and Testing: Documenting the feature and creating comprehensive
  tests are important to the stability and usability of the feature.
  Contributors who play significant roles in these works are qualified for
  co-authorship.
