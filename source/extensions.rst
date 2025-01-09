.. _extensions_list:

Extensions and optional dependencies
************************************

PySCF is structured into three hierarchical levels of features: core
modules, forge modules, and extensions. In addition, some features are provided
as optional dependencies, enhancing PySCF functionality.

* **Core Modules**: These modules provide fundamental functionalities such as
  basic DFT computations and integral evaluations. They serve as the
  foundational building blocks for other packages and features. The APIs in
  these modules are relatively stable, and performance optimization is not the
  primary focus.

* **PySCF Forge**: This repository includes newly developed functionalities that
  are in the staging phase. Features in the forge may eventually be integrated
  into the core modules. For further details, refer to :ref:`pyscf_forge`.

* **Extensions**: Extensions are not intended to become part of the core
  modules. They are designed to accommodate experimental features, specific
  optimizations, or libraries that require complex compilation environments or
  dependencies. These extensions rely on the core modules for functionality. The
  PySCF organization commits to maintaining these extensions to ensure
  compatibility with other components of PySCF.

* **Optional Dependencies**: Many features and modules are developed by the
  community using the PySCF APIs. The key distinction between optional
  dependencies and extensions is that the PySCF organization does not maintain
  these optional dependencies.

List of Extensions
==================

(To be listed here:
 Names, brief descriptions, and installation instructions)

Summary
-------

=================== =========================================================
Project             URL
=================== =========================================================
cornell-shci        https://github.com/pyscf/cornell-shci
dispersion          https://github.com/pyscf/dispersion
dmrgscf             https://github.com/pyscf/dmrgscf
doci                https://github.com/pyscf/doci
fciqmc              https://github.com/pyscf/fciqmc
forge               https://github.com/pyscf/pyscf-forge
icmpspt             https://github.com/pyscf/icmpspt
mbd                 https://github.com/pyscf/mbd
naive-hci           https://github.com/pyscf/naive-hci
nao                 https://github.com/pyscf/nao
properties          https://github.com/pyscf/properties
qsdopt              https://github.com/pyscf/qsdopt
rt                  https://github.com/pyscf/rt
semiempirical       https://github.com/pyscf/semiempirical
shciscf             https://github.com/pyscf/shciscf
tblis               https://github.com/pyscf/pyscf-tblis
gpu4pyscf           https://github.com/pyscf/gpu4pyscf
=================== =========================================================

List of Optional Dependencies
=============================
(To be listed here:
 Names, brief descriptions, and links to the dependencies)

Summary
-------

=================== =========================================================
Project             URL
=================== =========================================================
cas_ac0             https://github.com/CQCL/pyscf-ac0
ccpy                https://github.com/piecuch-group/ccpy
cppe                https://github.com/maxscheurer/cppe
ebcc                https://github.com/BoothGroup/ebcc
dftd3               https://github.com/dftd3/simple-dftd3
pyqmc               https://github.com/WagnerGroup/pyqmc
zquatev             https://github.com/sunqm/zquatev
geometric           https://github.com/leeping/geomeTRIC
pyberny             https://github.com/jhrmnn/pyberny
=================== =========================================================
