# PySCF API Reference

This manual documents the exposed API functions, modules in PySCF.
For learning how to use PySCF, see the User Guide 



```{toctree}
:caption: Python APIs
:maxdepth: 1

general

naming conventions

global_configurations
  - PYSCF_DISABLE_EVAL environment variable

logger

gto

chkfile

SCF
   how to custom SCF class, which methods need to be overridden

XC and NumInt overriding
  - which methods required by fock build in SCF
  - which methods required by gradients computation
  - higher order derivatives:
    - energy_nuclear_gradients
```

```{toctree}
:caption: C API
:maxdepth: 1

C extensions are not supposed to be accessed directly. There is no guarantee of
stable C APIs
```
