# PySCF API Reference

This manual describes PySCF interfaces, their expected behavior, and extension
points. For guidance on running calculations, see the [User Guide](../user/index.md).

```{toctree}
:caption: Python APIs
:maxdepth: 1

general
global_configurations
logger
NumPy helpers <numpy_helper>
gto
Checkpoint files <chkfile>
scf
pbc
xc_numint
```

## C API

PySCF's C extensions are intended for internal use. Their interfaces are not
guaranteed to remain stable.
