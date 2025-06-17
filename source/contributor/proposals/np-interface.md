# Introducing a universal NumPy interface for various NumPy-compatible libraries

A general NumPy wrapper can be developed in PySCF to provide a universal interface
to access NumPy functions. This interface would enable automatic switching to
`jax.numpy` or `torch` functions depending on the runtime context.

This interface can be imported as

```python
import pyscf.numpy as fnp
```

