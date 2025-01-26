.. _api_v3:

# PySCF v3 API proposals

.. This page lists proposed changes to the API in PySCF-v3

## Removing TaggedNPArray

### Passing additional keyword arguments

### Returning multiple values instead of the TaggedNPArray

### Global configuration to control the return style

---

## Introducing a universal NumPy interface for various NumPy-compatible libraries

A general NumPy wrapper can be developed in PySCF to provide a universal interface
to access NumPy functions. This interface would enable automatic switching to
`jax.numpy` or `torch` functions depending on the runtime context.

This interface can be imported as::

  import pyscf.numpy as fnp


---

## Allowing extensions written in C++
