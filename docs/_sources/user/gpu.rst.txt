.. _user_gpu:

GPU Acceleration (GPU4PySCF)
****************************

*Modules*: :py:mod:`gpu4pyscf`

.. module:: GPU4PySCF
   :synopsis: GPU4PySCF
.. sectionauthor:: Xiaojie Wu <wxj6000@gmail.com>.

Introduction
============

Modern GPUs accelerate quantum chemistry calculation significantly, but also have an advantage in cost saving `[1]`_.
Some of basic PySCF modules, such as SCF and DFT, are accelerated with GPU via a plugin package
GPU4PySCF (See the end of this page for the supported functionalities). For the density fitting scheme,
GPU4PySCF on A100-80G can be 1000x faster than PySCF on single-core CPU. The speedup of direct SCF scheme is relatively low.

.. _[1]: https://arxiv.org/abs/2404.09452

Installation
============
The binary package of GPU4PySCF is released based on the CUDA version.

.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

  * - CUDA version
    - GPU4PySCF
    - cuTensor
  * - CUDA 11.x
    - ``pip3 install gpu4pyscf-cuda11x``
    - ``pip3 install cutensor-cu11``
  * - CUDA 12.x
    - ``pip3 install gpu4pyscf-cuda12x``
    - ``pip3 install cutensor-cu12``

Usage of GPU4PySCF
==================
The GPU4PySCF APIs are designed to maintain compatibility with PySCF. The
classes and methods in GPU4PySCF are named identically to those in PySCF,
ensuring a familiar interface for users. However, GPU4PySCF classes do not
directly inherit from PySCF classes.

PySCF objects and GPU4PySCF objects can be converted to each other using the `to_gpu` and `to_cpu` methods.
The conversion process can automatically, recursively translate all attributes between GPU and CPU instances.
For example, numpy arrays on the CPU are converted into CuPy arrays on the GPU, and vice versa.
If certain attributes are exclusive to either the CPU or GPU, these attributes will be appropriately handled.
They are omitted or specifically converted, depending on the target platform.

There are two approaches to execute the computation on GPU.

1. Directly import GPU4PySCF classes and methods::

    import pyscf
    from gpu4pyscf.dft import rks

    atom ='''
    O       0.0000000000    -0.0000000000     0.1174000000
    H      -0.7570000000    -0.0000000000    -0.4696000000
    H       0.7570000000     0.0000000000    -0.4696000000
    '''

    mol = pyscf.M(atom=atom, basis='def2-tzvpp')
    mf = rks.RKS(mol, xc='LDA').density_fit()

    e_dft = mf.kernel()  # compute total energy
    print(f"total energy = {e_dft}")

    g = mf.nuc_grad_method()
    g_dft = g.kernel()   # compute analytical gradient

    h = mf.Hessian()
    h_dft = h.kernel()   # compute analytical Hessian

2. Convert PySCF object to the corresponding GPU4PySCF object with :func:`to_gpu`::

    import pyscf
    from pyscf.dft import rks

    atom ='''
    O       0.0000000000    -0.0000000000     0.1174000000
    H      -0.7570000000    -0.0000000000    -0.4696000000
    H       0.7570000000     0.0000000000    -0.4696000000
    '''

    mol = pyscf.M(atom=atom, basis='def2-tzvpp')
    mf = rks.RKS(mol, xc='LDA').density_fit().to_gpu()  # move PySCF object to GPU4PySCF object
    e_dft = mf.kernel()  # compute total energy


When the GPU task is done, the GPU4PySCF object can be converted into the corresponding PySCF object via :func:`mf.to_cpu()`.

In GPU4PySCF, wavefunctions, density matrices, and other array data are stored in CuPy arrays.
To transfer these data to NumPy arrays on the CPU, the `.get()` method of the CuPy array can be invoked.
For more detailed information on handling CuPy array conversions, please refer to the `CuPy APIs` documentation.

.. Cupy APIs: https://docs.cupy.dev/en/stable/user_guide/index.html

GPU4PySCF and PySCF Hybrid Programming
======================================
GPU4PySCF allows for seamless integration with existing PySCF programs, enabling
a hybrid approach that leverages both CPU and GPU resources in the program. This
integration is facilitated through the use of `to_gpu()` and `to_cpu()`
functions, which convert PySCF instances between CPU and GPU.

For instance, we can perform DFT calculations on GPU to obtain a set of DFT
orbitals followed by orbital localization using the Boys method on the CPU::

    import pyscf
    from pyscf import lo
    mol = pyscf.M(atom = '''
    O       0.0000000000    -0.0000000000     0.1174000000
    H      -0.7570000000    -0.0000000000    -0.4696000000
    H       0.7570000000     0.0000000000    -0.4696000000
    ''', basis='def2-tzvpp')

    # Perform DFT computation on GPU
    mf = mol.RKS(xc='b3lyp').to_gpu().run()

    # Transfer the computation back to CPU and continue the tasks on the CPU
    mf = mf.to_cpu()
    loc_orb = lo.Boys(mol, mf.mo_coeff[:,[2,3,4]]).kernel()

**GPU Implementation Availability**: The `to_gpu` method is implemented for
almost all methods in PySCF. However, the actual availability of GPU4PySCF
implementations for specific modules may vary. If a GPU4PySCF module is
available, `to_gpu` will return a GPU4PySCF instance. Otherwise, it will raise a
`NotImplementedError`.

Functionalities supported by GPU4PySCF
======================================
.. list-table::
   :widths: 25 25 25 25
   :header-rows: 1

  * - Method
    - SCF
    - Gradient
    - Hessian
  * - direct SCF
    - O
    - GPU
    - CPU
  * - density fitting
    - O
    - O
    - O
  * - LDA
    - O
    - O
    - O
  * - GGA
    - O
    - O
    - O
  * - mGGA
    - O
    - O
    - O
  * - hybrid
    - O
    - O
    - O
  * - unrestricted
    - O
    - O
    - O
  * - PCM solvent
    - GPU
    - GPU
    - FD
  * - SMD solvent
    - GPU
    - GPU
    - FD
  * - dispersion correction
    - CPU*
    - CPU*
    - FD
  * - nonlocal correlation
    - O
    - O
    - NA
  * - ECP
    - CPU
    - CPU
    - CPU
  * - MP2
    - GPU
    - CPU
    - CPU
  * - CCSD
    - GPU
    - CPU
    - NA
