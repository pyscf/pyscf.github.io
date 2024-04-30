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
GPU4PySCF

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
GPU4PySCF APIs are designed to be compatible with PySCF. When supported, high-level functions and objects are named the same as PySCF. But, GPU4PySCF classes do not directly inherit from PySCF class.
PySCF objects and GPU4PySCF objects can be converted into each other by :func:`to_gpu` and :func:`to_cpu`. In the conversion, the numpy arrays will be converted into cupy array. And the functions will be omitted if they are not supported with GPU acceleration.

One can use the two modes to accelerate the calculations: directly use GPU4PySCF object::

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

Alternatively, one can convert PySCF object to the corresponding GPU4PySCF object with :func:`to_gpu` since PySCF 2.5.0 ::

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
Then, more sophisticated methods in PySCF can apply. One can also convert the individual CuPy array to numpy array with `Cupy APIs`_.

.. Cupy APIs: https://docs.cupy.dev/en/stable/user_guide/index.html

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
