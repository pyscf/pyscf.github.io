.. _overview:

An overview of PySCF
********************

Python-based simulations of chemistry framework (PySCF) is a general-purpose
electronic structure platform designed from the ground up to emphasize code
simplicity, so as to facilitate new method development and enable flexible
computational workflows. The package provides a wide range of tools to support
simulations of finite-size systems, extended systems with periodic boundary
conditions, low-dimensional periodic systems, and custom Hamiltonians, using
mean-field and post-mean-field methods with standard Gaussian basis functions.
To ensure ease of extensibility, PySCF uses the Python language to implement
almost all of its features, while computationally critical paths are
implemented with heavily optimized C routines. Using this combined Python/C
implementation, the package is as efficient as the best existing C or Fortran-
based quantum chemistry programs.


How to cite
===========
Bibtex entries::

  @article{PYSCF,
    title = {PySCF: the Python‐based simulations of chemistry framework},
    author = {Qiming Sun and Timothy C. Berkelbach and Nick S. Blunt and George H. Booth and Sheng Guo and Zhendong Li and Junzi Liu and James D. McClain and Elvira R. Sayfutyarova and Sandeep Sharma and Sebastian Wouters and Garnet Kin‐Lic Chan},
    year = {2017},
    journal = {Wiley Interdisciplinary Reviews: Computational Molecular Science},
    volume = {8},
    number = {1},
    pages = {e1340},
    doi = {10.1002/wcms.1340},
    url = {https://onlinelibrary.wiley.com/doi/abs/10.1002/wcms.1340},
    eprint = {https://onlinelibrary.wiley.com/doi/pdf/10.1002/wcms.1340},
  }

  @article{PYSCF,
    title = {Recent developments in the PySCF program package},
    author = {Sun,Qiming  and Zhang,Xing  and Banerjee,Samragni  and Bao,Peng  and Barbry,Marc  and Blunt,Nick S.  and Bogdanov,Nikolay A.  and Booth,George H.  and Chen,Jia  and Cui,Zhi-Hao  and Eriksen,Janus J.  and Gao,Yang  and Guo,Sheng  and Hermann,Jan  and Hermes,Matthew R.  and Koh,Kevin  and Koval,Peter  and Lehtola,Susi  and Li,Zhendong  and Liu,Junzi  and Mardirossian,Narbe  and McClain,James D.  and Motta,Mario  and Mussard,Bastien  and Pham,Hung Q.  and Pulkin,Artem  and Purwanto,Wirawan  and Robinson,Paul J.  and Ronca,Enrico  and Sayfutyarova,Elvira R.  and Scheurer,Maximilian  and Schurkus,Henry F.  and Smith,James E. T.  and Sun,Chong  and Sun,Shi-Ning  and Upadhyay,Shiv  and Wagner,Lucas K.  and Wang,Xiao  and White,Alec  and Whitfield,James Daniel  and Williamson,Mark J.  and Wouters,Sebastian  and Yang,Jun  and Yu,Jason M.  and Zhu,Tianyu  and Berkelbach,Timothy C.  and Sharma,Sandeep  and Sokolov,Alexander Yu.  and Chan,Garnet Kin-Lic},
    journal = {The Journal of Chemical Physics},
    volume = {153},
    number = {2},
    pages = {024109},
    year = {2020},
    doi = {10.1063/5.0006074},
    URL = {https://doi.org/10.1063/5.0006074},
    eprint = {https://doi.org/10.1063/5.0006074},
  }

In addition, if you use Libcint to compute integrals, please cite the following paper:

"Libcint: An efficient general integral library for Gaussian basis functions", Q. Sun, J. Comp. Chem. 36, 1664 (2015).

Features
========

.. include:: FEATURES

* Interface to integral package `Libcint <https://github.com/sunqm/libcint>`_

* Interface to DMRG `CheMPS2 <https://github.com/SebWouters/CheMPS2>`_

* Interface to DMRG `Block <https://github.com/sanshar/Block>`_

* Interface to FCIQMC `NECI <https://github.com/ghb24/NECI_STABLE>`_

* Interface to XC functional library `XCFun <https://github.com/dftlibs/xcfun>`_

* Interface to XC functional library `Libxc <http://www.tddft.org/programs/octopus/wiki/index.php/Libxc>`_

* Interface to tensor contraction library `TBLIS <https://github.com/devinamatthews/tblis>`_

* Interface to Heat-bath Selected CI program `Dice <https://sanshar.github.io/Dice/>`_

* Interface to geometry optimizer `Pyberny <https://github.com/jhrmnn/pyberny>`_

.. * Interface to `pyWannier90 <https://github.com/hungpham2017/pyWannier90>`_


.. include:: version.rst

.. include:: benchmark.rst
    
