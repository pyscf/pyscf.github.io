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


.. _version:

Version history
===============

===============  ==========
---------------  ----------
1.7.6a1_         2020-12-17
1.7.5_           2020-10-04
1.7.4_           2020-08-02
1.7.3_           2020-06-10
1.7.2_           2020-05-13
1.7.1_           2020-02-29
1.7.0_           2020-01-01
1.6.6_           2020-01-01
1.6.5            2019-11-17
1.6.4            2019-09-14
1.6.3            2019-07-28
1.6.2            2019-06-17
1.6.1            2019-03-15
1.6              2018-12-31
1.6 beta         2018-11-26 (feature frozen)
1.6 alpha        2018-08-15
1.5.5_           2018-12-31
1.5.4            2018-11-16
1.5.3            2018-09-06
1.5.2            2018-08-15
1.5.1            2018-07-01
1.5              2018-06-08
1.5 beta         2018-04-15 (feature frozen)
1.5 alpha        2018-03-21
1.4.7_           2018-04-15
1.4.4            2018-03-20
1.4              2017-10-05
1.4 beta         2017-08-22 (feature frozen)
1.4 alpha        2017-06-24
1.3.5_           2017-08-12
1.3.4            2017-08-09
1.3.3            2017-07-05
1.3.2            2017-06-05
1.3.1            2017-05-14
1.3              2017-04-25
1.3 beta         2017-02-15 (feature frozen)
1.3 alpha 2      2017-01-04
1.3 alpha 1      2016-12-04
1.2.3            2017-04-24
1.2.2            2017-02-15
1.2.1            2017-01-26
1.2              2016-11-07
1.2 beta         2016-09-13 (feature frozen)
1.2 alpha        2016-08-05
1.1              2016-06-04
1.1 beta         2016-04-11 (feature frozen)
1.1 alpha 2      2016-03-08
1.1 alpha 1      2016-02-08
1.0              2015-10-07
1.0 rc1          2015-09-07
1.0 beta 1       2015-08-02 (feature frozen)
1.0 alpha 2      2015-07-03
1.0 alpha 1      2015-04-07
===============  ==========

.. _1.7.6a1: https://github.com/pyscf/pyscf/releases/tag/v1.7.6a1
.. _1.7.5: https://github.com/pyscf/pyscf/releases/tag/v1.7.5
.. _1.7.4: https://github.com/pyscf/pyscf/releases/tag/v1.7.4
.. _1.7.3: https://github.com/pyscf/pyscf/releases/tag/v1.7.3
.. _1.7.2: https://github.com/pyscf/pyscf/releases/tag/v1.7.2
.. _1.7.1: https://github.com/pyscf/pyscf/releases/tag/v1.7.1
.. _1.7.0: https://github.com/pyscf/pyscf/releases/tag/v1.7.0
.. _1.6.6: https://github.com/pyscf/pyscf/releases/tag/v1.6.6
.. _1.5.5: https://github.com/pyscf/pyscf/releases/tag/v1.5.5
.. _1.4.7: https://github.com/pyscf/pyscf/releases/tag/v1.4.7
.. _1.3.5: https://github.com/pyscf/pyscf/releases/tag/v1.3.5


.. _benchmark:

Benchmark
=========

.. ========= =========================
.. ========= =========================
.. CPU         Intel i5 @ 3.1 GB
.. Memory      16 GB DDR3 @ 1333 MHz
.. OS          Debian 6.0
.. BLAS        MKL 10.3
.. Compiler    gcc-4.4
.. ========= =========================
.. 
.. N2, on 1 CPU core:
.. 
.. ================ ========= =============
..  Basis            cc-pVTZ   ANO-Roos-TZ
.. ================ ========= =============
..  HF               0.31 s    3.10 s
..  density fit HF   0.96 s    1.58 s
..  B3LYP            1.12 s    4.32 s
..  MP2              0.05 s    0.24 s
..  CASSCF(4,4)      0.50 s    2.24 s
..  CCSD             1.61 s    7.87 s
.. ================ ========= =============
.. 
.. Benzene, on 1 CPU core:

========= ==============================
Platform
========= ==============================
CPU       4 Intel E5-2670 @ 2.6 GB
Memory    64 GB DDR3
OS        Custom Redhat 6.6
BLAS      MKL 11.0
Compiler  Intel 13.0
========= ==============================

Benzene, on 16 CPU cores

================ ========= ========= =============
 Basis            6-31G**   cc-pVTZ   ANO-Roos-TZ
================ ========= ========= =============
 HF               0.55 s     5.76 s   389.1 s
 density fit HF   3.56 s     7.61 s    13.8 s
 B3LYP            3.84 s    11.44 s   360.2 s
 MP2              0.21 s     4.66 s   115.9 s
 CASSCF(6,6)      2.88 s    34.73 s   639.7 s
 CCSD             18.24 s   477.0 s   6721 s 
================ ========= ========= =============

C60, on 16 CPU cores

================ ========= =========
 Basis            6-31G**   cc-pVTZ 
================ ========= =========
 HF               1291 s    189 m
 SOSCF (newton)             77 m
 density fit HF   316.7 s   43.3 m
================ ========= =========

Fe(II)-porphyrin (FeC20H12N4), on 16 CPU cores

================ ========= ========= =========
 Basis            cc-pVDZ   cc-pVTZ   cc-pVQZ
================ ========= ========= =========
 SOSCF (newton)   193.4 s   20.1 m    127.1 m
 CASSCF(10,10)    1808 s    241 m
 CASSCF(11,8)     763.8 s   150.3 m   1280 m
================ ========= ========= =========

