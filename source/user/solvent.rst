
.. _user_solvent:

Solvation models
****************

*Modules*: :mod:`solvent`

.. module:: solvent
   :synopsis: Solvation models and solvent effects
.. sectionauthor:: Qiming Sun <osirpt.sun@gmail.com>.

Introduction
============

Solvation model allows the quantum chemistry calculations to include the
interactions between solvents and the quantum solute. Solvents can be treated
implicitly, known as continuum solvents, and explicitly. For continuum solvents,
we implemented the ddCOSMO (domain-decomposition COSMO solvation model). For
the explicit solvent environment, we provided the interface to call the
polarizable embedding library CPPE.


ddCOSMO
=======

Self-consistent solvents for ground state
-----------------------------------------

Solvent model can be applied on to an SCF object::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RKS(xc='b3lyp').DDCOSMO().run()

In regular MCSCF (CASCI or CASSCF), and post-SCF calculations, the setup for
self-consistent solvent is similar::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mc = mol.CASCI(4, 4).DDCOSMO().run()

  mp2_model = mol.MP2().DDCOSMO().run()


Solvent for excited states
--------------------------
When combining to TDDFT or other methods of excited states, solvent can be
modelled in the manner of self-consistency (fast solvent) or single shot (slow
solvent). Below we use TDDFT to demonstrate the treatments of fast solvent
and slow solvent.

In vertical excitations, the solvent almost does not respond to the change of
electronic structure. It should be viewed as slow solvent. The calculation
can be started with an SCF with fully relaxed solvent and followed by a regular
TDDFT method::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RHF().ddCOSMO().run()
  td = mf.TDA().run()

In the diabatic excitations, we would like to let the solvent rapidly responds
to the electronic structure of excited states. The entire solvation system
should converge to an equilibrium state between solvent and the excited state of
the solute. In this scenario, solvent model should be applied on to the excited
state methods::

  mf = mol.RHF().ddCOSMO().run()
  td = mf.TDA().ddCOSMO()
  td.with_solvent.equilibrium_solvation = True
  td.kernel()

Please note that the flag ``equilibrium_solvation`` needs to be set to ``True``
in this case. PySCF by default assumes the slow solvent model for TDDFT.

In the complicated procedure which involves for example electronic states from
different states (typically in the MCSCF calculations with state-average or
state-specific approximations), PySCF ddCOSMO implementation allows to input a
density matrix and freeze the solvent equilibrated against the input density
matrix::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mc = mol.CASCI(4, 4).DDCOSMO()
  mc.fcisolver.nstates = 5
  mc.with_solvent.state_id = 1  # Slow solvent wrt the first excited state
  mc.run()

The slow solvent does not have to be corresponding to a particular state. It can
be even the solvent from a different geometry or an artificial quantum state of
solute::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  scf_dm = mol.RHF().run().make_rdm1()

  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    1.035   -1.082
       H  0.   -1.035   -1.082''',
                basis='6-31g*', verbose=4)
  mc = mol.CASCI(4, 4).DDCOSMO(dm=scf_dm).run()


Solvent parameters
------------------
The default solvent in the ddCOSMO module is water. When studying other types of
solvents, you can consider to modify the dielectric parameter ``eps`` using the
constants listed below::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RHF().ddCOSMO()
  mf.with_solvent.eps = 32.613   # methanol                          
  mf.run()

These dielectric constants are obtained from https://gaussian.com/scrf/.
More dataset can be found in Minnesota Solvent Descriptor Database
(https://comp.chem.umn.edu/solvation)

================================== ====================
Solvent                            dielectric constant
================================== ====================
Water                              78.3553
Acetonitrile                       35.688
Methanol                           32.613
Ethanol                            24.852
IsoQuinoline                       11.00
Quinoline                          9.16
Chloroform                         4.7113
DiethylEther                       4.2400
Dichloromethane                    8.93
DiChloroEthane                     10.125
CarbonTetraChloride                2.2280
Benzene                            2.2706
Toluene                            2.3741
ChloroBenzene                      5.6968
NitroMethane                       36.562
Heptane                            1.9113
CycloHexane                        2.0165
Aniline                            6.8882
Acetone                            20.493
TetraHydroFuran                    7.4257
DiMethylSulfoxide                  46.826
Argon                              1.430
Krypton                            1.519
Xenon                              1.706
n-Octanol                          9.8629
1,1,1-TriChloroEthane              7.0826
1,1,2-TriChloroEthane              7.1937
1,2,4-TriMethylBenzene             2.3653
1,2-DiBromoEthane                  4.9313
1,2-EthaneDiol                     40.245
1,4-Dioxane                        2.2099
1-Bromo-2-MethylPropane            7.7792
1-BromoOctane                      5.0244
1-BromoPentane                     6.269
1-BromoPropane                     8.0496
1-Butanol                          17.332
1-ChloroHexane                     5.9491
1-ChloroPentane                    6.5022
1-ChloroPropane                    8.3548
1-Decanol                          7.5305
1-FluoroOctane                     3.89
1-Heptanol                         11.321
1-Hexanol                          12.51
1-Hexene                           2.0717
1-Hexyne                           2.615
1-IodoButane                       6.173
1-IodoHexaDecane                   3.5338
1-IodoPentane                      5.6973
1-IodoPropane                      6.9626
1-NitroPropane                     23.73
1-Nonanol                          8.5991
1-Pentanol                         15.13
1-Pentene                          1.9905
1-Propanol                         20.524
2,2,2-TriFluoroEthanol             26.726
2,2,4-TriMethylPentane             1.9358
2,4-DiMethylPentane                1.8939
2,4-DiMethylPyridine               9.4176
2,6-DiMethylPyridine               7.1735
2-BromoPropane                     9.3610
2-Butanol                          15.944
2-ChloroButane                     8.3930
2-Heptanone                        11.658
2-Hexanone                         14.136
2-MethoxyEthanol                   17.2
2-Methyl-1-Propanol                16.777
2-Methyl-2-Propanol                12.47
2-MethylPentane                    1.89
2-MethylPyridine                   9.9533
2-NitroPropane                     25.654
2-Octanone                         9.4678
2-Pentanone                        15.200
2-Propanol                         19.264
2-Propen-1-ol                      19.011
3-MethylPyridine                   11.645
3-Pentanone                        16.78
4-Heptanone                        12.257
4-Methyl-2-Pentanone               12.887
4-MethylPyridine                   11.957
5-Nonanone                         10.6
AceticAcid                         6.2528
AcetoPhenone                       17.44
a-ChloroToluene                    6.7175
Anisole                            4.2247
Benzaldehyde                       18.220
BenzoNitrile                       25.592
BenzylAlcohol                      12.457
BromoBenzene                       5.3954
BromoEthane                        9.01
Bromoform                          4.2488
Butanal                            13.45
ButanoicAcid                       2.9931
Butanone                           18.246
ButanoNitrile                      24.291
ButylAmine                         4.6178
ButylEthanoate                     4.9941
CarbonDiSulfide                    2.6105
Cis-1,2-DiMethylCycloHexane        2.06
Cis-Decalin                        2.2139
CycloHexanone                      15.619
CycloPentane                       1.9608
CycloPentanol                      16.989
CycloPentanone                     13.58
Decalin-mixture                    2.196
DiBromomEthane                     7.2273
DiButylEther                       3.0473
DiEthylAmine                       3.5766
DiEthylSulfide                     5.723
DiIodoMethane                      5.32
DiIsoPropylEther                   3.38
DiMethylDiSulfide                  9.6
DiPhenylEther                      3.73
DiPropylAmine                      2.9112
e-1,2-DiChloroEthene               2.14
e-2-Pentene                        2.051
EthaneThiol                        6.667
EthylBenzene                       2.4339
EthylEthanoate                     5.9867
EthylMethanoate                    8.3310
EthylPhenylEther                   4.1797
FluoroBenzene                      5.42
Formamide                          108.94
FormicAcid                         51.1
HexanoicAcid                       2.6
IodoBenzene                        4.5470
IodoEthane                         7.6177
IodoMethane                        6.8650
IsoPropylBenzene                   2.3712
m-Cresol                           12.44
Mesitylene                         2.2650
MethylBenzoate                     6.7367
MethylButanoate                    5.5607
MethylCycloHexane                  2.024
MethylEthanoate                    6.8615
MethylMethanoate                   8.8377
MethylPropanoate                   6.0777
m-Xylene                           2.3478
n-ButylBenzene                     2.36
n-Decane                           1.9846
n-Dodecane                         2.0060
n-Hexadecane                       2.0402
n-Hexane                           1.8819
NitroBenzene                       34.809
NitroEthane                        28.29
n-MethylAniline                    5.9600
n-MethylFormamide-mixture          181.56
n,n-DiMethylAcetamide              37.781
n,n-DiMethylFormamide              37.219
n-Nonane                           1.9605
n-Octane                           1.9406
n-Pentadecane                      2.0333
n-Pentane                          1.8371
n-Undecane                         1.9910
o-ChloroToluene                    4.6331
o-Cresol                           6.76
o-DiChloroBenzene                  9.9949
o-NitroToluene                     25.669
o-Xylene                           2.5454
Pentanal                           10.0
PentanoicAcid                      2.6924
PentylAmine                        4.2010
PentylEthanoate                    4.7297
PerFluoroBenzene                   2.029
p-IsoPropylToluene                 2.2322
Propanal                           18.5
PropanoicAcid                      3.44
PropanoNitrile                     29.324
PropylAmine                        4.9912
PropylEthanoate                    5.5205
p-Xylene                           2.2705
Pyridine                           12.978
sec-ButylBenzene                   2.3446
tert-ButylBenzene                  2.3447
TetraChloroEthene                  2.268
TetraHydroThiophene-s,s-dioxide    43.962
Tetralin                           2.771
Thiophene                          2.7270
Thiophenol                         4.2728
trans-Decalin                      2.1781
TriButylPhosphate                  8.1781
TriChloroEthene                    3.422
TriEthylAmine                      2.3832
Xylene-mixture                     2.3879
z-1,2-DiChloroEthene               9.2
================================== ====================


Polarizable embedding
=====================
To use polarizable embedding model for mean-field calculations, one would need
to first generate potential data for the input of CPPE library. The best way to
generate potential files is with `PyFraME <https://gitlab.com/FraME-projects/PyFraME>`_.
You can directly throw in a pdb file, select the QM region and how to parametrize
different parts of the environment (with either pre-defined potentials, or with LoProp).
Some guidance is also provided in the Tutorial Review paper about PE, section 4:
https://onlinelibrary.wiley.com/doi/full/10.1002/qua.25717
Therein, the format of the potential file is also explained (it’s the same
format as used in the original Dalton pelib implementation).

With the generated potential file, one can carry out the polarizable embedding
calculations::

  import pyscf
  mol = pyscf.M(atom='''
       C  0.    0.      -0.542
       O  0.    0.       0.677
       H  0.    0.935   -1.082
       H  0.   -0.935   -1.082''',
                basis='6-31g*', verbose=4)
  mf = mol.RKS(xc='b3lyp')
  mf = pyscf.solvent.PE(mf, 'potfile')
  mf.run()


References
==========

.. bibliography:: ref_solvent.bib
   :style: unsrt
