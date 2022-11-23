.. _user_md:

Molecular Dynamics
******************

*Modules*: :mod:`md`

Introduction
============

Using any available method for computing gradients, molecular dynamics can be run to evolve the nuclei along their potential energy surfaces.
The MD module in PySCF supports the Velocity-Verlet integration scheme in the NVE ensemble. The positions are evolved according to:

.. math:: 
    \mathbf{r}(t_{i+1}) = \mathbf{r}(t_i) + (\delta t) \mathbf{v}(t_i) + \frac{1}{2}(\delta t)^2 \mathbf{a}(t_i)

and the velocities are evolved as 

.. math:: 
    \mathbf{v}(t_{i+1}) = \mathbf{v}(t_i) + (\delta t) \frac{\mathbf{a}(t_{i+1}) + \mathbf{a}(t_i)}{2}

Note that to compute the acceleration, the most common isotope masses are used (rather than the average isotope mass).
A simple example of a Born-Oppenheimer MD simulation with CASSCF is given in :source:`examples/md/00-simple_nve.py`

.. literalinclude:: ../../examples/md/00-simple_nve.py

This will create 2 files, BOMD.md.energies and BOMD.md.xyz that will contain the energy and structures for each frame along the trajectory. 
The BOMD.md.energies will look like

.. code:: 

    time          Epot                 Ekin                 Etot
    0.00  -1.497083671814E+02  0.000000000000E+00  -1.497083671814E+02
    5.00  -1.497083678220E+02  6.403560409065E-07  -1.497083671817E+02
   10.00  -1.497083697402E+02  2.557810719094E-06  -1.497083671824E+02
   15.00  -1.497083729250E+02  5.741509150521E-06  -1.497083671835E+02
   20.00  -1.497083773584E+02  1.017340316838E-05  -1.497083671850E+02
   25.00  -1.497083830154E+02  1.582841563588E-05  -1.497083671870E+02
   30.00  -1.497083898640E+02  2.267458853778E-05  -1.497083671894E+02
   35.00  -1.497083978654E+02  3.067323581032E-05  -1.497083671922E+02
   40.00  -1.497084069752E+02  3.977943396938E-05  -1.497083671958E+02
   45.00  -1.497084171413E+02  4.994201038141E-05  -1.497083671993E+02

and the BOMD.md.xyz will look like

.. code:: 

    2
    MD Time 0
    O           0.00000        0.00000        0.00000
    O           0.00000        0.00000        1.20000
    2
    MD Time 5
    O          -0.00000       -0.00000       -0.00001
    O           0.00000        0.00000        1.20001
    2
    MD Time 10
    O          -0.00000        0.00000       -0.00002
    O           0.00000       -0.00000        1.20002
    2
    MD Time 15
    O           0.00000        0.00000       -0.00006
    O          -0.00000       -0.00000        1.20006
    2
    MD Time 20
    O           0.00000        0.00000       -0.00010
    O          -0.00000       -0.00000        1.20010
    2
    MD Time 25
    O           0.00000        0.00000       -0.00015
    O          -0.00000       -0.00000        1.20015
    2
    MD Time 30
    O           0.00000        0.00000       -0.00022
    O          -0.00000       -0.00000        1.20022
    2
    MD Time 35
    O           0.00000        0.00000       -0.00030
    O          -0.00000       -0.00000        1.20030
    2
    MD Time 40
    O           0.00000        0.00000       -0.00039
    O          -0.00000       -0.00000        1.20039
    2
    MD Time 45
    O           0.00000        0.00000       -0.00050
    O          -0.00000       -0.00000        1.20050

