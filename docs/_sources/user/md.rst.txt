.. _user_md:

Molecular Dynamics
******************

*Modules*: :py:mod:`pyscf.md`

Introduction
============

Using any available method for computing gradients, molecular dynamics can be run to evolve the nuclei along their potential energy surfaces.
The MD module in PySCF supports both NVE and NVT ensemble integrations. A simple example of an NVE Born-Oppenheimer MD simulation with CASSCF is given in :source:`examples/md/00-simple_nve.py`

.. literalinclude:: ../../examples/md/00-simple_nve.py

This will create 2 files, BOMD.md.data and BOMD.md.xyz that will contain the energy and structures for each frame along the trajectory. 
The BOMD.md.data will look like

.. code:: 

  time          Epot                 Ekin                 Etot               T
      0.00  -1.497083671804E+02  0.000000000000E+00  -1.497083671804E+02   0.0000
      5.00  -1.497083678218E+02  6.402095343905E-07  -1.497083671816E+02   0.0674
     10.00  -1.497083697399E+02  2.557532001000E-06  -1.497083671823E+02   0.2692
     15.00  -1.497083729245E+02  5.741058994089E-06  -1.497083671835E+02   0.6043
     20.00  -1.497083773578E+02  1.017280141026E-05  -1.497083671850E+02   1.0708
     25.00  -1.497083830147E+02  1.582769986352E-05  -1.497083671870E+02   1.6660
     30.00  -1.497083898631E+02  2.267377586371E-05  -1.497083671893E+02   2.3866
     35.00  -1.497083978644E+02  3.067236427680E-05  -1.497083671921E+02   3.2285
     40.00  -1.497084069741E+02  3.977851890862E-05  -1.497083671956E+02   4.1870
     45.00  -1.497084171401E+02  4.994098786544E-05  -1.497083671991E+02   5.2567

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


Integrators
===========

The NVE ensemble uses the Velocity-Verlet integration scheme. The positions are evolved according to:

.. math:: 
    \mathbf{r}(t_{i+1}) = \mathbf{r}(t_i) + (\delta t) \mathbf{v}(t_i) + \frac{1}{2}(\delta t)^2 \mathbf{a}(t_i)

and the velocities are evolved as 

.. math:: 
    \mathbf{v}(t_{i+1}) = \mathbf{v}(t_i) + (\delta t) \frac{\mathbf{a}(t_{i+1}) + \mathbf{a}(t_i)}{2}

Note that to compute the acceleration, the most common isotope masses are used (rather than the average isotope mass).

The NVT ensemble only has the Berendsen thermostat implemented. The Berendsen thermostat updates the velocities by weakly coupling the system to a heat bath at some temperature (:math:`T_0`). The temperature is therefore corrected such that the deviation exponentially decays with some time constant :math:`\tau`. :math:`\tau` can also be understood as the strength of the coupling between the system and bath. 

.. math::
   \frac{dT}{dt} = \frac{T_0 - T}{\tau}

At each time step, the velocities are rescaled according to 

.. math::
   \mathbf{v}_\mathrm{scaled}(t) = \mathbf{v}(t)\lambda

where :math:`\lambda` is given as

.. math::
   \lambda = \sqrt{1 + \frac{\delta t}{\tau}\left(\frac{T_0}{T} - 1\right)}

Note that we limit :math:`\lambda` to reasonable values. That is, if :math:`\lambda` is less 0.9, we set it to 0.9; and if :math:`\lambda` is greater than 1.1, we set it to 1.1. See :source:`examples/md/04-mb_nvt_berendson.py` for an example of how to use the Berendson theromstat.

Initial Velocity
================

The initial velocity for the system can be set through the `veloc` kwarg in the construction of the integrator such as

.. code::

  integrator = pyscf.md.NVE(hf, dt=5, steps=10, veloc=init_veloc)


This can be useful for restarting a dynamics simulation with a prior velocity, or starting a simulation with a specific temperature. Additionally, the MD module also contains a function to sample the velocities from a Maxwell-Boltzmann distribution at a given temperature (in Kelvin). 

.. code::

  init_veloc = md.distributions.MaxwellBotlzmannVelocity(mol, T=300)


The Maxwell-Botlzmann distribution at a specific temperature :math:`T` for a particle of mass :math:`m` is given as

.. math::
   f(v)dv = \sqrt{\frac{m}{2\pi k T}}\exp\left(-\frac{mv^2}{2kT}\right)

Since all of the distributions use the same pseudo-random number generator, the specific seed used can be set with the `md.set_seed()` function. The seed can also be set within the `PYSCF_CONFIG_FILE` as the `SEED` variable.
