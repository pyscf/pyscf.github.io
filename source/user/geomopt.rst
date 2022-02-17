.. _user_geomopt:

Geometry optimization
*********************

*Modules*: :mod:`geomopt`

Basics
------

PySCF implements geometry optimization via
interfaces to `geomeTRIC <https://github.com/leeping/geomeTRIC>`_ and
`PyBerny <https://github.com/jhrmnn/pyberny>`_, and through PySCF extension
`qsdopt <https://github.com/pyscf/qsdopt>`_(see :ref:`installing`
for installation instructions).

There are two ways to invoke geometry optimization.
The first is to import the :func:`optimize` function 
from the respective modules, i.e., :mod:`pyscf.geomopt.geometric_solver` 
and :mod:`pyscf.geomopt.berny_solver`::

    from pyscf import gto, scf
    mol = gto.M(atom='N 0 0 0; N 0 0 1.2', basis='ccpvdz')
    mf = scf.RHF(mol)

    # geometric
    from pyscf.geomopt.geometric_solver import optimize
    mol_eq = optimize(mf, maxsteps=100)
    print(mol_eq.atom_coords())

    # pyberny
    from pyscf.geomopt.berny_solver import optimize
    mol_eq = optimize(mf, maxsteps=100)
    print(mol_eq.atom_coords())


The second way is to create an :func:`optimizer` from the :class:`Gradients` class::

    # geometric
    mol_eq = mf.Gradients().optimizer(solver='geomeTRIC').kernel()
    print(mol_eq.atom_coords())

    # pyberny
    mol_eq = mf.Gradients().optimizer(solver='berny').kernel()
    print(mol_eq.atom_coords())

For the ``geomeTRIC`` backend, the convergence criteria are 
controlled by the following parameters::

    conv_params = { # These are the default settings
        'convergence_energy': 1e-6,  # Eh
        'convergence_grms': 3e-4,    # Eh/Bohr
        'convergence_gmax': 4.5e-4,  # Eh/Bohr
        'convergence_drms': 1.2e-3,  # Angstrom
        'convergence_dmax': 1.8e-3,  # Angstrom
    }
    mol_eq = optimize(mf, **conv_params)
    mol_eq = mf.Gradients().optimizer(solver='geomeTRIC').kernel(conv_params)

For the ``PyBerny`` backend, the convergence criteria are
controlled by the following parameters::

    conv_params = {  # These are the default settings
        'gradientmax': 0.45e-3,  # Eh/[Bohr|rad]
        'gradientrms': 0.15e-3,  # Eh/[Bohr|rad]
        'stepmax': 1.8e-3,       # [Bohr|rad]
        'steprms': 1.2e-3,       # [Bohr|rad]
    }
    mol_eq = optimize(mf, **conv_params)
    mol_eq = mf.Gradients().optimizer(solver='berny').kernel(conv_params)


Constraints
-----------

``geomeTRIC`` supports user defined constraints. The constraints can 
be specified in a text file with the format described
`here <https://github.com/leeping/geomeTRIC/blob/master/examples/constraints.txt>`_.
One needs to pass the name of this file to PySCF::

    params = {"constraints": "constraints.txt",}
    mol_eq = optimize(mf, **params)
    mol_eq = mf.Gradients().optimizer(solver='geomeTRIC').kernel(params)

The geometry optimization can also be carried out based on
custom energy and gradient functions:

.. literalinclude:: /../examples/geomopt/02-as_pyscf_method.py

Transition state optimization
--------------
The PySCF extension `qsdopt <https://github.com/pyscf/qsdopt>`_ performs transition
state optimizations through
`quadratic steepest descent method <https://aip.scitation.org/doi/10.1063/1.467721>`_
This is a second order method that requires computation of the hessian at some steps during
the optimization process. The following is a minimal usage example::

    from pyscf imort gto, scf
    from pyscf.qsdopt.qsd_optimizer import QSD

    mol = gto.M(atom='''
    O 0 0 0
    H 0 0 1.2
    H 0, 0.5, -1.2''',
    basis='minao', verbose=0, unit="Bohr")
    mf = scf.RHF(mol)

    optimizer = QSD(mf, stationary_point="TS")
    optimizer.kernel()

Several keyword arguments can be passed to `kernel`:

- hess_update_freq: Frequency for numerical reevaluation of hessian. = 0 evaluates the numerical
  hessian in the first iteration and is updated with an BFGS rule, unless approaching a trap region,
  where it is reevaluated. (Default: 0)
- numhess_method: Method for evaluating numerical hessian. Forward and central differences are available
  with "forward" and "central", respectively. (Default: "forward")
- max_iter: Maximum number of optimization steps. (Default: 100)
- step: Maximum norm between two optimization steps. (Default: 0.1)
- hmin: Minimum distance between current geometry and stationary point of the quadratic form
  to consider convergence reached. (Default: 1e-6)
- gthres: Gradient norm to consider convergence. (Default: 1e-5)

.. literalinclude:: /../examples/geomopt/16-ethane_transition_state.py

Excited states
--------------
For excited-state geometry optimizations, the state to be optimized 
needs to be specified in the respective :class:`Gradients` objects:

.. literalinclude:: /../examples/geomopt/12-excited_states.py

For examples of state-specific and state-averaged CASSCF geometry optimizations,
see :source:`examples/geomopt/12-mcscf_excited_states.py`.

Callback
--------

Callback functions can be invoked at each optimization step.
The following example shows how to add charge analysis for each 
geometry during the optimization.

.. literalinclude:: /../examples/geomopt/20-callback.py
