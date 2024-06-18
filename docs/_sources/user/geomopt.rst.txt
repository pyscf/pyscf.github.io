.. _user_geomopt:

Geometry optimization
*********************

*Modules*: :py:mod:`pyscf.geomopt`

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
    print(mol_eq.tostring())

    # pyberny
    from pyscf.geomopt.berny_solver import optimize
    mol_eq = optimize(mf, maxsteps=100)
    print(mol_eq.tostring())


The second way is to create an :func:`optimizer` from the :class:`Gradients` class::

    # geometric
    mol_eq = mf.Gradients().optimizer(solver='geomeTRIC').kernel()
    print(mol_eq.tostring())

    # pyberny
    mol_eq = mf.Gradients().optimizer(solver='berny').kernel()
    print(mol_eq.tostring())

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
be specified in a text file. Its format can be found
in the `online documentation <https://geometric.readthedocs.io/en/latest/constraints.html#input-format>`_
or in the `template file <https://github.com/leeping/geomeTRIC/blob/master/examples/constraints.txt>`_.
Then the name of the constraints file needs to be passed to PySCF::

    params = {"constraints": "constraints.txt",}
    mol_eq = optimize(mf, **params)
    mol_eq = mf.Gradients().optimizer(solver='geomeTRIC').kernel(params)

The geometry optimization can also be carried out based on
custom energy and gradient functions:

.. literalinclude:: /../examples/geomopt/02-as_pyscf_method.py

Transition state optimization
-----------------------------
Transition state optimization are available in ``geomeTRIC`` and ``qsdopt``.

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

When using `geomeTRIC` to search transition state, you can simply add the
keyword `transition` in geomeTRIC input configuration to trigger the TS search
module::

    params = {'transition': True}
    mol_ts = mf.Gradients().optimizer(solver='geomeTRIC').kernel(params)

    params = {'transition': True, 'trust': 0.02, 'tmax': 0.06}
    mol_ts = mf.Gradients().optimizer(solver='geomeTRIC').kernel(params)

Keywords `trust` and `tmax` control the trust region in TS search. They are
not mandatory input settings. The default settings in
geomeTRIC library are chosen to maximize job success rates. Tuning trust region
may bring small improvements to the TS optimization performance.
For more detailed discussions of available options,
we refer to the `online documentation of geomeTRIC transition states <https://geometric.readthedocs.io/en/latest/transition.html>`_.

Except the initial step, geomeTRIC library does not support the ability to read
analytical Hessian from quantum chemistry calculations at runtime.
The initial Hessian can be fed to the optimizer after enabling the keyword
`hessian`::

    params = {'transition': True, 'hessian': True}
    mol_ts = mf.Gradients().optimizer(solver='geomeTRIC').kernel(params)

The interface implemented in PySCF can generate a temporary file to save the
initial Hessian and pass to geomeTRIC. This keyword will be ignored if the
analytical Hessian of the underlying method was not implemented.

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
