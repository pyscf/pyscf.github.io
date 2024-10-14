# How to use PySCF 

This page provides an introduction to the generic organization of PySCF
and typical workflows.

## Modules, classes, and the kernel method

Similar to [NumPy](https://numpy.org/) or [SciPy](https://scipy.org/), PySCF is
a collection of modules, such as `gto` (for defining molecules with Gaussian type orbitals), 
`scf` (for self-consistent field calculations), or `cc` (for coupled-cluster calculations). 
Modules must be imported to be used,
```python
from pyscf import gto, scf, cc
```

Modules provide access to both functions and classes, where the latter is more commonly
used to define a calculation. For example, the `gto` module provides the `gto.Mole` class,
the `scf` module provides the `scf.RHF` class (and others, such as `scf.UHF`, etc.), and
the `cc` module provides the `cc.CCSD` class.

Performing a calculation in PySCF typically involves
importing a module, instantiating a class provided by that module with some arguments, 
and executing the functions of that class. For example,
```python
from pyscf import scf  # import module
myhf = scf.RHF(...)  # instantiate class
e_hf = myhf.kernel()  # execute kernel() method to do the calculation
```

Every class has the `kernel()` method, which serves as the driver of
the calculation, although many classes provide an alias to the `kernel()` method,
such as the `build()` method of the `gto.Mole` class. 

Once an object is created, you can always call `kernel()` to start or restart
a calculation.  The return value of the kernel method depends on the class.

The instance of one class is commonly passed as an argument to instantiate the
next class in a workflow. For example, the instance of the molecular structure
class is passed to instantiate a Hartree-Fock class, whose instance is passed
to instantiate a coupled-cluster class,
```python
from pyscf import gto, scf, cc
mymol = gto.Mole(...)
mymol.build()  # returns mymol
myscf = scf.RHF(mymol)
e_hf = mymol.kernel()
mycc = cc.CCSD(myscf)
e_corr, t1, t2 = mycc.kernel()
```

Chained calculations, like the one above, can also be performed more concisely
using [Stream methods](#stream-methods) or [Scanners](#scanners), as described
in the following sections.

## Stream methods

To unify the return value of different methods and thus allow chaining calculations together, 
PySCF includes three "stream methods". A stream method of an object only returns the object
itself. The three stream methods are described below.

1.  The `set` method updates object attributes. For example,
    ```python
    mf = scf.RHF(mol).set(conv_tol=1e-5)
    ```
    is identical to two lines of statements,
    ```python
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-5
    ```

2.  The `run` method calls the `kernel` method.  Arguments passed to the
    `run` method will be passed to the kernel method.  If keyword arguments are 
    given, `run` will first call the `set` method to update the attributes and 
    then execute the `kernel` method.  For example,
    ```python
    mf = scf.RHF(mol).run(dm_init, conv_tol=1e-5)
    ```
    is identical to three lines of statements,
    ```python
    mf = scf.RHF(mol)
    mf.conv_tol = 1e-5
    mf.kernel(dm_init)
    ```

3.  The `apply` method passes the current object (as the first argument) to the
    given function/class and returns a new object.  If arguments and keyword
    arguments are given, they will all be passed to the function/class. For
    example,
    ```python
    mc = mol.apply(scf.RHF).run().apply(mcscf.CASSCF, 6, 4, frozen=4)
    ```  
    is identical to,
    ```python
    mf = scf.RHF(mol)
    mf.kernel()
    mc = mcscf.CASSCF(mf, 6, 4, frozen=4)
    ```
    Note that the `apply()` method does not call the `kernel()` method.

In addition to these three stream methods, many regular class methods also
return the object (especially those that do not have any particular values to
return). Such methods can therefore be used in streams. For
example,
```python  
dm = gto.M(atom='H 0 0 0; H 0 0 1') \
  .apply(scf.RHF) \
  .dump_flags() \
  .run() \
  .make_rdm1()
```
This code works because the `dump_flags()` method simply prints information and then
returns the object.

## Scanners

A scanner is a function that takes a `Mole` (or `Cell`) object as input and
returns the energy or nuclear gradients at a chosen level of theory. A scanner
can be considered as a shortcut function for a sequence of statements, which
includes the initialization of a required calculation model with possible
precomputing, updating the attributes based on the settings of the referred
object, calling the kernel function, and finally returning results. 

For example, consider the following conventional script to perform a potential
energy surface scan of the dissociation of the hydrogen molecule using CCSD,
```python
for r in (1.0, 1.1, 1.2):
    mol = gto.M(atom=f"H 0 0 0; H 0 0 {r}")
    mf = scf.RHF(mol).run()
    mycc = cc.CCSD(mf).run()
    print(mycc.e_tot)
```

This can be simplified using the `as_scanner()` method,
```python
cc_scanner = gto.M().apply(scf.RHF).apply(cc.CCSD).as_scanner()
for r in (1.0, 1.1, 1.2):
    print(cc_scanner(gto.M(atom=f"H 0 0 0; H 0 0 {r}")))
```

There are two types of scanners available in the package: energy scanners and
nuclear gradients scanners. An energy scanner, like the example above, only
returns the energy of the given molecular structure while the nuclear gradients
scanner returns the nuclear gradients.

A scanner is a special derived object of the calling class. Most methods that are
defined in the calling are also accessible through the scanner object. For example,
```python
mf_scanner = gto.M().apply(scf.RHF).as_scanner()
mf_scanner(gto.M(atom='H 0 0 0; H 0 0 1.2'))
mf_scanner.analyze()
dm1 = mf_scanner.make_rdm1()

mf_grad_scanner = mf_scanner.nuc_grad_method().as_scanner()
mf_grad_scanner(gto.M(atom='H 0 0 0; H 0 0 1.2'))
```
As shown in this example, the scanner behaves very similarly to an RHF 
class object, except that the scanner does not need the `kernel` or `run`
methods to run a calculation.  Given a molecule structure, the scanner
automatically checks and updates the necessary object dependencies and passes the
work flow to the `kernel` method.  The computational results are held in the
scanner object the same way as in the regular class object.

To make the behavior of scanner objects uniform for all levels of theory, 
two attributes (`e_tot` and `converged`) are defined for all energy scanners,
and three attributes (`e_tot`, `de`, and `converged`) are defined for
all nuclear gradients scanners.



## Class and function behaviors

Classes are designed to hold only the final results (such as energies and wavefunction parameters) 
and the control parameters (such as the convergence threshold and the maximum number of iterations).
Intermediate quantities are **not** saved in the class. 

After calling the `kernel()` or `run()` method, results will be generated and
saved as attributes of the object.  For example,
```python
from pyscf import gto, scf, cc
mol = gto.M(atom='H 0 0 0; H 0 0 1.1', basis='ccpvtz')
mf = scf.RHF(mol).run()
mycc = ccsd.CCSD(mf).run()
print(mycc.e_tot)
print(mycc.e_corr)
print(mycc.t1.shape)
print(mycc.t2.shape)
```

Many useful functions are defined at both the class level (as methods) and
the module level. For example, 
```python
myhf = scf.RHF(mol)
vj, vk = myhf.get_jk(mol, dm)  # class method
vj, vk = scf.hf.get_jk(mol, dm)  # module function
```
Note that some module functions may require the class object as the first argument,
```python
e_hf = myhf.kernel(conv_tol=1e-5)  # class method
e_hf = scf.hf.kernel(mymf, conv_tol=1e-5)  # module function
```

In PySCF, most functions and classes are **pure**, which means that no
intermediate status is held within the classes, and the arguments of the
methods and functions are immutable during calculations.  Pure functions can be
called any number of times in arbitrary order and their return values should
always be the same.

:::{warning}

 Exceptions to "pure" function behavior are often indicated with an underscore at the
end of the function name,
```python
mcscf.state_average_(mc) 
# the attributes of the mc object may be changed 
# or overwritten by state_average_
```
Be careful when you see functions or methods ending with an underscore!

:::

## Global configurations

Default behaviors in PySCF can be controlled by using global configurations.
A global configuration file is a Python script that contains PySCF configurations.
When PySCF is imported in a Python program (or Python interpreter), the
package will preload the global configuration file to set default values.
For example, the configuration file below detects the available memory in the
operating system at runtime and sets the maximum memory for PySCF,
```{code-block} python
:caption: ~/.pyscf_conf.py

import psutil
MAX_MEMORY = int(psutil.virtual_memory().available / 1e6)
```

By setting `MAX_MEMORY` in the global configuration file, you don't need to set
the `max_memory` attribute in every script. The dynamically determined
`MAX_MEMORY` will be loaded during the program initialization automatically.

There are two ways to identify a global configuration file.
The first is to create a configuration file `.pyscf_conf.py` in your home directory or
in the current working directory.  The second is to set the environment variable
`PYSCF_CONFIG_FILE` to the configuration file (absolute) path.
The environment variable `PYSCF_CONFIG_FILE` has higher priority than
the configuration file found in the home or working directories.
If the environment variable `PYSCF_CONFIG_FILE` is available, PySCF will
use its configurations. If `PYSCF_CONFIG_FILE` is not set or the file it points 
to does not exist, PySCF will look for the file `.pyscf_conf.py` in the
home and working directories.  If no configuration file is found, PySCF 
will use the built-in configurations which are generally conservative.

Global configurations are set in the `pyscf.__config__` module, which
is then imported and used by PySCF,
```python
from pyscf import __config__
MAX_MEMORY = getattr(__config__, 'MAX_MEMORY')
```

Available configurations can be found by reading the source code of PySCF
and its modules. For example, generic configuration parameters include `DEBUG`, `MAX_MEMORY`,
`TMPDIR`, `ARGPARSE`, `VERBOSE`, and `UNIT`, and specific configuration parameters
for a Hartree-Fock calculation can be found at the top of the file,
```{code-block} python
:caption: pyscf/scf/hf.py

from pyscf import __config__

WITH_META_LOWDIN = getattr(__config__, 'scf_analyze_with_meta_lowdin', True)
PRE_ORTH_METHOD = getattr(__config__, 'scf_analyze_pre_orth_method', 'ANO')
MO_BASE = getattr(__config__, 'MO_BASE', 1)
TIGHT_GRAD_CONV_TOL = getattr(__config__, 'scf_hf_kernel_tight_grad_conv_tol', True)
MUTE_CHKFILE = getattr(__config__, 'scf_hf_SCF_mute_chkfile', False)
```
For example, you can choose to change the default behavior associated with the use of
meta Lowdin population analysis,
```{code-block} python
:caption: ~/.pyscf_conf.py

scf_analyze_with_meta_lowdin = False
```

