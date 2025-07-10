# Overview 

PySCF-forge is a staging ground for new methods or features that are expected
to be made available in the main PySCF distribution. Not all new features are
appropriate for PySCF-forge; experimental or niche features are more
appropriately distributed as [](./extensions).

## How to install PySCF-forge

PySCF-forge is distributed as a 
[separate GitHub repository](https://github.com/pyscf/pyscf-forge).
Like PySCF, it can be installed via pip,
```
pip install pyscf-forge
# or
pip install git+https://github.com/pyscf/pyscf-forge
```
or cloned and built from source.  As long as PySCF-forge can be found in your
`PYTHONPATH`, its modules can be imported as if they were native to PySCF.
Alternatively, you can define the `PYSCF_EXT_PATH` environment variable to
point to your local `pyscf-forge` directory,
```
export PYSCF_EXT_PATH=/path/to/pyscf-forge
```
PySCF will read the `PYSCF_EXT_PATH` environment variable and load modules from
this path at runtime. For more details of the `PYSCF_EXT_PATH` environment
variable, see the
[Extensions](./extensions.md#environment-variable-pyscf-ext-path) section.  For
example,
```python
from pyscf import gto, mcpdft  # mcpdft is in pyscf-forge
mol = gto.M(
    atom = 'O 0 0 0; O 0 0 1.2',
    basis = 'ccpvdz',
    spin = 2)
mf = mol.RHF().run ()
# Calling the MC-PDFT method provided by the pyscf-forge modules
mc = mcpdft.CASCI(mf, 'tPBE', 6, 8).run()
```

## How to develop in PySCF-forge 

If you are developing new features or modifying the code in PySCF-forge, an
editable installation is recommended.  By configuring the package in editable
mode, you can modify existing modules and add new features to PySCF-forge.
After cloning the library to your local repository, there are two ways to
enable the editable installation.

### Method 1. Using pip for editable installation

Install the package with the following pip command,
```
pip install --no-deps -e /path/to/pyscf-forge
```
This command creates a `.pth` file in `~/.local/lib/python3.*/site-packages/`
or other Python runtime paths. It is recommended to use this method with Python
virtual environments.

### Method 2. Setting an environment variable

As described above, you can also define the `PYSCF_EXT_PATH` environment
variable to point to your local PySCF-forge directory,
```
export PYSCF_EXT_PATH=/path/to/pyscf-forge
```

### Adding new features: An example

Suppose you want to create a new module in PySCF-forge that provides a new
feature called `abc`.  You can follow these steps to add the module,

1. Install PySCF-forge in editable installation mode.

2. Create a folder named `pyscf-forge/pyscf/abc`.  Due to the editable
installation mode, this folder can be readily imported as a regular `pyscf`
module.
```
>>> from pyscf import abc
>>> print(abc)
<module 'pyscf.abc' (namespace)>

>>> print(abc.__path__)
_NamespacePath(['path/to/pyscf-forge/pyscf/abc'])
```

3. Add Python code files to the `pyscf-forge/pyscf/abc` directory.  This
process is similar to developing new methods for PySCF.  For example, you
can add the following Python files into the `pyscf-forge/pyscf/abc` folder,
```
pyscf-forge
├── ...
└── pyscf
    ├── ...
    └── abc
        ├── __init__.py
        ├── rabc.py 
        └── uabc.py
```

### Path Conflicts

What if you want to add a feature to a module that already exists in PySCF?
For example, if you want to add a new method, like `pp_rpa.py`, to the
`pyscf/tdscf` folder, this could conflict with the existing `pyscf.tdscf`
module in the PySCF core repository.  Adding features to existing modules
requires more a slightly more complex configuration.

To import the `pp_rpa` module from the PySCF-forge repository, you will need to
modify the `__init__.py` file of the PySCF core module.  Add the following line
of code to modify the `__path__` attribute of the `pyscf.tdscf` module,
```python
__path__ = __import__('pkgutil').extend_path(__path__, __name__)
```
This command extends the search path of the `tdscf` module, resulting in the
`__path__` attribute being set to
```
['/path/to/pyscf/pyscf/tdscf',
 '/path/to/pyscf-forge/pyscf/tdscf']
```
This configuration allows Python to locate and load the new `pp_rpa.py` module
from the extended directory in PySCF-forge.

Note that the PySCF core `tdscf` folder already contains an `__init__.py` file.
To avoid overwriting the existing `__init__.py` file in PySCF during the
installation of PySCF-forge, you should not add an `__init__.py` file in the
`pyscf-forge/pyscf/tdscf` directory.

The structure of PySCF and PySCF-forge can be organized as follows,
```
pyscf
├── ...
└── pyscf
    ├── ...
    └── tdscf
        ├── __init__.py  // modify the __path__ attribute in pyscf core module
        ├── rhf.py
        ├── rks.py
        └── ...

pyscf-forge
├── ...
└── pyscf
    ├── ...
    └── tdscf  // no __init__.py file in pyscf-forge
        └── pp_rpa.py
```

When installing the PySCF-forge wheels using `pip install` in the normal
installation mode, the `pp_rpa.py` file will be added to the `pyscf/tdscf`
folder, integrating seamlessly as part of the regular PySCF module.  After this
standard installation, there is no need to adjust the `__path__` attribute, as
all features and modules are located within the same directory.
