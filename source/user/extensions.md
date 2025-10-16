# Extensions 

Due to its flexibility and modular nature, many extensions have been developed
based on PySCF.  As a courtesy, repositories for some of these extensions are
hosted by the [PySCF organization's GitHub](https://github.com/pyscf). However,
the PySCF organization provides no support for these extensions, nor a
guarantee of future compatibility. Questions should be addressed to the authors
of the extension, and issues can be raised on the GitHub repository of the
extension.

Once an extension module has been correctly installed (see [How to install
extensions](#how-to-install-extensions)), you can use it as a regular submodule
of PySCF, e.g.,
```python
import pyscf
from pyscf.semiempirical import MINDO3
mol = pyscf.M(atom='N 0 0 0; N 0 0 1')
MINDO3(mol).run()
```

 A list of PySCF extensions is given below.

| Project       | URL                                                                              |
| ------------- | -------------------------------------------------------------------------------- |
| cas_ac0       | [github.com/CQCL/pyscf-ac0](https://github.com/CQCL/pyscf-ac0)           |
| cornell-shci  | [github.com/pyscf/cornell-shci](https://github.com/pyscf/cornell-shci)   |
| ccpy          | [github.com/piecuch-group/ccpy](https://github.com/piecuch-group/ccpy)   |
| cppe          | [github.com/maxscheurer/cppe](https://github.com/maxscheurer/cppe)       |
| dftd3         | [github.com/pyscf/dftd3](https://github.com/pyscf/dftd3)                 |
| dispersion    | [github.com/pyscf/dispersion](https://github.com/pyscf/dispersion)       |
| dmrgscf       | [github.com/pyscf/dmrgscf](https://github.com/pyscf/dmrgscf)             |
| doci          | [github.com/pyscf/doci](https://github.com/pyscf/doci)                   |
| fciqmc        | [github.com/pyscf/fciqmc](https://github.com/pyscf/fciqmc)               |
| forge         | [github.com/pyscf/pyscf-forge](https://github.com/pyscf/pyscf-forge)     |
| icmpspt       | [github.com/pyscf/icmpspt](https://github.com/pyscf/icmpspt)             |
| mbd           | [github.com/pyscf/mbd](https://github.com/pyscf/mbd)                     |
| naive-hci     | [github.com/pyscf/naive-hci](https://github.com/pyscf/naive-hci)         |
| nao           | [github.com/pyscf/nao](https://github.com/pyscf/nao)                     |
| properties    | [github.com/pyscf/properties](https://github.com/pyscf/properties)       |
| pyqmc         | [github.com/WagnerGroup/pyqmc](https://github.com/WagnerGroup/pyqmc)     |
| qsdopt        | [github.com/pyscf/qsdopt](https://github.com/pyscf/qsdopt)               |
| rt            | [github.com/pyscf/rt](https://github.com/pyscf/rt)                       |
| semiempirical | [github.com/pyscf/semiempirical](https://github.com/pyscf/semiempirical) |
| shciscf       | [github.com/pyscf/shciscf](https://github.com/pyscf/shciscf)             |
| zquatev       | [github.com/sunqm/zquatev](https://github.com/sunqm/zquatev)             |
| tblis         | [github.com/pyscf/pyscf-tblis](https://github.com/pyscf/pyscf-tblis)     |

## How to install extensions

Many extension modules (e.g., the semiempirical module) can be
installed using pip's extra dependency mechanism,
```bash
pip install pyscf[semiempirical]
```

Although not recommended, all extension modules can be installed,
```bash
pip install pyscf[all]
```

## How to install extensions (advanced)

Based on the technique of namespace packages specified in [PEP
420](https://www.python.org/dev/peps/pep-0420/>), PySCF has developed a couple
of methods to install the extension modules.

### Pip command

For pip version newer than 19.0, projects hosted on GitHub can be installed on
the command line,
```bash
pip install git+https://github.com/pyscf/semiempirical
```

A particular release on GitHub can be installed with the release URL you can
look up on GitHub,
```bash
pip install https://github.com/pyscf/semiempirical/archive/v0.1.0.tar.gz
```

### Pip command for local paths

If you wish to load an extension module developed in a local directory, you can
use the local install mode of pip. Use of a Python virtual environment is
recommended to avoid polluting the system default Python runtime environment,
for example,
```bash
python -m venv /home/abc/pyscf-local-env
source /home/abc/pyscf-local-env/bin/activate
git clone https://github.com/pyscf/semiempirical /home/abc/semiempirical
pip install -e /home/abc/semiempirical
```

### Environment variable `PYSCF_EXT_PATH`

You can place the location of each extension module (or a file that contains
these locations) in this environment variable. The PySCF library will parse the
paths defined in this environment variable, and load the relevant submodules.
For example,
```bash
git clone https://github.com/pyscf/semiempirical /home/abc/semiempirical
git clone https://github.com/pyscf/doci /home/abc/doci
git clone https://github.com/pyscf/dftd3 /home/abc/dftd3
echo /home/abc/doci > /home/abc/.pyscf_ext_modules
echo /home/abc/dftd3 >> /home/abc/.pyscf_ext_modules
export PYSCF_EXT_PATH=/home/abc/semiempirical:/home/abc/.pyscf_ext_modules
```

Using this definition of `PYSCF_EXT_PATH`, the three extension
submodules (semiempirical, doci, dftd3) are loaded when PySCF is
imported, and you don't have to use a Python virtual environment.


## Examples

### Installing DMRG solvers

The Density Matrix Renormalization Group (DMRG) theory is a method for solving ab initio quantum chemistry problems.
PySCF can be used with three implementations of DMRG: [Block](https://sanshar.github.io/Block), [block2](https://block2.readthedocs.io/en/latest), and [CheMPS2](http://sebwouters.github.io/CheMPS2/index.html).

Installing [Block](https://sanshar.github.io/Block/build.html) requires a C++11
compiler. Alternatively, you can download the precompiled Block binary from https://sanshar.github.io/Block/build.html .

`block2` can be easily installed via `pip install block2` or `pip install block2-mpi`.
To building block2 from source, please refer to the installation guide https://block2.readthedocs.io/en/latest/user/installation.html .

Before using Block or CheMPS2, you need create a configuration file
`pyscf/dmrgscf/settings.py` (as shown by the `settings.py.example`) to
store the path where the DMRG solver was installed.

### Installing TBLIS

[TBLIS](https://github.com/devinamatthews/tblis) provides a native algorithm for
performing tensor contraction for arbitrarily high-dimensional tensors. Unlike
the implementation in `numpy.einsum`, which transforms tensors into matrices by
permuting the indices, invokes BLAS for matrix contractions, and then
back-permutes the results, TBLIS eliminates the need for tensor transpositions
and data movements. TBLIS eliminates the need to transpose tensors and data
movements. This results in improved performance in various correlated quantum
chemistry methods in PySCF, such as the coupled cluster methods.

The `tblis-einsum` extension provides an interface to TBLIS, which offers an
efficient implementation for tensor contractions in the style of `numpy.einsum`.
The `tblis-einsum` plugin can be installed using the command:
```bash
$ pip install pyscf-tblis
```
