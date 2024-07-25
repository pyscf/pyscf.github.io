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

| Project         |   URL                                       |
| --------------- | ------------------------------------------- |
| cas_ac0         | https://github.com/CQCL/pyscf-ac0           |
| cornell-shci    | https://github.com/pyscf/cornell-shci       |
| dftd3           | https://github.com/pyscf/dftd3              |
| dmrgscf         | https://github.com/pyscf/dmrgscf            |
| doci            | https://github.com/pyscf/doci               |
| fciqmc          | https://github.com/pyscf/fciqmc             |
| icmpspt         | https://github.com/pyscf/icmpspt            |
| mbd             | https://github.com/pyscf/mbd                |
| naive-hci       | https://github.com/pyscf/naive-hci          |
| nao             | https://github.com/pyscf/nao                |
| qsdopt          | https://github.com/pyscf/qsdopt             |
| rt              | https://github.com/pyscf/rt                 |
| semiempirical   | https://github.com/pyscf/semiempirical      |
| shciscf         | https://github.com/pyscf/shciscf            |
| zquatev         | https://github.com/sunqm/zquatev            |
| tblis           | https://github.com/pyscf/pyscf-tblis        |

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

