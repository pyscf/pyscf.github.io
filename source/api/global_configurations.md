# Global configurations

Global configurations set defaults shared by PySCF calculations, such as
memory, temporary-file locations, and method settings. PySCF loads them through
`pyscf.__config__` when the package is imported.

## Configuration files

A PySCF configuration file is a Python script that assigns values to
configuration variables. PySCF searches for a file in this order:

1. The path given by the `PYSCF_CONFIG_FILE` environment variable.
2. `.pyscf_conf.py` in the current working directory.
3. `.pyscf_conf.py` in the home directory.

Only the first existing file is loaded; the files are not merged. If
`PYSCF_CONFIG_FILE` points to a missing file, PySCF continues searching the
working and home directories. If no file is found, it uses the built-in and
environment-derived defaults.

For example, a `.pyscf_conf.py` can contain:

```python
MAX_MEMORY = 4000
VERBOSE = 4
scf_hf_SCF_conv_tol = 1e-10
```

These assignments become attributes of `pyscf.__config__`.

### Method-specific defaults

Individual modules read settings from `pyscf.__config__`. For
example, `scf_hf_SCF_conv_tol` sets the default `SCF.conv_tol`, and
`gto_mole_Mole_cart` sets the default `Mole.cart`.

Object attributes do not automatically have corresponding global configuration
names. Look for the configuration variable names in the relevant module's source
code, for example:

```python
from pyscf import __config__
conv_tol = getattr(__config__, 'scf_hf_SCF_conv_tol', 1e-9)
```

Here, `scf_hf_SCF_conv_tol` is the configurable parameter.

### When changes take effect

PySCF normally loads the configuration once, when the package is first
imported. Many modules then copy its settings into module constants, class
defaults, or function defaults.

The variables provided by `pyscf.__config__` should be treated as immutable
values during runtime. Changing its attributes at runtime does not reliably
update settings already read by other modules. The configuration file should be
prepared before importing PySCF.


## Environment variables and shared defaults

The following environment variables set initial defaults before the
configuration file is executed. Assignments in that file can override them.

| Environment variable | Configuration name | Default and meaning |
| --- | --- | --- |
| `PYSCF_MAX_MEMORY` | `MAX_MEMORY` | 4000 MB; default calculation memory budget |
| `PYSCF_TMPDIR` | `TMPDIR` | System temporary directory; location for scratch files |
| `PYSCF_ARGPARSE` | `ARGPARSE` | `False`; enable command-line argument parsing by the molecule builder |
| `PYSCF_DISABLE_EVAL` | `DISABLE_EVAL` | `False`; disable expression evaluation in supported input parsers |

`PYSCF_MAX_MEMORY` is not a hard limit on process memory. A calculation can
exceed it, particularly when in-core operation is forced.


## Disabling Python `eval` in input parsers

Some geometry, basis-set, ECP and pseduo-potential parsers accept Python
expressions in numeric fields and use `eval()` to evaluate them. This allows
expressions in input strings alongside ordinary numbers.

Set `PYSCF_DISABLE_EVAL=1` before importing PySCF, or put `DISABLE_EVAL = True`
in the configuration file, to disable this behavior. By setting this, the input
string must be numeric values instead of expressions.
