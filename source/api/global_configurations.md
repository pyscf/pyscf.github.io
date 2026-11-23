# Global configurations

Global settings define defaults shared by PySCF calculations, such as
memory limits, temporary-file locations, and method settings. PySCF loads them through
`pyscf.__config__` when the package is imported.

## Configuration files

A PySCF configuration file is a Python script that assigns values to
configuration variables. PySCF searches for a file in this order:

1. The path given by the `PYSCF_CONFIG_FILE` environment variable.
2. `.pyscf_conf.py` in the current working directory.
3. `.pyscf_conf.py` in the home directory.

PySCF loads only the first file it finds; it does not combine the files. If
`PYSCF_CONFIG_FILE` points to a missing file, PySCF continues searching the
working and home directories. If no file is found, it uses the built-in
defaults and any settings supplied through environment variables.

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

Not every object attribute has a matching global setting. To find the name
of a supported setting, look in the relevant module's source code. For example:

```python
from pyscf import __config__
conv_tol = getattr(__config__, 'scf_hf_SCF_conv_tol', 1e-9)
```

Here, `scf_hf_SCF_conv_tol` is the name to use in the configuration file.

### When changes take effect

PySCF normally loads the configuration once, when the package is first
imported. Many modules keep their own copies of these settings as defaults.

Prepare the configuration file before importing PySCF, and avoid changing
`pyscf.__config__` while the program is running. Other modules may have
already read the settings, so later changes may not take effect.


## Environment variables and defaults

The following environment variables set defaults before PySCF loads the
configuration file. Settings in that file take precedence.

| Environment variable | Configuration name | Default and meaning |
| --- | --- | --- |
| `PYSCF_MAX_MEMORY` | `MAX_MEMORY` | 4000 MB; default calculation memory budget |
| `PYSCF_TMPDIR` | `TMPDIR` | System temporary directory; location for scratch files |
| `PYSCF_ARGPARSE` | `ARGPARSE` | `False`; enable command-line argument parsing by the molecule builder |
| `PYSCF_DISABLE_EVAL` | `DISABLE_EVAL` | `False`; disable expression evaluation in supported input parsers |

`PYSCF_MAX_MEMORY` does not strictly limit how much memory a calculation can
use. A calculation may exceed it, especially if it is forced to keep data
in memory rather than on disk.

## Temporary directories

To choose where PySCF stores scratch files, set the environment variable
`PYSCF_TMPDIR` before importing PySCF, or set `TMPDIR` in the PySCF
configuration file. If neither is set, PySCF uses Python's system temporary
directory, which can be selected through the environment variable `TMPDIR`.
The directory must already exist. PySCF does not create it.

Scratch files can contain temporary SCF checkpoints, DIIS vectors,
density-fitting integrals, AO-to-MO transformation intermediates,
and coupled-cluster intermediates.

PySCF removes its temporary files when they are closed or no longer needed.
If a calculation is interrupted, temporary files may be left behind.

The default SCF checkpoint is also temporary. It remains available after
`kernel()` returns, but is normally deleted when the SCF object is no longer
in use. To keep the checkpoint, assign a path before running the calculation:

```python
mf.chkfile = 'restart.chk'
mf.kernel()
```


## Disabling Python `eval` in input parsers

Some geometry, basis-set, ECP and pseudopotential parsers accept Python
expressions in numeric fields and use `eval()` to evaluate them. This allows
expressions in input strings alongside ordinary numbers.

Set `PYSCF_DISABLE_EVAL=1` before importing PySCF, or put `DISABLE_EVAL = True`
in the configuration file, to disable this behavior. Numeric fields must
then contain numbers rather than Python expressions.
