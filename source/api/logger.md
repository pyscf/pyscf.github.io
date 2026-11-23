# Logging API

`pyscf.lib.logger` provides the output and timing utilities used by PySCF
calculations. It is independent of Python's built-in `logging` module.

Use this module when adding diagnostics to a PySCF method so that messages
follow the calculation's output settings.

## Using a logger

By default, a new logger uses `sys.stdout` and verbosity `NOTE` (3) for the
output stream and verbosity level.

You can write messages through a `Logger` object:

```python
from pyscf.lib import logger
log = logger.new_logger(mol)
log.info('SCF energy = %.12f', mf.e_tot)
```

Or pass a method object directly to a logging function:

```python
logger.info(mf, 'SCF energy = %.12f', mf.e_tot)
```

When an ouptut file is specified for `mol.output`, `Mole.build()` opens it and
assigns the stream to `mol.stdout`. The logger writes to this stream but does
not automatically close it. To avoid an unclosed-file warning, you need to
explicitly close it with `mol.stdout.close()` after finishing all calculations.


## Verbosity levels

Larger values produce more output. A message at level `L` is normally shown
when `verbose >= L`. These numbers are specific to PySCF and differ from
Python's logging levels.

| Level | Value | Intended use |
| --- | --- | --- |
| `QUIET` | 0 | Suppress ordinary output |
| `ERROR` | 1 | Error diagnostics |
| `WARN` | 2 | Warnings |
| `NOTE` | 3 | Important results and summaries |
| `INFO` | 4 | Routine calculation progress |
| `DEBUG` | 5 | Diagnostic details and default timing output |
| `DEBUG1`–`DEBUG4` | 6–9 | Increasingly detailed diagnostics |

The corresponding calls are `error()`, `warn()`, `note()`, `info()`,
`debug()`, and `debug1()` through `debug4()`.
`error()` always writes to the current `sys.stderr`, even at `QUIET` level.

## Timer

Timers report durations in seconds.

Timer functions (`init_timer()`, `timer()`, and `timer_debug1()`) return
(`cpu_time`, `wall_time`). Passing both to timer functions can report CPU and
wall time:

```python
t0 = log.init_timer()
# Run the calculation to be timed.
t1 = log.timer('calculation', *t0)
```

The returned pair contains new timestamps, not the reported durations.

Calling `log.timer(label)` without a starting timestamp measures CPU time since
the previous timer call.

By default, timing messages appear at verbosity `DEBUG` (5) or higher. The
module-level setting `TIMER_LEVEL` controls this threshold. To show timing
messages at different verbosity level, e.g., at `INFO` (4), set
`logger.TIMER_LEVEL = logger.INFO`.

`timer_debug1()` outputs message at verbosity `DEBUG1` (6) or higher. It cannot
be configured.

Both timer functions update their timestamps and return them even when no
message is printed, so verbosity changes do not alter the timing sequence.
