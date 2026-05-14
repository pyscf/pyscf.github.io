# Dispersion Convention

Dispersion corrections in DFT can be enabled either:

1. implicitly through the `xc` functional string,
2. or explicitly through the `disp` keyword.

The explicit `disp` keyword always takes precedence over the dispersion information encoded in `xc`.

## Enabling Dispersion through `xc`

Some functionals encode the dispersion model directly in the XC name. For examples,
```
mf.xc = 'wb97x-d3bj'
mf.xc = 'wb97x-d4'
```
These automatically enable dispersion corrections.

When calling the dispersion module, the `xc` attribute is parsed and normalized to the internal representation `(disp_method, disp_version, with_3body)`. For example,

|xc           |Parsed dispersion        |
|-------------|-------------------------|
|'b3lyp'      |no dispersion            |
|'wb97x-d3bj' |('wb97x', 'd3bj', False) |
|'wb97x-d4'   |('wb97x', 'd4', True)    |

## Enabling Dispersion through `disp`

Dispersion can also be enabled explicitly through the `disp` attribute:
```
mf.xc = 'b3lyp'
mf.disp = 'd3bj'
```
This corresponds to `('b3lyp', 'd3bj', False)`.

The disp keyword overrides any dispersion setting implied by xc. For example, 
```
mf.xc = 'wb97x-d3bj'
mf.disp = 'd4'
```
This uses D4 dispersion instead of D3(BJ).

## 2-Body vs 3-Body Convention
Internally, the normalized representation separates the dispersion version from
the 3-body flag. By default, D4 includes the 3-body contribution, whereas D3
does not. They can both be adjusted by suffix 2b and atm:
* suffix 2b disables the 3-body term;
* suffix atm enables the ATM 3-body term.

For example,
```
mf.disp = 'd3bj2b'   # D3(BJ), pairwise only
mf.disp = 'd3bjatm'  # D3(BJ) + ATM
```

## Custom Dispersion Parameterization
The dispersion version can be specified using the syntax
```
disp_version:method
```

For example, 
```
mf.disp = 'd4:pbe0'
```
employs DFTD4 dispersion and passes `pbe0` as the DFT functional name to the underlying dispersion library.
In the internal representation (returned by `parse_disp` function), this
setting leads to `('pbe0', 'd4', True)`.
Another example,
```
mf.disp = 'd3bj:pbe'
```
corresponds to `('pbe', 'd3bj', False)`.

## Examples
* No dispersion
```
mf.xc = 'pbe'
mf.disp = None
```

* Implicit dispersion from XC
```
mf.xc = 'wb97x-d3bj'
mf.disp = None
```
This setting enables `d3bj` type of dispersion parameterized for the `wb97x` functional.

* Explicit dispersion
```
mf.xc = 'b3lyp'
mf.disp = 'd3bj'
```
This setting leads to `d3bj` type of dispersion parameterized for the `b3lyp` functional.

* Three-body ATM term
```
mf.xc = 'b3lyp'
mf.disp = 'd3bjatm'
```
This setting provides `d3bj` dispersion with 3-body contribution for the `b3lyp` functional.

* Explicit dispersion version
```
mf.xc = 'scan'
mf.disp = 'd4:pbe'
```
This setting ignores the `xc` attribute and uses `d4` dispersion parameterized
for the `pbe` functional.
￼

