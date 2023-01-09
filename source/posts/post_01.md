---
blogpost: true
date: February 1, 2021
author: James Smith
location: World
category: Tutorial
tags: HF, DFT, MCSCF
language: English
---

# Geometry Opt in PySCF

Here's how we can run a CASSCF geometry optimization of N<sub>2</sub>.

```python
from pyscf import gto, scf, mcscf
from pyscf.geomopt.geometric_solver import optimize

# Single Point Energy Calc.
mol = gto.M(atom='N 0 0 0; N 0 0 1.2', basis='ccpvdz')
mf = scf.RHF(mol)
mc = mcscf.CASSCF(mf, 4, 4)
mc.kernel()

# Set up geometry optimization
conv_params = {
    'convergence_energy': 1e-4,  # Eh
    'convergence_grms': 3e-3,    # Eh/Bohr
    'convergence_gmax': 4.5e-3,  # Eh/Bohr
    'convergence_drms': 1.2e-2,  # Angstrom
    'convergence_dmax': 1.8e-2,  # Angstrom
}

mol_eq = mc.Gradients().optimizer(solver='geomeTRIC').kernel(conv_params)
```

## Here's how we can include an image

We have to use relative paths if our post is a `Markdown` file. It's not the end of the world, it's just something to be aware of.

Here's the `markdown` code to include the image below (pretty easy right?):
```markdown
![](../images/post_02/fe_pdi_lowres.png)
```

![](../images/post_02/fe_pdi_lowres.png)


## Embed Using iframe and HTML

Here's the HTML code to imbed the a molecule in the post:

```html
<iframe style="width: 800px; height: 600px;" frameborder="0" src="https://embed.molview.org/v1/?mode=balls&bg=white"></iframe>
```

<iframe style="width: 800px; height: 600px;" frameborder="0" src="https://embed.molview.org/v1/?mode=balls&bg=white"></iframe>
