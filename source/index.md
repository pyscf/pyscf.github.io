---
html_theme.sidebar_secondary.remove: true
---

# Quantum chemistry with Python 

<div style="background-color: #013243; padding: 5px; border-radius: 2px; font-family: monospace; font-size: 20px; color: white; display: flex; justify-content: center; align
items: center; height: 60px; text-align: center; margin-bottom: 20px;">

  <p style="margin: 0;">
    <a href="https://github.com/pyscf/pyscf/issues/3016">3rd PySCF Developers Meeting</a> is scheduled for **Aug 27, 2026, in University of Chicago**!
  </p>

</div>

```{gallery-grid}
:grid-columns: 1 2 2 3

- header: "{fab}`python;pst-color-primary` Primarily Python"
  content: "Easy to read and write, for quick development with fewer mistakes."
- header: "{fas}`bolt;pst-color-primary` Fast"
  content: "Uses NumPy, SciPy, and custom C/C++ code."
- header: "{fa}`code-branch;pst-color-primary` Free and open-source"
  content: "Accessible and community-driven. Distributed on [GitHub](https://github.com/pyscf/pyscf) under the [Apache-2.0](https://github.com/pyscf/pyscf/tree/master?tab=Apache-2.0-1-ov-file#readme) license."
- header: "{fas}`gears;pst-color-primary` Modular"
  content: "Easily integrated into complex workflows or other software packages."
- header: "{fa}`list-ul;pst-color-primary` Comprehensive"
  content: "Extensive collection of electronic structure methods, for molecules and periodic solids."
- header: "{fas}`lightbulb;pst-color-primary` Extensible"
  content: "See the list of [Extensions](user/extensions) to PySCF."
```

<!--- 
The Python-based Simulations of Chemistry Framework (PySCF) is an open-source
collection of electronic structure modules powered by Python. The package
provides a simple, lightweight, and efficient platform for quantum chemistry
calculations and methodology development.
PySCF can be used to simulate the properties of molecules, crystals, and
custom Hamiltonians using mean-field and post-mean-field methods.
To ensure ease of extensibility, almost all of the features in PySCF are
implemented in Python, while computationally critical parts are
implemented and optimized in C. Using this combined Python/C
implementation, the package is as efficient as the best existing C or Fortran
based quantum chemistry programs.
In addition to its core libraries, PySCF supports a rich
ecosystem of :ref:`installing_extproj`.
--->


```{toctree}
:hidden:

quickstart
user/index
contributor/index
About <about>
```
