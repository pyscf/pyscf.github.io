# How to Contribute to the PySCF Docs

## Requirements

Pip installable packages

- pyscf
- sphinx
- ablog
- myst-parser
- sphinx-material
- sphinxcontrib-bibtex
- nbsphinx

> Pro-tip: Use `pyscf-doc/requirements.txt` to install all prerequisite packages at once.

## How to build the docs locally
All sphinx related sources files (i.e. `.rst` and `.md`) are contained in `pyscf-doc/source`.

```bash
cd source
../scripts/material_style_patch.sh apply
../scripts/generate_api_docs.sh
ablog build
ablog serve # this will open a tab in your default browser
```

> :warning: :warning: :warning: PySCF must be accessible in your current Python environment when you run `../scripts/generate_api_docs.sh`.

> :warning: :warning: :warning: Running `ablog build` will be slow after you've generated the API docs. There are two hack-y strategies to speed these up: 1) Follow the instructions above and after running `ablog build` you can delete `source/api_docs`. The HTML files will still exist in `source/_website` so they'll still show up when you serve the website, but `sphinx` will no longer need to generate them every time. 2) You can skip the `../scripts/generate_api_docs.sh` command above and deal with the broken link.

## How to push changes

> :warning: :warning: :warning: Run the following before `git add`-ing any files. This is a temporary workaround until the upstream branch of `pyscf-doc` switches to using the `source/conf.py` and `source/index.rst` the are correctly setup for the `sphinx-material` theme.


```bash
cd source
../scripts/material_style_patch.sh revert
```

## Adding Blog Posts

Create a new `.md` file in `pyscf-doc/source/posts` and add the following header (modified for your post):

```
---
blogpost: true
date: February 1, 2021
author: James Smith
location: World
category: Tutorial
tags: HF, DFT, MCSCF
language: English
---
```

If you want to write a post in `.rst` that's fine too! Just use the following in your header:

```
:blogpost: true
:date: Oct 10, 2020
:author: Nabil Freij
:location: World
:category: Manual
:language: English
```