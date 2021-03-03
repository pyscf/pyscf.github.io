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
cd sources
ablog build
ablog serve
```

If you want to regenerate the API Docs, run the following:

```bash
cd source
../scripts/generate_api_docs.sh
```

> :warning: Warning: PySCF must be accessible in your current Python environment when you run `../scripts/generate_api_docs.sh`.

## How to push changes (and make them viewable online)

```bash

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