<div align="left">
  <img src="https://github.com/pyscf/pyscf-doc/blob/master/logo/pyscf-logo.png" height="80px"/>
</div>

PySCF documentation
===================

Installation
------------

* Prerequisites
    - [sphinx](https://pypi.org/project/Sphinx/)
    
    - [sphinxcontrib-bibtex](https://pypi.org/project/sphinxcontrib-bibtex/)

    - Set `PYTHONPATH` to include the PySCF source directory; otherwise, uncomment `sys.path.append(os.path.abspath('path_to_pyscf'))` in [source/conf.py](source/conf.py).

* Make HTML pages

        make html

    The resulting html files can be found in the \"build/html\" directory.

How to contribute
-----------------

1.  Add a rst file \"your\_method.rst\" in the [source/theory](source/theory/) directory in which one describes the basic theory and usage of the method.
2.  Reference \"theory/your\_method.rst\" in the \"toctree\" section in [source/theory.rst](source/theory.rst).
3.  Add a rst file \"your\_module.rst\" in the [source/modules](source/modules/) directory in which one lists the examples and the member classes and functions of the module (this is done by autodoc). (In the \"\_\_init\_\_.py\" file of each module, one should include a simple usage section. See [pyscf.dft.\_\_init\_\_.py](https://github.com/pyscf/pyscf/blob/master/pyscf/dft/__init__.py) as an example.)
4.  Reference \"your\_module.rst\" in the \"toctree\" section in [source/modules.rst](source/modules.rst).
5.  Optionally, one could also add a rst file \"your\_method\_develop.rst\" in the  [source/develop](source/develop/) directory in which one gives guidelines for further development of the module. Also reference \"your\_method\_develop.rst\" in [source/develop.rst](source/develop.rst). 
