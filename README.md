<div align="left">
  <img src="https://github.com/pyscf/pyscf-doc/blob/master/logo/pyscf-logo.png" height="80px"/>
</div>

PySCF documentation
===================


Installation
------------

* Prerequisites
    - sphinx-build

* Make HTML

        make html


How to contribute
-----------------

1.  Add a rst file *\"your\_method.rst\"* in the [source/theory](source/theory/) directory in which you describe the theory of your method.
2.  Insert *\"theory/your\_method.rst\"* in the \"toctree section\" in [source/theory/theory.rst](source/theory/theory.rst).
3.  Add a rst file *\"your\_module.rst\"* in the [source/](source/) directory in which you list the example files and the classes and functions in your module.
4.  Reference *\"your\_module.rst\"* properly in [source/index.rst](source/index.rst).
