<div align="left">
  <img src="https://github.com/pyscf/pyscf-doc/blob/master/logo/pyscf-logo.png" height="80px"/>
</div>

PySCF documentation
===================


Installation
------------

* Prerequisites
    - sphinx-build

    - Set `PYTHONPATH` to include PySCF source directory; otherwise, edit [source/conf.py](source/conf.py) by adding `sys.path.insert(0, os.path.abspath('path_to_pyscf_source'))`.

* Make HTML pages

        make html

    The resulting html files can be found in the \"build/html\" directory.

How to contribute
-----------------

1.  Add a rst file *\"your\_method.rst\"* in the [source/theory](source/theory/) directory in which you describe the theory of your method.
2.  Insert *\"theory/your\_method.rst\"* in the \"toctree\" section in [source/theory.rst](source/theory.rst).
3.  Add a rst file *\"your\_module.rst\"* in the [source/modules](source/modules/) directory in which you list the example files and the classes and functions in your module. (In the *\"\_\_init\_\_.py\"* file of each module, one should include a simple usage section. See *\"pyscf/dft/\_\_init\_\_.py\"* as an example.)
4.  Reference *\"your\_module.rst\"* properly in [source/modules.rst](source/modules.rst).
