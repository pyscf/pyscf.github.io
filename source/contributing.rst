How to contribute
*****************

Standards for core modules
==========================
In PySCF-2.0, features are split into core modules and extension projects.
Extension projects are hosted at https://github.com/pyscf.
The standards for features in the core modules are

* **Documentation**.
  Core modules need to be documented. Documentation should include the
  user-guide, the API manual and a few runnable examples. The document repo is
  https://github.com/pyscf/pyscf-doc.

* **Unittests and test coverage**.
  Unittests should be created for all APIs of the core modules.

* **Minimal compilation requirements**.
  PySCF core modules are compiled and released in a general Python package.
  These modules need to be compatible with as many platforms and system
  environments as possible. If your code requires specific Python versions,
  specific compiler sets, or specific libraries (e.g. MPI), it is recommended to
  develop the features in the extension project.


How to make an extension
========================

1. Download the template (https://github.com/pyscf/extension-template), run the
   script in the extension template::

    $ curl -L https://github.com/pyscf/extension-template/archive/refs/tags/0.1.0.tar.gz | tar xvzf -
    $ mv extension-template new-project
    $ cd new-project
    $ make new-project

2. Check whether files like `README.md`, `setup.py`, `__init__.py` are correctly updated.

3. In the new project folder, initialize git repo::

    $ git init
    $ git add .
    $ git commit -m "An example of new extension package"

4. If your program needs to compile C/C++ code, you can edit the setup.py to
   include your code files. For example::

    SO_EXTENSIONS = {
      'pyscf.new_feature.new_feature': ['pyscf/new_feature/new_feature.cc']
    }

5. Release the extension project on https://github.com/pyscf.

   * Open an issue on https://github.com/pyscf/pyscf about the new project.
   * When the issue is approved, a repo will be created under https://github.com/pyscf/new-project.
     You can run the commands below to sync your local repo to the remote one::

      $ git remote add origin https://github.com/pyscf/new_project
      $ git push

   * Open a PR on https://github.com/pyscf/pyscf to update the pyscf setup.py.
     In the PR, you can put the new extension project in the field `extras_requires`.


How to release a new version
============================
Since PySCF-2.0, a PR-review-release pipeline was created to release the program
on pypi, conda and docker. Here are the steps to release a new version:

1. Update the `__version__` variable in the file `pyscf/__init__.py`

2. In the file CHANGELOG, shortly summarize the new features, the bugfixes or
   other relevant changes

3. Open a PR on https://github.com/pyscf/pyscf. When the PR is merged, the
   release pipeline will be automatically triggered. In the action page
   https://github.com/pyscf/pyscf/actions/workflows/publish.yml, check whether
   all jobs are correctly executed.
