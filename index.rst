.. PeleMP documentation master file, created by
   sphinx-quickstart on Fri Sep  2 13:09:44 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PeleMP's documentation!
==================================

`PeleMP` is repository of the multiphysics capabilities for use with other `Pele` codes:

* `PeleC <https://amrex-combustion.github.io/PeleC/>`_
* `PeleLM <https://amrex-combustion.github.io/PeleLM/>`_
* `PeleLMeX <https://amrex-combustion.github.io/PeleLMeX/>`_

The official project `homepage <https://amrex-combustion.github.io/>`_, and can be obtained via `GitHub <https://github.com/AMReX-Combustion/PeleMP>`_.
The documentation pages appearing here are distributed with the code in the ``Docs`` folder as "restructured text" files.
The html is built automatically with certain pushes to the `PeleMP` GibHub repository.
A local version can also be built as follows ::

  cd ${PELEMP_HOME}/Docs
  cmake CMakeList.txt
  make sphinx

where ``PELEMP_HOME`` is the location of your clone of the `PeleMP` repository. To view the local pages, point your web browser at the file ``${PELEMP_HOME}/Docs/sphinx/html/index.html``.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   GettingStarted.rst
   Models.rst


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
