.. include:: ../links.inc

.. _manual-install:
.. _standard-instructions:

Install via :code:`pip` or :code:`conda`
========================================

.. hint::
   If you're unfamiliar with Python, we recommend using our :ref:`installers`
   instead.

AVNI version |version| requires Python version |min_python_version| or higher. If you
need to install Python, please see :ref:`install-python`.

Installing AVNI with all dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We suggest to install AVNI into its own ``conda`` environment.

The dependency stack is large and may take a long time (several tens of
minutes) to resolve on some systems via the default ``conda`` solver. We
therefore highly recommend using `mamba <https://mamba.readthedocs.io/>`__
instead, a ``conda`` replacement that is **much** faster.

Run in your terminal:

.. code-block:: console

    $ conda install --channel=conda-forge --name=base mamba
    $ mamba create --override-channels --channel=conda-forge --name=avni avni

This will create a new ``conda`` environment called ``avni`` (you can adjust
this by passing a different name via ``--name``) and install all
dependencies into it.


Installing a minimal AVNI with core functionality only
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you only need AVNI's core functionality including 2D plotting (but
**without 3D visualization**), install via :code:`pip`:

.. code-block:: console

   $ pip install avni

or via :code:`conda`:

.. code-block:: console

   $ conda create --strict-channel-priority --channel=conda-forge --name=avni avni-base

This will create a new ``conda`` environment called ``avni`` (you can adjust
this by passing a different name via ``--name``).

Installing AVNI for other scenarios
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`advanced_setup` page has additional
tips and tricks for special situations (servers, notebooks, CUDA, installing
the development version, etc). The :ref:`contributing` has additional
installation instructions for (future) contributors to AVNI (e.g, extra
dependencies for running our tests and building our documentation).

Python IDEs
===========

Most users find it convenient to write and run their code in an `Integrated
Development Environment`_ (IDE). Some popular choices for scientific
Python development are:

- `Visual Studio Code`_ (often shortened to "VS Code" or "vscode") is a
  development-focused text editor that supports many programming languages in
  addition to Python, includes an integrated terminal console, and has a rich
  ecosystem of packages to extend its capabilities. Installing
  `Microsoft's Python Extension
  <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`__ is
  enough to get most Python users up and running. VS Code is free and
  open-source.

- `PyCharm`_ is an IDE specifically for Python development that provides an
  all-in-one installation (no extension packages needed). PyCharm comes in a
  free "community" edition and a paid "professional" edition, and is
  closed-source.
