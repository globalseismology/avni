.. include:: ../links.inc

.. _standard-instructions:
.. _install_ref:

Installation
============

..
   .. hint::
..    If you're unfamiliar with Python, we recommend using our :ref:`installers`
..    instead.

AVNI version |version| requires Python version |min_python_version| or higher.
For the best experience, please considering using Anaconda as a virtual environment and
package manager for Python (see :ref:`install-python`). Follow the instructions below to
install AVNI with Anaconda.

Dependencies
~~~~~~~~~~~~

AVNI is built on top of Fortran and standard Python libraries - as such,
the following projects are required dependencies of AVNI:

* `gfortran Fortran compiler <https://gcc.gnu.org/wiki/GFortran>`_ - Legacy
  Fortran code is interfaced with AVNI using the Fortran compiler.
* `NumPy <https://pypi.org/project/numpy/>`_ - NumPy arrays provide a core
  foundation for AVNI's data array access.
* `pandas <https://pandas.pydata.org>`_ - Data analysis and manipulation routines
  that are required for AVNI's I/O routines.
* `scipy <https://scipy.org>`_ - Scientific routines that perform standard calculations.
* `configobj <https://pypi.org/project/configobj/>`_ - Configuration files
* `xarray <https://docs.xarray.dev/en/stable/>`_ - Labelled multi-dimensional arrays
  and datasets in Python
* `matplotlib <https://pypi.org/project/matplotlib/>`_ - Used for colormaps and
  2D plotting.
* `fortranformat <https://pypi.org/project/fortranformat/>`_ - Some I/O routines
  in AVNI write files in format specifiers that are easily specified in Fortran.
* `pint <https://pypi.org/project/Pint/>`_ - Unit system for conversion of state
  variables from SI to other systems.
* `h5py <https://www.h5py.org>`_ - Pythonic interface to the HDF5 binary data format
  that is used by AVNI is writing files in reference seismic data format (RSDF) files.

.. Optional Dependencies
.. ~~~~~~~~~~~~~~~~~~~~~
.. AVNI includes several optional dependencies for visualization and reading a
.. variety of additional file formats, including:

Installing Fortran compiler
^^^^^^^^^^^^^^^^^^^^^^^^^^^
A prerequisite for AVNI is the `gfortran Fortran compiler <https://gcc.gnu.org/wiki/GFortran>`_.
Once Fortran is installed, please set the environment variable ``F90`` in your
command line shell to ``gfortran``.

Run in your terminal:

.. code-block:: console

    $ echo $SHELL
    $ export F90=gfortran # if SHELL is bash
    $ setenv F90 gfortran # if SHELL is csh or tcsh

Installing AVNI with all dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
We suggest to install AVNI into its own ``conda`` environment.

The dependency stack is moderate but may take some time (several minutes)
to resolve on some systems via the default ``conda`` solver.

Run in your terminal:

.. code-block:: console

    $ conda create --channel conda-forge --name avni python==3.6
    $ conda activate avni
    $ curl --remote-name https://raw.githubusercontent.com/globalseismology/avni/main/requirements_base.txt
    $ pip install --user -r requirements_base.txt

This will create a new ``conda`` environment called ``avni`` (you can adjust
this by passing a different name via ``--name``) and install all dependencies
into it.


Run the following installation commands to then install AVNI into this environment.

From PyPI
~~~~~~~~~

.. image:: https://img.shields.io/pypi/v/avni?logo=python&logoColor=white
   :target: https://pypi.org/project/avni/

AVNI can be installed from `PyPI <https://pypi.org/project/avni/>`_ using
``pip``::

    $ pip install avni

.. To install all the additional packages that extend AVNI, install using ``pip``
..
    with::

..     pip install avni[all]

.. From Anaconda
.. ~~~~~~~~~~~~~

.. .. image:: https://img.shields.io/conda/vn/conda-forge/avni.svg?logo=conda-forge&logoColor=white
..    :target: https://anaconda.org/conda-forge/avni

..
    To install this package with ``conda`` run::

..     conda install -c conda-forge avni




.. Installing a minimal AVNI with core functionality only
.. ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. If you only need AVNI's core functionality including 2D plotting (but **without
.. 3D visualization**), install via :code:`pip`:

.. .. code-block:: console

..    $ pip install avni

.. or via :code:`conda`:

.. .. code-block:: console

..    $ conda create --strict-channel-priority --channel=conda-forge --name=avni
..    avni-base

.. This will create a new ``conda`` environment called ``avni`` (you can adjust
.. this by passing a different name via ``--name``).

Source / Developers
~~~~~~~~~~~~~~~~~~~

Alternatively, you can install the latest version from GitHub by visiting `AVNI
GitHub`_, and downloading the source (cloning) by running::

    $ git clone https://github.com/globalseismology/avni.git
    $ cd avni
    $ python -m pip install -e .


The latest documentation for the ``main`` branch of AVNI can be found at `AVNI`_.

Installing AVNI for other scenarios
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The :ref:`advanced_setup` page has additional tips and tricks for special
situations (servers, notebooks, CUDA, installing the development version, etc).
The :ref:`contributing` has additional installation instructions for (future)
contributors to AVNI (e.g, extra dependencies for running our tests and building
our documentation).

Python IDEs
===========

Most users find it convenient to write and run their code in an `Integrated
Development Environment`_ (IDE). Some popular choices for scientific Python
development are:

- `Visual Studio Code`_ (often shortened to "VS Code" or "vscode") is a
  development-focused text editor that supports many programming languages in
  addition to Python, includes an integrated terminal console, and has a rich
  ecosystem of packages to extend its capabilities. Installing `Microsoft's
  Python Extension
  <https://marketplace.visualstudio.com/items?itemName=ms-python.python>`__ is
  enough to get most Python users up and running. VS Code is free and
  open-source.

- `PyCharm`_ is an IDE specifically for Python development that provides an
  all-in-one installation (no extension packages needed). PyCharm comes in a
  free "community" edition and a paid "professional" edition, and is
  closed-source.


.. note::

    A more comprehensive testing suite is available after cloning the source
    repository. For details on how to clone and test the AVNI source, please see
    our `contributing guide`_.