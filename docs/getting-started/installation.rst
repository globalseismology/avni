.. _install_ref:

Installation
============

AVNI is supported on Python versions 3.6+. For the best experience, please considering using Anaconda as a virtual
environment and package manager for Python and following the instructions to
install AVNI with Anaconda.


Dependencies
~~~~~~~~~~~~

AVNI is built on top of the Visualization Toolkit (VTK) and NumPy - as such,
the following projects are required dependencies of AVNI:

* `vtk <https://pypi.org/project/vtk/>`_ - AVNI directly inherits types from the VTK library.
* `NumPy <https://pypi.org/project/numpy/>`_ - NumPy arrays provide a core foundation for AVNI's data array access.
* `pillow <https://pypi.org/project/Pillow/>`_ - PIL fork used for saving screenshots.


Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~
AVNI includes several optional dependencies for visualization and reading a variety of additional file formats, including:

* `cmocean <https://pypi.org/project/cmocean/>`_ - Colormaps for Oceanography.
* `colorcet <https://colorcet.holoviz.org/>`_ - Perceptually accurate 256-color colormaps for use with Python.
* `trame <https://github.com/Kitware/trame>`_ - Used for client and server-side rendering in Jupyter.
* `matplotlib <https://pypi.org/project/matplotlib/>`_ - Used for colormaps and 2D plotting with :class:`avni.ChartMPL`.
* `meshio <https://pypi.org/project/meshio/>`_ - Input/Output for many mesh formats.


PyPI
~~~~

.. image:: https://img.shields.io/pypi/v/avni?logo=python&logoColor=white
   :target: https://pypi.org/project/avni/

AVNI can be installed from `PyPI <https://pypi.org/project/avni/>`_
using ``pip``::

    pip install avni

To install all the additional packages that extend AVNI, install using
``pip`` with::

    pip install avni[all]

.. asciinema:: 507562

Anaconda
~~~~~~~~

.. image:: https://img.shields.io/conda/vn/conda-forge/avni.svg?logo=conda-forge&logoColor=white
   :target: https://anaconda.org/conda-forge/avni

To install this package with ``conda`` run::

    conda install -c conda-forge avni

.. asciinema:: 507565



Source / Developers
~~~~~~~~~~~~~~~~~~~

Alternatively, you can install the latest version from GitHub by visiting
`AVNI <https://github.com/avni/avni>`_, and downloading the source
(cloning) by running::

    git clone https://github.com/globalseismology/avni.git
    cd avni
    python -m pip install -e .


The latest documentation for the ``main`` branch of AVNI can be found at
`dev.avni.org <https://dev.avni.org>`_.


Test Installation
~~~~~~~~~~~~~~~~~

You can test your installation by running an example:

.. code:: python

    >>> from avni import demos
    >>> demos.plot_wave()

See other examples and demos:

.. code:: python

    >>> from avni import examples
    >>> from avni import demos

    List all available examples.

    >>> print(dir(examples))

    List all available demos.


.. note::

    A more comprehensive testing suite is available after cloning the source
    repository. For details on how to clone and test the AVNI source, please
    see our `Contributing Guide`_ and specifically, the `Testing`_ section.

.. _Contributing Guide: https://github.com/avni/avni/blob/main/CONTRIBUTING.rst
.. _Testing: https://github.com/avni/avni/blob/main/CONTRIBUTING.rst#user-content-testing


