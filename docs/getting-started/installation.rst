.. include:: ../links.inc

.. _install_ref:

Installation
============

AVNI is supported on Python versions 3.6+. For the best experience, please considering using Anaconda as a virtual
environment and package manager for Python and following the instructions to
install AVNI with Anaconda.


Dependencies
~~~~~~~~~~~~

AVNI is built on top of NumPy - as such,
the following projects are required dependencies of AVNI:

* `NumPy <https://pypi.org/project/numpy/>`_ - NumPy arrays provide a core foundation for AVNI's data array access.

Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~
AVNI includes several optional dependencies for visualization and reading a variety of additional file formats, including:

* `matplotlib <https://pypi.org/project/matplotlib/>`_ - Used for colormaps and 2D plotting.


From PyPI
~~~~~~~~~

.. image:: https://img.shields.io/pypi/v/avni?logo=python&logoColor=white
   :target: https://pypi.org/project/avni/

AVNI can be installed from `PyPI <https://pypi.org/project/avni/>`_
using ``pip``::

    pip install avni

To install all the additional packages that extend AVNI, install using
``pip`` with::

    pip install avni[all]

From Anaconda
~~~~~~~~~~~~~

.. image:: https://img.shields.io/conda/vn/conda-forge/avni.svg?logo=conda-forge&logoColor=white
   :target: https://anaconda.org/conda-forge/avni

To install this package with ``conda`` run::

    conda install -c conda-forge avni

.. asciinema:: 507565



Source / Developers
~~~~~~~~~~~~~~~~~~~

Alternatively, you can install the latest version from GitHub by visiting
`AVNI GitHub`_, and downloading the source
(cloning) by running::

    git clone https://github.com/globalseismology/avni.git
    cd avni
    python -m pip install -e .


The latest documentation for the ``main`` branch of AVNI can be found at `AVNI Docs`_.

.. note::

    A more comprehensive testing suite is available after cloning the source
    repository. For details on how to clone and test the AVNI source, please
    see our `contributing guide`_.

