.. include:: ../links.inc

Getting Started
***************

This guide is here to help you start creating interactive 3D plots with AVNI
with the help of our examples and tutorials.


..
    tab-set::

..    .. tab-item:: JupyterLab

..       Here's a quick demo of AVNI running within `JupyterLab
..       <https://jupyter.org/>`_.

..       .. raw:: html

..          <video width="100%" height="auto" controls autoplay muted> <source
..            src="https://tutorial.avni.org/_static/avni_jupyterlab_demo.mp4"
..            type="video/mp4" style="margin-left: -220px; margin-right: -10.5%">
..            Your browser does not support the video tag.  </video>

..    .. tab-item:: IPython

..       Here's a quick demo of AVNI running within a terminal using `IPython
..       <https://ipython.org/>`_.

..       .. raw:: html

..          <video width="100%" height="auto" controls autoplay muted> <source
..            src="https://tutorial.avni.org/_static/avni_ipython_demo.mp4"
..            type="video/mp4"> Your browser does not support the video tag.
..            </video>


.. toctree::
   :hidden:

   cite
   installation
   advanced

  .. why
  .. installers
  .. updating
  .. check_installation

Installation
============
The first prerequisite for installing AVNI is Python itself. If you don’t
have Python yet and want the simplest way to get started, we recommend you use
the `Anaconda Distribution <https://www.anaconda.com/>`_. The second
prerequisite is the `gfortran Fortran compiler <https://gcc.gnu.org/wiki/GFortran>`_.
Once Fortran is installed, please set the environment variable `F90` in your
command line shell to `gfortran`.

.. .. grid:: 2

..     .. grid-item-card:: Working with conda?
..        :class-title: avni-card-title

..        AVNI is available on `conda-forge
..        <https://anaconda.org/conda-forge/avni>`_.

..        .. code-block:: bash

..           conda install -c conda-forge avni


..     .. grid-item-card:: Prefer pip?
..        :columns: auto
..        :class-title: avni-card-title

..        AVNI can be installed via pip from `PyPI`_.

..        .. code-block:: bash

..           pip install avni


.. when https://github.com/executablebooks/sphinx-design/issues/66 is fixed,
   prepend |cloud-arrow-down|   to the "Download installers" button text
   and |wrench|   to the "Setup instructions" button text

.. grid:: 2


    .. grid-item-card::
        :text-align: center

        .. rst-class:: font-weight-bold mb-0

            Install via ``pip``

        .. rst-class:: card-subtitle text-muted mt-0

            For Advanced Users

        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        .. image:: ../_static/avni_installer_console.png
           :alt: Terminal Window

        **Already familiar with Python?**
        Follow our advanced setup instructions for ``pip``!
        +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        .. button-ref:: install_ref
            :ref-type: ref
            :color: primary
            :shadow:
            :class: font-weight-bold

            Setup Instructions

..
    .. grid::

..    .. grid-item-card:: In-depth instructions?
..       :link: install_ref
..       :link-type: ref
..       :class-title: avni-card-title

..       Installing a specific version? Installing from source? Check the
..       :ref:`install_ref` page.



Support
=======

For general questions about the project, its applications, or about software
usage, please create a discussion in `AVNI Forum`_
where the community can collectively address your questions. You are also send one of the developers an email.
The project support team can be reached at `avni@globalseismology.org`_.


Citing AVNI
===========


If you are using AVNI in your scientific research, please help our scientific
visibility by citing our work! Head over to :ref:`cite` page to learn more
about citing AVNI.
