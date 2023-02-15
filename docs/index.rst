.. title:: AVNI

.. raw:: html

    <div class="banner">
        <a href="./examples/index.html"><center><img src="_static/logo_avni_color.png" alt="avni" width="30%"/></a>
        <h2>Analysis and Visualization toolkit for plaNetary Inferences</h2>
    </div>

|


Overview
********

AVNI is a Python library for analyzing and interpreting reference Earth models and data sets. The development version has legacy routines in C and Fortran as well, which are interfaced with Python subroutines. Some installation files as well as applets and API access require registration on `our homepage`_.

* Open-source Python package with APIs to handle intensive queries
* Introduce HDF5 storage formats for planetary models and processed seismic data
* Interactive web-based visualization tools for data and model exploration
* Formulate and benchmark solvers for rapid data validation of models


.. _our homepage: http://globalseismology.org


.. |tweet| image:: https://img.shields.io/twitter/url.svg?style=social&url=http%3A%2F%2Fshields.io
   :target: https://twitter.com/intent/tweet?text=Check%20out%20this%20project%20for%20Earth%20and%20other%20planets%20in%20Python&url=https://avni.globalseismology.org&hashtags=3D,visualization,Python,vtk,mesh,plotting,AVNI

Share this project on Twitter: |tweet|


.. |Hubzero| image:: https://img.shields.io/badge/Launch-Hubzero-orange.svg
   :target: https://geodynamics.org/tools/avninotebooks
   :alt: Launch on Hubzero

Want to test-drive AVNI? Check out our live examples on Hubzero: |Hubzero|

.. grid::

   .. grid-item-card:: AVNI is a CIG affiliated project
      :link: https://geodynamics.org
      :class-title: avni-card-title

      .. image:: _static/logos/CIG_logo_with_text.png
         :target: https://numfocus.org/sponsored-projects/affiliated-projects
         :alt: CIG affiliated projects
         :height: 60px



.. toctree::
   :hidden:

   Getting Started <getting-started/index>
   Install <install/index>
   Documentation <overview/index>
   API Reference <python_reference>
   Get help <overview/get_help>
   Development <overview/development>


Brief Examples
**************
Here are some brief interactive examples that demonstrate how you
might want to use AVNI:


.. jupyter-execute::
   :hide-code:

   # Configure for panel
   import avni
   print(avni.__version__)

   
Plot Volumetric Data
~~~~~~~~~~~~~~~~~~~~
Plot the :math:`3d_{xy}` orbital of a hydrogen atom.

.. note::
   This example requires `sympy <https://www.sympy.org/>`_.

.. jupyter-execute::

   import avni
   print(avni.__version__)

Status
******

.. |pypi| image:: https://img.shields.io/pypi/v/avni?logo=python&logoColor=white
   :target: https://pypi.org/project/avni/

.. |MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT

.. |issues| image:: https://img.shields.io/badge/GitHub-Issues-green?logo=github
   :target: https://github.com/globalseismology/avni/issues

.. |python| image:: https://img.shields.io/badge/python-3.6+-blue.svg
   :target: https://www.python.org/downloads/
   
.. |discuss| image:: https://img.shields.io/badge/GitHub-Discussions-green?logo=github
   :target: https://github.com/globalseismology/avni/discussions

+----------------------+----------------+-------------+
| Deployment           | |pypi|         |             |
+----------------------+----------------+-------------+
| License              | |MIT|          |             |
+----------------------+----------------+-------------+
| Community            | |discuss|      |  |issues|   |
+----------------------+----------------+-------------+


Indices and tables
*************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
