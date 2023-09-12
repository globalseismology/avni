.. title:: AVNI
.. include:: links.inc
.. include:: buttons.inc
.. raw:: html

    <div class="banner">
        <a href="./index.html"><center><img src="_static/logos/logo_avni_color_withname.png" alt="avni" width="50%"/></a>
        <h2>Analysis and Visualization toolkit for plaNetary Inferences</h2>
    </div>

.. toctree::
   :hidden:

   Getting Started <getting-started/index>
   Documentation <overview/index>
   API Reference <python_reference>
   Get help <overview/get_help>
   Development <overview/development>

**Want to test-drive AVNI?** Check out our live examples on Hubzero: |Hubzero|

AVNI is a software ecosystem for analyzing and interpreting planetary models and data sets that were initially
designed for the three-dimensional reference Earth model project |REM3D|.
The codes are primarily written in Python with interfaces to legacy routines in C and Fortran.
Some installation files as well as applets and API access require registration on the |register| link.

* Open-source Python package with APIs to handle intensive queries |GNU GPL|
* Introduce HDF5 storage formats for planetary models and processed seismic data
* Interactive web-based visualization tools for data and model exploration
* Formulate and benchmark solvers for rapid data validation of models

*Share this project on Twitter*: |tweet|

.. grid::

   .. grid-item-card:: AVNI is a CIG affiliated project
      :link: https://geodynamics.org
      :class-title: avni-card-title

      .. image:: _static/logos/CIG_logo_with_text.png
         :target: https://geodynamics.org
         :alt: CIG affiliated projects
         :height: 60px

.. jupyter-execute::

   import avni
   print(avni.__version__)

Links
^^^^^

+----------------------+----------------+-------------+
| Website              | |home|         |  |Github|   |
+----------------------+----------------+-------------+
| Deployment           | |pypi_button|  |  |Hubzero|  |
+----------------------+----------------+-------------+
| License              | |GNU GPL|      |             |
+----------------------+----------------+-------------+
| Community            | |discuss|      |  |issues|   |
+----------------------+----------------+-------------+

Applications
^^^^^^^^^^^^

.. frontpage gallery is added by a conditional in _templates/layout.html
