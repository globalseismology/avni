.. title:: AVNI

.. raw:: html

    <div class="banner">
        <a href="./examples/index.html"><center><img src="_static/logo_avni_color.png" alt="avni" width="30%"/></a>
        <h2>Analysis and Visualization toolkit for plaNetary Inferences</h2>
    </div>

|


Overview
********
AVNI is...

* *Pythonic VTK*: a high-level API to the `Visualization Toolkit`_ (VTK)
* mesh data structures and filtering methods for spatial datasets
* 3D plotting made simple and built for large/complex data geometries

.. _Visualization Toolkit: https://vtk.org

AVNI is a helper library for the Visualization Toolkit (VTK) that
takes a different approach on interfacing with VTK through NumPy and
direct array access.  This package provides a Pythonic,
well-documented interface exposing VTK's powerful visualization
backend to facilitate rapid prototyping, analysis, and visual
integration of spatially referenced datasets.

This module can be used for scientific plotting for presentations and
research papers as well as a supporting module for other mesh
dependent Python modules.

.. |tweet| image:: https://img.shields.io/twitter/url.svg?style=social&url=http%3A%2F%2Fshields.io
   :target: https://twitter.com/intent/tweet?text=Check%20out%20this%20project%20for%20Earth%20and%20other%20planets%20in%20Python&url=https://avni.globalseismology.org&hashtags=3D,visualization,Python,vtk,mesh,plotting,AVNI

Share this project on Twitter: |tweet|


.. |binder| image:: https://static.mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/pyvista/pyvista-examples/master
   :alt: Launch on Binder

Want to test-drive AVNI? Check out our live examples on MyBinder: |binder|

.. grid::

   .. grid-item-card:: AVNI is a NumFOCUS affiliated project
      :link: https://numfocus.org/sponsored-projects/affiliated-projects
      :class-title: pyvista-card-title

      .. image:: https://raw.githubusercontent.com/numfocus/templates/master/images/numfocus-logo.png
         :target: https://numfocus.org/sponsored-projects/affiliated-projects
         :alt: NumFOCUS affiliated projects
         :height: 60px


.. toctree::
   :hidden:

   getting-started/index
   user-guide/index
   examples/index
   api/index
   extras/index


Brief Examples
**************
Here are some brief interactive examples that demonstrate how you
might want to use AVNI:


.. jupyter-execute::
   :hide-code:

   # Configure for panel
   import avni
   print(avni.__version__)


Maps and Geoscience
~~~~~~~~~~~~~~~~~~~
Download the surface elevation map of Mount St. Helens and plot it.

.. jupyter-execute::

   import avni
   print(avni.__version__)


Finite Element Analysis
~~~~~~~~~~~~~~~~~~~~~~~
Plot the 'X' component of elastic stress of a 3D notch specimen.

.. jupyter-execute::

   import avni
   print(avni.__version__)

Simple Point Cloud with Numpy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Easily integrate with NumPy and create a variety of geometries and plot
them.  You could use any geometry to create your glyphs, or even plot
the points directly.

.. jupyter-execute::

   import avni
   print(avni.__version__)

Plot a Spline
~~~~~~~~~~~~~
Generate a spline from an array of NumPy points.

.. jupyter-execute::

   import avni
   print(avni.__version__)

Boolean Operations on Meshes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Subtract a sphere from a cube mesh.

.. jupyter-execute::

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

Translating
***********
The recommended way for new contributors to translate AVNI's
documentation is to join the translation team on Transifex.

There is a `pyvista translation page`_ for pyvista (main) documentation.

1. Login to transifex_ service.
2. Go to `pyvista translation page`_.
3. Click ``Request language`` and fill form.
4. Wait acceptance by transifex pyvista translation maintainers.
5. (After acceptance) Translate on transifex.
6. We can host the translated document in `GitHub Pages`_ by creating `GitHub repository`_.
7. Translation is backed up in `pyvista-doc-translations`_.

Details can be found here: https://docs.transifex.com/getting-started-1/translators

.. _`pyvista translation page`: https://www.transifex.com/tkoyama010/pyvista-doc/
.. _Transifex: https://www.transifex.com/
.. _`GitHub Pages`: https://pyvista.github.io/pyvista-docs-dev-ja/index.html
.. _`GitHub repository`: https://github.com/pyvista/pyvista-docs-dev-ja
.. _`pyvista-doc-translations`: https://github.com/pyvista/pyvista-doc-translations


Status
******

.. |pypi| image:: https://img.shields.io/pypi/v/avni?logo=python&logoColor=white
   :target: https://pypi.org/project/avni/

.. |MIT| image:: https://img.shields.io/badge/License-MIT-yellow.svg
   :target: https://opensource.org/licenses/MIT

.. |PyPIact| image:: https://img.shields.io/pypi/dm/pyvista.svg?label=PyPI%20downloads
   :target: https://pypi.org/project/pyvista/

.. |discuss| image:: https://img.shields.io/badge/GitHub-Issues-green?logo=github
   :target: https://github.com/globalseismology/avni/issues

.. |python| image:: https://img.shields.io/badge/python-3.6+-blue.svg
   :target: https://www.python.org/downloads/

+----------------------+----------------+-------------+
| Deployment           | |pypi|         |  |PyPIact|  |
+----------------------+----------------+-------------+
| License              | |MIT|          |             |
+----------------------+----------------+-------------+
| Community            | |discuss|      |             |
+----------------------+----------------+-------------+


Indices and tables
*************

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
