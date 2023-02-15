Getting Started
***************

This guide is here to help you start creating interactive 3D plots with AVNI
with the help of our examples and tutorials.


.. tab-set::

   .. tab-item:: JupyterLab

      Here's a quick demo of AVNI running within `JupyterLab
      <https://jupyter.org/>`_.

      .. raw:: html

         <video width="100%" height="auto" controls autoplay muted> <source
           src="https://tutorial.avni.org/_static/avni_jupyterlab_demo.mp4"
           type="video/mp4" style="margin-left: -220px; margin-right: -10.5%">
           Your browser does not support the video tag.  </video>

   .. tab-item:: IPython

      Here's a quick demo of AVNI running within a terminal using `IPython
      <https://ipython.org/>`_.

      .. raw:: html

         <video width="100%" height="auto" controls autoplay muted> <source
           src="https://tutorial.avni.org/_static/avni_ipython_demo.mp4"
           type="video/mp4"> Your browser does not support the video tag.
           </video>


.. toctree::
   :hidden:

   why
   authors
   installation
   connections
   external_examples


Installation
============
The only prerequisite for installing AVNI is Python itself. If you don’t
have Python yet and want the simplest way to get started, we recommend you use
the `Anaconda Distribution <https://www.anaconda.com/>`_.

.. grid:: 2

    .. grid-item-card:: Working with conda?
       :class-title: avni-card-title

       AVNI is available on `conda-forge
       <https://anaconda.org/conda-forge/avni>`_.

       .. code-block:: bash

          conda install -c conda-forge avni


    .. grid-item-card:: Prefer pip?
       :columns: auto
       :class-title: avni-card-title

       AVNI can be installed via pip from `PyPI
       <https://pypi.org/project/avni>`__.

       .. code-block:: bash

          pip install avni


.. grid::

   .. grid-item-card:: In-depth instructions?
      :link: install_ref
      :link-type: ref
      :class-title: avni-card-title

      Installing a specific version? Installing from source? Check the
      :ref:`install_ref` page.


First Steps
===========
We've provided a variety of resources for you to get used to AVNI's API
through a range of examples and tutorials.


.. grid::

   .. grid-item-card:: Tutorial
      :link: https://tutorial.avni.org/tutorial.html
      :class-title: avni-card-title

      Probably the best way for you to get used to AVNI is to visit our
      dedicated `tutorial <https://tutorial.avni.org/tutorial.html>`_.

..
   This code is used in the plot in the card.

.. avni-plot::
   :include-source: False
   :context:

   >>> bunny_cpos = [( 0.14826, 0.275729,  0.4215911),
   ...               (-0.01684, 0.110154, -0.0015369),
   ...               (-0.15446, 0.939031, -0.3071841)]


.. grid:: 2

   .. grid-item-card:: Why AVNI?
      :link: why_avni
      :link-type: ref
      :class-title: avni-card-title

      Learn more about why we created AVNI as an interface to the
      Visualization Toolkit (VTK).

      .. code:: python

         import avni
         mesh = avni.read('bunny.stl')
         mesh.plot()

      .. avni-plot::
         :include-source: False
         :context:

         from avni import examples
         mesh = examples.download_bunny()
         mesh.plot(cpos=bunny_cpos, anti_aliasing='ssao')


   .. grid-item-card:: Authors & Citation
      :link: authors_ref
      :link-type: ref
      :class-title: avni-card-title

      Using AVNI in your research? Please consider citing or acknowledging
      us.  We have a `JOSS Publication`_!

      .. image:: ../images/user-generated/joss.png
         :target: https://joss.theoj.org/papers/10.21105/joss.01450

.. grid::

   .. grid-item-card:: See AVNI in External Efforts
      :link: external_examples
      :link-type: ref
      :class-title: avni-card-title

      Take a look at third party projects using AVNI.


Support
=======

For general questions about the project, its applications, or about software
usage, please create a discussion in `avni/discussions`_
where the community can collectively address your questions. You are also
welcome to join us on Slack_ or send one of the developers an email.
The project support team can be reached at `info@avni.org`_.

.. _avni/discussions: https://github.com/avni/avni/discussions
.. _Slack: http://slack.avni.org
.. _info@avni.org: mailto:info@avni.org


Citing AVNI
==============

There is a `paper about AVNI <https://doi.org/10.21105/joss.01450>`_!

If you are using AVNI in your scientific research, please help our scientific
visibility by citing our work! Head over to :ref:`citation_ref` to learn more
about citing AVNI.

.. _JOSS Publication: https://joss.theoj.org/papers/10.21105/joss.01450
