.. include:: ../links.inc

.. _documentation_overview:

Documentation overview
======================

.. note::

   If you haven't already installed AVNI, please take a look
   at our :ref:`installation guides<installers>`. Please also kindly find some
   resources for :doc:`learn_python` if you need to.


The documentation for AVNI is divided into three main sections:

1. The :doc:`tutorials/index` page contains narrative tutorials in the Jupyter Notebook format.
   These executable notebooks contain explanations, sample
   code, and expected output for the most common AVNI analysis tasks. The
   emphasis is on thorough explanations that get new users up to speed quickly,
   at the expense of covering only a limited number of topics.

2. The :doc:`scripts/index` page provides scripts that can be executed from the
   command line. These are either in Python or BASH/TCSH shell formats. The
   emphasis is on some useful applications that do not need the interactivity of the tutorials.

3. The :doc:`API reference <../python_reference>` provides documentation for
   the classes, functions and methods in the AVNI codebase. This is the
   same information that is rendered when running
   :samp:`help(avni.{<function_name>})` in an interactive Python session, or
   when typing :samp:`avni.{<function_name>}?` in an IPython session or Jupyter
   notebook.

4. The :doc:`frequenty asked questions (FAQ) <faq>` provides answers to several
   common questions about the AVNI codebase. Feel feel free to reach out to us
   to get :doc:`more help <get_help>`

.. important::
   Checkout the :doc:`project governance <governance>`,
   :doc:`design philosophy <design_philosophy>` and :doc:`code of conduct <CODE_OF_CONDUCT>`
   pages for a broad overview of the project expectations and outlook.

.. toctree::
   :hidden:

   governance
   design_philosophy
   tutorials/index
   scripts/index
   CODE_OF_CONDUCT
