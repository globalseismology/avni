.. include:: ../links.inc

.. _advanced_setup:

Advanced setup
==============

Working with Jupyter Notebooks and JupyterLab
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you like using Jupyter notebooks, you should also update the "base"
conda environment to include the ``nb_conda_kernels`` package; this will
make it easier to use AVNI in Jupyter Notebooks launched from the
Anaconda GUI:

.. code-block:: console

    $ conda install --name=base nb_conda_kernels

When using AVNI within IPython or a Jupyter notebook, we strongly
recommend using the Qt matplotlib backend for fast and correct rendering. On
Linux, for example, Qt is the only matplotlib backend for which 3D rendering
will work correctly. On macOS, certain matplotlib functions might not work as
expected on backends other than Qt. Enabling Qt can be accomplished when
starting IPython from a terminal:

.. code-block:: console

    $ ipython --matplotlib=qt

or in a Jupyter Notebook, you can use the "magic" command:

.. code-block:: ipython

    In [1]: %matplotlib qt

This will create separate pop-up windows for each figure, and has the advantage
that the 3D plots will retain rich interactivity (so, for example, you can
click-and-drag to rotate cortical surface activation maps).

If you are creating a static notebook or simply prefer Jupyter's inline plot
display, AVNI will work with the standard "inline" magic:

.. code-block:: ipython

    In [1]: %matplotlib inline

but some functionality will be lost. For example, PyVista scenes will still
pop-up a separate window, but only one window at a time is possible, and
interactivity within the scene is limited in non-blocking plot calls.

.. admonition:: |windows| Windows
  :class: note

  If you are using AVNI on Windows through IPython or Jupyter, you might
  also have to use the IPython magic command ``%gui qt`` (see `here
  <https://github.com/ipython/ipython/issues/10384>`_). For example:

  .. code-block:: ipython

     In [2]: %gui qt

If you installed the ``nb_conda_kernels`` package into your ``base``
environment (as recommended), you should be able to launch ``avni``-capable
notebooks from within the Anaconda Navigator GUI without having to explicitly
switch to the ``avni`` environment first; look for ``Python [conda env:avni]``
when choosing which notebook kernel to use. Otherwise, be sure to activate the
``avni`` environment before launching the notebook.

If you use another Python setup and you encounter some difficulties please
report them on the `AVNI Forum`_ or on the `GitHub issues page`_ to get
assistance.

It is also possible to interact with the 3D plots without installing Qt by using
the notebook 3d backend:

.. code-block:: ipython

   In [1]: import avni
   In [2]: avni.viz.set_3d_backend("notebook")


The notebook 3d backend requires PyVista to be installed along with other packages,
please follow :ref:`standard-instructions`.

Using the development version
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See :ref:`installing_main` for how to do a one-time update to the latest
development version of AVNI. If you plan to contribute to AVNI, or
just prefer to use git rather than pip to make frequent updates, there are
instructions for installing from a ``git clone`` in the :ref:`contributing`.


.. _other-py-distros:

Other Python distributions
^^^^^^^^^^^^^^^^^^^^^^^^^^

While the `Anaconda`_ Python distribution provides many conveniences, other
distributions of Python should also work with AVNI.  In particular,
`Miniconda`_ is a lightweight alternative to Anaconda that is fully compatible;
like Anaconda, Miniconda includes the ``conda`` command line tool for
installing new packages and managing environments; unlike Anaconda, Miniconda
starts off with a minimal set of around 30 packages instead of Anaconda's
hundreds. See the `installation instructions for Miniconda`_ for more info.
A similar alternative is `MiniForge`_, which uses the ``conda-forge`` channel
as the default source for package installation (saving you the trouble of
typing ``--channel=conda-forge`` with each ``conda install`` command).

.. warning::

    If you have the ``PYTHONPATH`` or ``PYTHONHOME`` environment variables set,
    you may run into difficulty using Anaconda. See the
    `Anaconda troubleshooting guide`_ for more information. Note that it is
    easy to switch between ``conda``-managed Python installations and the
    system Python installation using the ``conda activate`` and ``conda
    deactivate`` commands, so you may find that after adopting Anaconda it is
    possible (indeed, preferable) to leave ``PYTHONPATH`` and ``PYTHONHOME``
    permanently unset.


It is also possible to use a system-level installation of Python (version
|min_python_version| or higher) and use ``pip`` to install AVNI and its
dependencies, using the provided requirements file:

.. code-block:: console

    $ curl --remote-name https://raw.githubusercontent.com/globalseismology/avni/main/requirements.txt
    $ pip install --user -r requirements.txt

Other configurations will probably also work, but we may be unable to offer
support if you encounter difficulties related to your particular Python
installation choices.


