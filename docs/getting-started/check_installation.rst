.. include:: ../links.inc

.. _testing-installation:


Testing your installation
=========================

To make sure AVNI was installed correctly, type the following command in a
terminal::

    python -c "import avni; avni.sys_info()"

.. hint::
   If you installed AVNI using one of our installers, enter the above
   command in the **Prompt**.

This should display some system information along with the versions of AVNI and
its dependencies. Typical output looks like this::

    Platform:      Linux-5.0.0-1031-gcp-x86_64-with-glibc2.2.5
    Python:        3.8.1 (default, Dec 20 2019, 10:06:11)  [GCC 7.4.0]
    Executable:    /home/travis/virtualenv/python3.8.1/bin/python
    CPU:           x86_64: 2 cores
    Memory:        7.8 GB

    avni:           0.1.dev0
    numpy:         1.19.0.dev0+8dfaa4a {blas=openblas, lapack=openblas}
    scipy:         1.5.0.dev0+f614064
    matplotlib:    3.2.1 {backend=QtAgg}

    sklearn:       0.22.2.post1
    numba:         0.49.0
    nibabel:       3.1.0
    cupy:          Not found
    pandas:        1.0.3
    dipy:          1.1.1
    pyvista:       0.25.2 {pyvistaqt=0.1.0}
    vtk:           9.0.0
    qtpy:          2.0.1 {PySide6=6.2.4}


.. dropdown:: If you get an error...
    :color: danger
    :icon: alert-fill

    .. rubric:: If you see an error like:

    ::

        Traceback (most recent call last):
          File "<string>", line 1, in <module>
        ModuleNotFoundError: No module named 'avni'

    This suggests that your environment containing AVNI is not active. If you
    followed the setup for 3D plotting/source analysis (i.e., you installed to a
    new ``avni`` environment instead of the ``base`` environment) try running
    ``conda activate avni`` first, and try again. If this works, you might want
    to set your terminal to automatically activate the ``avni`` environment each
    time you open a terminal::

        echo conda activate avni >> ~/.bashrc    # for bash shells
        echo conda activate avni >> ~/.zprofile  # for zsh shells

If something else went wrong during installation and you can't figure it out,
check out the :ref:`advanced_setup` instructions to see if your problem is
discussed there. If not, the `AVNI Forum`_ is a good resources for
troubleshooting installation problems.

.. highlight:: python

