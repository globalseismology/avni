.. include:: ../links.inc

Overview of the AVNI tools suite
================================

AVNI is an open-source Python module for processing, analysis, and
visualization of functional neuroimaging data (EEG, MEG, sEEG, ECoG, and
fNIRS). There are several related or interoperable software packages that you
may also want to install, depending on your analysis needs.

Related software
^^^^^^^^^^^^^^^^

- AVNI was the initial stage of this project,
  providing a set of interrelated command-line and GUI programs focused on
  computing cortically constrained Minimum Norm Estimates from MEG and EEG
  data. These tools were written in C by Matti Hämäläinen, and are
  documented `here <AVNI manual_>`_. See :ref:`install_avni_c` for installation
  instructions.

- AVNI reimplements the functionality of AVNI, extends considerably the
  analysis and visualization capabilities, and adds support for additional data
  types like functional near-infrared spectroscopy (fNIRS). AVNI is
  collaboratively developed and has more than 200 contributors.

- :ref:`AVNI MATLAB <avni_matlab>` provides a MATLAB interface to the .fif file
  format and other AVNI data structures, and provides example MATLAB
  implementations of some of the core analysis functionality of AVNI. It is
  distributed alongside AVNI, and can also be downloaded from the `AVNI-MATLAB
  git repository`_.

- :ref:`AVNIPP <avni_cpp>` provides core AVNI functionality implemented in
  C++ and is primarily intended for embedded and real-time applications.

There is also a growing ecosystem of other Python packages that work alongside
AVNI, including packages for:

.. note:: Something missing?
    :class: sidebar

    If you know of a package that is related but not listed here, feel free to
    :ref:`make a pull request <contributing>` to add it to this list.

- a graphical user interface for AVNI (`AVNILAB`_)
- easily importing MEG data from the Human Connectome Project for
  use with AVNI (`AVNI-HCP`_)
- managing AVNI projects so that they comply with the `Brain
  Imaging Data Structure`_ specification (`AVNI-BIDS`_)
- automatic bad channel detection and interpolation (`autoreject`_)
- convolutional sparse dictionary learning and waveform shape estimation
  (`alphaCSC`_)
- independent component analysis (ICA) with good performance on real data
  (`PICARD`_)
- phase-amplitude coupling (`pactools`_)
- representational similarity analysis (`rsa`_)
- microstate analysis (`microstate`_)
- connectivity analysis using dynamic imaging of coherent sources (DICS)
  (`conpy`_)
- general-purpose statistical analysis of M/EEG data (`eelbrain`_)
- post-hoc modification of linear models (`posthoc`_)
- a python implementation of the Preprocessing Pipeline (PREP) for EEG data
  (`pyprep`_)
- automatic multi-dipole localization and uncertainty quantification with
  the Bayesian algorithm SESAME (`sesameeg`_)
- GLM and group level analysis of near-infrared spectroscopy data (`avni-nirs`_)
- High-level EEG Python library for all kinds of EEG inverse solutions (`invertmeeg`_)
- All-Resolutions Inference (ARI) for statistically valid circular inference
  and effect localization (`AVNI-ARI`_)


What should I install?
^^^^^^^^^^^^^^^^^^^^^^

If you intend only to perform ERP, ERF, or other sensor-level analyses,
:ref:`AVNI <standard-instructions>` is all you need. If you prefer to
work with
shell scripts and the Unix command line, or prefer MATLAB over Python, probably
all you need is :doc:`AVNI <avni_c>` — the AVNI MATLAB toolbox is distributed
with it — although note that the C tools and the MATLAB toolbox are less
actively developed than the AVNI module, and hence are considerably less
feature-complete.



Getting help
^^^^^^^^^^^^

Help with installation is available through the `AVNI Forum`_. See the
:ref:`help` page for more information.

