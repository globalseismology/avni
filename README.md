# AVNI - A toolkit for analyzing and interpreting reference Earth models and data sets

<img src="docs/logos/avnilogo.png" width="256">

## About

AVNI is a Python library for analyzing and interpreting reference Earth models and data sets. The development version has legacy routines in C and Fortran as well, which are interfaced with Python subroutines. Some installation files as well as applets and API access require registration on our [homepage](http://globalseismology.org/register).

Homepage: [avni.globalseismology.org](http://avni.globalseismology.org)

Documentation: readthedocs [html](http://avni.readthedocs.io), [pdf](https://media.readthedocs.org/pdf/avni/latest/avni.pdf)

Source code: [github](https://github.com/geodynamics/avni)

Requests/Bug Reports: [issues](https://github.com/geodynamics/avni/issues)

Contact: **avni@globalseismology.org**

Frequently Asked Questions: [FAQ](docs/FAQ.md)

A core team maintains the public repository and releases versions after benchmarking; if you do not see activity on Github, that does not mean improvements or bug fixes are not underway! We provide APIs that interface with heavy, legacy codes hosted our servers so that AVNI installation remains light to serve various applications. Contact the AVNI team at **avni@globalseismology.org** with any questions or suggestions.

[![Documentation Status](https://readthedocs.org/projects/avni/badge/?version=latest)](https://avni.readthedocs.io/en/latest/?badge=latest) [![codecov](https://codecov.io/gh/globalseismology/avni/branch/devel/graph/badge.svg?token=NTCVjCUfJm)](https://codecov.io/gh/globalseismology/avni) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![PyPI version](https://badge.fury.io/py/avni.svg)](https://badge.fury.io/py/avni) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/avni.svg?style=popout)
![Build Status](https://github.com/globalseismology/avni/actions/workflows/test-and-build-docs-on-pr-self-hosted.yml/badge.svg)

## Requirements

AVNI needs the following python modules installed in the system.
* Python 3.6+
* Python modules provided in [environment.yml](docs/conda/environment.yml)

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Basic Installation

1. Install the [Anaconda Python Distribution](https://www.continuum.io/downloads). We recommend the 64-bit Python 3.7 version.
2. Create a new `avni` environment for AVNI and activate it. You may use the list of
packages that work in tandem with AVNI and have been tested, provided as [environment.yml](docs/conda/environment.yml). Please enterin the following in a terminal:
```
conda env create --name avni -f docs/conda/environment.yml
source $CONDA_ROOT/etc/profile.d/conda.sh
conda activate avni
```
where `CONDA_ROOT` is the location of Anaconda directory (e.g. `/home/user/anaconda3`)

3. (Optional) You may install some mapping toolboxes (Basemap data) with the following options to add topography at all resolutions:
`conda install -c conda-forge basemap-data-hires`
Basemap installation may give errors for the PROJ4 library. We have found it useful to specify the location of the library using the following in our .tcshrc shell configuration
`setenv PROJ_LIB $CONDA_PREFIX/share/proj`
where `CONDA_PREFIX` is the location of Anaconda installation, which should be automatically initialized with step 2.

4. Install the stable version of AVNI and its additional dependencies with
`pip install avni`

5. Create an account on our [Web site](http://globalseismology.org/register) to be allotted an API key for submitting queries to our servers. Enter the key when AVNI is
initialized when run for the first time.

## Start Here

To begin, the user may want to look at these examples to begin to understand what tools are available in AVNI and how workflows may be constructed. Such examples are available as [Jupyter Notebooks](examples/Notebooks) and [Shell scripts](examples/Scripts).

## Contributing to our software

AVNI is a community project that lives by the participation of its
members — i.e., including you! It is our goal to build an inclusive and
participatory community so we are happy that you are interested in
participating! We have collected a set of guidelines and advice on how to get
involved in the community and keep them in the
[CONTRIBUTING.md](CONTRIBUTING.md) file in the repository. Please follow the contributor covenant code of conduct specified in the [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) in all interactions and deliberations.

## Authors

The *Principal Developers* that maintain AVNI are listed in [AUTHORS](AUTHORS). See also the list of [Github contributors](https://github.com/geodynamics/avni/contributors) who participate in this project.

## Versioning

Current software version is provided in the file [version](avni/version.py) file. We use [SemVer](http://semver.org/) for versioning as explained in [versioning_conventions](docs/versioning_conventions.md). For the versions available, see the [tags on this repository](https://github.com/globalseismology/avni/tags).

## License

This software is published under the GNU GPL v3 license - see the [LICENSE](LICENSE) file for details. Please write to **avni@globalseismology.org** or see our [FAQ](docs/FAQ.md) for additional clarifications.

## Acknowledgments

* Funded by the [National Science Foundation](http://nsf.gov) and the [Computational Infrastructure for Geodynamics](http://geodynamics.org).
* We thank Göran Ekström, Adam Dziewonski, and other members of the open source and Earth Science community for their input and non-Github contributions.
* Computational resources maintained by Princeton University accessible through the [AVNI website](http://avni.globalseismology.org). We thank the system administrators for assistance.

<img src="docs/logos/NSF.png" width="100"> &nbsp; &nbsp; &nbsp; &nbsp; <img src="docs/logos/CIG_logo.png" width="200"> &nbsp; &nbsp; &nbsp; &nbsp; <img src="docs/logos/PU-standard.png" width="325">