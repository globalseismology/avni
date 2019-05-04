# REM3D - A toolkit for analyzing and interpreting reference Earth models and data sets

<img src="docs/rem3dlogo.png" width="256">

## About

REM3D is a Python library for analyzing and interpreting reference Earth models and data sets. The development version has legacy routines in C and Fortran as well, which are interfaced with Python subroutines. Some installation files as well as applets and API access require registration on our [homepage](http://rem3d.org/login/register).

Homepage: [rem3d.org](http://rem3d.org)

Documentation: readthedocs [html](http://rem3d.readthedocs.io), [pdf](https://media.readthedocs.org/pdf/rem3d/latest/rem3d.pdf)

Source code: [github](https://github.com/globalseismology/rem3d)

Requests/Bug Reports: [issues](https://github.com/globalseismology/rem3d/issues)

Frequently Asked Questions: [FAQ](docs/FAQ.md)

A core team maintains the public repository and releases versions after benchmarking; if you do not see activity on Github, that does not mean improvements or bug fixes are not underway! We provide APIs that interface with heavy, legacy codes hosted our servers so that REM3D installation remains light to serve various applications. Contact the REM3D team at **info@rem3d.org** with any questions or suggestions.

[![Codacy Badge](https://api.codacy.com/project/badge/Grade/110c5a409f60485f83d442b8834eba2c)](https://www.codacy.com?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=globalseismology/rem3d&amp;utm_campaign=Badge_Grade) [![Build Status](https://travis-ci.com/globalseismology/rem3d.svg?token=Z1JjFn7SrxG1nGGE9y1u&branch=devel)](https://travis-ci.com/globalseismology/rem3d) [![Documentation Status](https://readthedocs.org/projects/rem3d/badge/?version=latest)](https://rem3d.readthedocs.io/en/latest/?badge=latest) [![codecov](https://codecov.io/gh/globalseismology/rem3d/branch/devel/graph/badge.svg?token=NTCVjCUfJm)](https://codecov.io/gh/globalseismology/rem3d) [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0) [![PyPI version](https://badge.fury.io/py/rem3d.svg)](https://badge.fury.io/py/rem3d) [![Gitter](https://badges.gitter.im/rem3d/community.svg)](https://gitter.im/rem3d/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge) ![PyPI - Python Version](https://img.shields.io/pypi/pyversions/rem3d.svg?style=popout)

## Requirements

REM3D needs the following python modules installed in the system.
* Python 3.6+
* Python modules: NumPy, SciPy, Matplotlib, Basemap, Pandas, netCDF4

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Basic Installation

1. Install the [Anaconda Python Distribution](https://www.continuum.io/downloads). We recommend the 64-bit Python 3.7 version.
2. Create a new environment for REM3D and activate it
`conda create --name rem3d --clone base`  
Activate the environment using
`source $CONDA_ROOT/etc/profile.d/conda.csh`  
`conda activate rem3d`  
where `CONDA_ROOT` is the location of Anaconda directory (e.g. `/home/user/anaconda3`)

3. Install some mapping toolboxes (Basemap and its data) with the following options to add topography at all resolutions:  
`conda install netcdf4`  
`conda install -c conda-forge basemap`  
`conda install -c conda-forge basemap-data-hires`  
Basemap installation may give errors for the PROJ4 library. We have found it useful to specify the location of the library using the following in our .tcshrc shell configuration
`setenv PROJ_LIB $CONDA_PREFIX/share/proj`  
where `CONDA_PREFIX` is the location of Anaconda installation, which should be automatically initialized with step 2.
4. Install the stable version of REM3D and its additional dependencies with
`pip install rem3d`  


### Advanced installation (for developers)

Please note the license terms below specific to developers. If you want to use REM3D's development routines, you will need to do the following additional steps. This assumes that you have a free account at <https://github.com>

1. Create an account on our [Web site](http://rem3d.org/login/register) and request access rights from [here](http://rem3d.org/join-us/github).
2. Please accept the invitation through e-mail. Clone the development branch from the REM3D git repository through the terminal. This will create a directory that contains a folder called `rem3d`.  
`git clone -b devel https://username:password@github.com/globalseismology/rem3d.git`  
3. Install REM3D module by opening a terminal window, navigating to `rem3d` directory and entering  
`pip install -e . --user`  
This lets you to keep working on files inside the Github folder without recompiling the codes.

We maintain 3 major branches for our client libraries and these are relevant to public releases. Read/write access to these are restricted to main administrators and developers:
* `devel` — Active development occurs on this branch or in branches that spin off from this.
* `public` — Development for bug fixes happens here. We also bump versions and update the changelog on this branch.
* `master` — We squash commits from the release branch into single release commits on this branch as well as tagging releases.

New branches may be created for individual projects and the relevant team. Please clone the `devel` branch to build upon the latest codes  
`git checkout -b new_branch devel`  
You can push this locally created branch to the remote `globalseismology/rem3d` with  
`git push -u origin new_branch`  

## Start Here

To begin, the user may want to look at these examples to begin to understand what tools are available in REM3D and how workflows may be constructed. Such examples are available as [Jupyter Notebooks](examples/Notebooks) and [Shell scripts](examples/Scripts).

## Contributing

Please read [contributing.md](docs/contributing.md) for details on our code of conduct, and the process for submitting pull requests to our public repository. Contributions to our private repository are made to custom branches and will not be merged to the public repository without approval of the whole REM3D team.

### Authors

* **Pritwiraj Moulik** - *Primary administrator* - [github](https://github.com/pmoulik)
* **Ved Lekic** - *Co-administrator* - [github](https://github.com/vedlekic)

See also the list of [Github contributors](https://github.com/globalseismology/rem3d/contributors) who participated in this project.

## License

The master branch of this project is licensed under the GNU GPL v3 or newer, is synced with public repositories and may be installed in any platform - see the [LICENSE](LICENSE) file for details. All other branches are copyrighted, must remain in our private organizational repository and may not be installed in platforms without permission from code contributors and repository administrators. Please write to **info@rem3d.org** or see our [FAQ](docs/FAQ.md) for more details.

### Versioning

We use [SemVer](http://semver.org/) for versioning as explained in [versioning_conventions](docs/versioning_conventions.md). For the versions available, see the [tags on this repository](https://github.com/globalseismology/rem3d/tags).

## Acknowledgments

* Funded by the National Science Foundation and the David and Lucile Packard Foundation.
* We thank Göran Ekström, Adam Dziewonski, other members of the open source and Earth Science community for their input and non-Github contributions.
* Computational resources maintained at the University of Maryland at College Park. We thank Todd Karowski and other system administrators for assistance.

<img src="docs/NSF.png" width="100"> <img src="docs/packard.png" width="200"> <img src="docs/UMD.png" width="100">