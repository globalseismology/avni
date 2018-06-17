# REM3D - A toolkit for analyzing and interpreting reference Earth models and data sets.

<img src="docs/rem3dlogo.png" width="256">

## About

REM3D is a Python library for analyzing and interpreting reference Earth models and data sets. The development version has modules in MATLAB, C and Fortran as well.

Homepage: [rem3d.org](https://maurya.umd.edu)

Documentation: readthedocs [html](http://rem3d.readthedocs.io), [pdf](https://media.readthedocs.org/pdf/rem3d/latest/rem3d.pdf)  

Source code: [github](https://github.com/globalseismology/rem3d) 

Contact the REM3D team at info@rem3d.org with any questions or suggestions.

## Requirements

REM3D needs the following python modules installed in the system.
* Python 2.7.x or Python 3.4+
* Python modules:
  NumPy, SciPy, Matplotlib, Basemap, Pandas, netCDF4
  
Users may also want to install MATLAB for some codes. We have tested our versions on MATLAB 2015.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Basic Installation

1. Install the [Anaconda Python Distribution](https://www.continuum.io/downloads). We recommend the 64-bit Python 2.7 version. 
2. Install some mapping toolboxes (pandas, Basemap and its data) with the following options to add topography at all resolutions:
`conda install -c anaconda pandas`
`conda install -c anaconda netcdf4`
`conda install -c conda-forge basemap`
`conda install -c conda-forge basemap-data-hires`
3. Install the stable version of REM3D and its additional dependencies with `pip install rem3d`. 

### Advanced installation (for developers)

If you want to use REM3D's development routines, you will need to do the following additional steps. This assumes that you have a free account at <https://github.com>

1. Create an account on our [Web site](https://maurya.umd.edu/login/register) and request access rights from [here](https://maurya.umd.edu/join-us/github).
2. Please accept the invitation through e-mail. Clone the development branch REM3D git repository and the appro from the terminal `git clone -b devel https://username:password@github.com/globalseismology/rem3d.git`. This will create a directory that contains a folder call `rem3d`.
3. Install REM3D module by opening a terminal window, navigating to `rem3d` directory and entering
`pip install -e . --user`. This lets you to keep working on files inside the Github folder without recompiling the codes..

## Start Here

To begin, the user may want to look at these examples to begin to understand
what tools are available in REM3D and how values are calculated. Such examples are available in the [examples](examples) folder.

## About scripting in Python

REM3D has the advantage of being adaptable and extendable in easy scripts. The downside might be that we do not
provide a graphical user interface. For those of you who are not familiar  with python, we suspect it will still be 
relatively easy to adapt the scripts for computations and plotting. 
Here are some specific features and pitfalls on Python:

* Python uses specific indentation. A script might fail if a code block is not indented correctly. We use four spaces and no tabs, mixing these can give trouble.
* Indices require square brackets and function or method calls parentheses (mainly different from Matlab).
* The first index of an array is 0 (e.g. x[0])
* Put dots after numbers to make them floats instead of integers (e.g. 5/3 will give 1 (Python 2.x rounds downward), while 5./3. will give 1.66666666667)

## Contributing

Please read [contributing.md](docs/contributing.md) for details on our code of conduct, and the process for submitting pull requests to us.

### Versioning

We use [SemVer](http://semver.org/) for versioning as explained in [versioning_conventions.md](docs/versioning_conventions.md). For the versions available, see the [tags on this repository](https://github.com/globalseismology/rem3d/tags). 

### Authors

* **Pritwiraj Moulik** - *Primary administrator* - [github](https://github.com/pmoulik)
* **Ved Lekic** - *Co-administrator* - [github](https://github.com/vedlekic)

See also the list of [contributors](https://github.com/globalseismology/rem3d/contributors) who participated in this project.

## License

This project is licensed under the GNU GPL v3 or newer - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc
