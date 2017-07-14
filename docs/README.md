# REM3D - A toolkit for analyzing and interpreting reference Earth models and data sets.

<img src="docs/rem3dlogo.png" width="256">

## About

REM3D is a Python library for analyzing and interpreting reference Earth models and data sets. The development version has modules in MATLAB, C and Fortran as well.

Homepage:   [rem3d.org](https://maurya.umd.edu)

Documentation: readthedocs [html](http://rem3d.readthedocs.io), [pdf](https://media.readthedocs.org/pdf/rem3d/latest/rem3d.pdf)  

Source code: [github](https://github.com/pmoulik/rem3d) 

Contact the REM3D team at info@rem3d.org with any questions or suggestions.

## Requirements

REM3D needs the following python modules installed in the system.
* Python 2.7.x or Python 3.4+
* Python modules:
  NumPy, SciPy, Matplotlib
  
Users may also want to install MATLAB for some codes. We have tested our versions on MATLAB 2015.

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Basic Installation

1. Install the [Anaconda Python Distribution]( https://www.continuum.io/downloads). We recommend the 64-bit Python 2.7 version. 
2. Install REM3D module using pip by opening a terminal window and entering
`sudo pip install rem3d` 

## Start Here

To begin, the user may want to look at these examples to begin to understand
what tools are available in REM3D and how values are calculated. Below is a
suggested order of examples that begin by introducing each of the user inputs
possible as well as each of the helpers involved with each example.

## About scripting in Python

REM3D has the advantage of being adaptable and extendable in easy scripts. The downside might be that we do not
provide a graphical user interface. For those of you who are not familiar  with python, we suspect it will still be 
relatively easy to adapt the scripts for computations and plotting. 
Here are some specific features and pitfalls on Python:

* Python uses specific indentation. A script might fail if a code block is not indented correctly. We use four spaces and no tabs, 
  mixing these can give trouble.
* Indices require square brackets and function or method calls parentheses (mainly different from Matlab).
* The first index of an array is 0 (e.g. x[0])
* Put dots after numbers to make them floats instead of integers (e.g. 5/3 will give 1 (Python 2.x rounds downward), while 5./3. will give 1.66666666667)

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

### Contributor Guide

REM3D developers may want to install the following softwares and extensions to facilitate proper code usage and documentation.
* sphinx
* numpydoc (sphinx extension)
* sphinxcontrib.bibtex (sphinx extension)
* sphinxcontrib.matlab (sphinx extension)

To install these, please enter the following in a terminal:
```
sudo pip install sphinx numpydoc
sudo pip install sphinxcontrib-bibtex sphinxcontrib-matlabdomain
```

### Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/pmoulik/rem3d/tags). 

### Authors

* **Pritwiraj Moulik** - *Primary administrator* - [github](https://github.com/pmoulik)
* **Ved Lekic** - *Co-administrator* - [github](https://github.com/vedlekic)

See also the list of [contributors](https://github.com/pmoulik/rem3d/contributors) who participated in this project.

## License

This project is licensed under the GNU GPL v3 or newer - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

* Hat tip to anyone who's code was used
* Inspiration
* etc

