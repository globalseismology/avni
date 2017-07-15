Welcome to the development wiki of REM3D.  The goal of this page is to provide instructions and best practices to guide developers and contributing scientists. It is intended as a developer reference, not as a replacement for the user manual.

Contributor Guide
=================

TPd2VQ developers may want to install the following softwares and extensions to facilitate proper code usage and documentation.
* sphinx
* numpydoc (sphinx extension)
* sphinxcontrib.bibtex for references (sphinx extension)
* sphinxcontrib.matlab for matlab codes (sphinx extension)
* recommonmark (writing in Markdown *.md rather than reStructured Text *.rst)

To install these, please enter the following in a terminal:
```
sudo pip install sphinx numpydoc
sudo pip install sphinxcontrib-bibtex sphinxcontrib-matlabdomain
sudo pip install recommonmark
```


Contributors should always work in the 'devel' branch. To clone that directory, type the following in a terminal.
```
git clone -b devel https://username:password@github.com/pmoulik/rem3d.git
```
Pull requests to the protected 'master' branch from the 'devel' branch will be attempted at regular intervals.

Please document code in the standard format (e.g. with three " in python routines) so that sphinx can auto build documentation in the ./doc folder. To make a local copy of the documentation, run ./doc/make_docs.sh on the command line.

Further reading
---------------

- [Software Carpentry Git Tutorial](https://swcarpentry.github.io/git-novice/index.html): For those who are not familiar with Git, we recommend attempting this tutorial, or even better, attending an in-person tutorial session if available in your area.

- [Best Practices](best_practices.md):  Accepted conventions for REM3D development are described here.  To create a product that is maintainable over the long term, it is important that contributing scientists follow these conventions.

- [Using Git for REM3D](git_for_REM3D.md):  The standard way of contributing to REM3D development is described here. This approach is recommended for most users, especially those not familiar with Git. These instructions are adequate for contributing to REM3D, but we recommend that all new users to Git attempt the tutorial above for more complete understanding of the workflow.

- [Versioning Conventions](versioning_conventions.md):  The conventions behind REM3D's version numbering system are explained here.

