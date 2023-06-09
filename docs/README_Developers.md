
Development Guide
=================

Welcome to the development wiki of AVNI.  The goal of this page is to provide instructions and best practices to guide developers and contributing scientists. It is intended as a developer reference, not as a replacement for the user manual.


Advanced installation
---------------------

Please note the license terms below specific to developers. If you want to use AVNI's development routines, you will need to do the following additional steps. This assumes that you have a free account at <https://github.com>

1. Create an account on our [Web site](http://globalseismology.org/login/register) and request access rights from [here](http://globalseismology.org/join-us/github).
2. Please accept the invitation through e-mail. Clone the development branch from the AVNI git repository through the terminal. This will create a directory that contains a folder called `avni`. Here `username` and `password` are your login credentials for Github.
`git clone -b devel https://username:password@github.com/globalseismology/avni.git`
3. Install AVNI module by opening a terminal window, navigating to `avni` directory and entering
`pip install -e . --user`
This lets you keep working on files inside the Github folder without recompiling the codes.


Get started
-----------

AVNI developers may want to install the following softwares and extensions to facilitate proper code usage and documentation.
* sphinx
* numpydoc (sphinx extension)
* sphinxcontrib.bibtex for references (sphinx extension)
* recommonmark (writing in Markdown *.md rather than reStructured Text *.rst)
* sphinx_rtd_theme (sphinx theme designed to look modern and be mobile-friendly)

This can easily be done using the [environment_develop.yml](conda/environment_devel.yml) to create the `avni_devel` environment by entering the following in a terminal:
```
conda create --yes --name avni_devel -f docs/conda/environment_devel.yml
source $CONDA_ROOT/etc/profile.d/conda.sh
conda activate avni_devel
```
where `CONDA_ROOT` is the location of Anaconda directory (e.g. `/home/user/anaconda3`)

Please document code in the standard format (e.g. with three """Comments """ in python routines) so that sphinx can auto build documentation in the ./doc folder. To make a local copy of the documentation, run ./doc/make_docs.sh on the command line. More detailed guidelines are available in the [best practices](best_practices.md) webpage.


Development Branches
--------------------
We maintain 3 major branches for our client libraries and these are relevant to public releases. Read/write access to these are restricted to main administrators and developers:
* `devel` — Active development occurs on this branch or in branches that spin off from this.
* `public` — Development for bug fixes happens here. We also bump versions and update the changelog on this branch.
* `master` — We squash commits from the release branch into single release commits on this branch as well as tagging releases.

Contributors and core developers should always start with the `devel` branch. To clone the AVNI repository, type the following in a terminal using your Github username and password.
```
git clone -b devel https://username:password@github.com/globalseismology/avni.git
```
Pull requests to the protected `master` branch from the `devel` branch will be made only by *Principal Developers* before public releases. More detailed guidelines are available in the [Git for AVNI](git_for_AVNI.md) webpage.

New branches may be created for individual projects and the relevant team. Please clone the `devel` branch to build upon the latest codes
`git checkout -b new_branch devel`
You can push this locally created branch to the remote `globalseismology/avni` with
`git push -u origin new_branch`


License
-------

The master branch of this project is licensed under the GNU GPL v3 or newer, is synced with the [CIG public repository](https://github.com/geodynamics/avni) and may be installed in any platform - see the [LICENSE](../LICENSE) file for details. All other branches are copyrighted, must remain in our [private Globalseismology repository](https://github.com/globalseismology/avni) and may not be installed in platforms without permission from code contributors and repository administrators. Please write to **avni@globalseismology.org** or see our [FAQ](FAQ.md) for more details.

Further reading
---------------

- [Software Carpentry Git Tutorial](https://swcarpentry.github.io/git-novice/index.html): For those who are not familiar with Git, we recommend attempting this tutorial, or even better, attending an in-person tutorial session if available in your area.

- [Using Git for AVNI](git_for_AVNI.md):  The standard way of contributing to AVNI development is described here. This approach is recommended for most users, especially those not familiar with Git. These instructions are adequate for contributing to AVNI, but we recommend that all new users to Git attempt the tutorial above for more complete understanding of the workflow.

- [Conda Environments](conda/README_conda.md):  AVNI can be installed with various combinations of packages. We have tested the sets of packages listed here which can be installed with [Anaconda Python Distribution](https://www.continuum.io/downloads). Please feel free to update these files or create new ones that work for specific applications.

- [Best Practices](best_practices.md):  Accepted conventions for AVNI development are described here.  To create a product that is maintainable over the long term, it is important that contributing scientists follow these conventions.

- [Versioning Conventions](versioning_conventions.md):  The conventions behind AVNI's version numbering system are explained here.

- [Public Release](public_release.md):  The workflow employed by admins for public releases to CIG repositories.
