**Table of Contents**

- [One time configurations for all repositories](#one-time-configurations-for-all-repositories)
  - [Open a free account at GitHub](#open-a-free-account-at-github)
  - [Configure ssh-keys for GitHub](#configure-ssh-keys-for-github)
  - [Make sure your version of Git is fine](#make-sure-your-version-of-git-is-fine)
  - [Install hub](#install-hub)
  - [Download our configuration script `config_repo`](#download-our-configuration-script-config_repo)
- [One-time configuration for each repository](#one-time-configuration-for-each-repository)
- [Regular daily usage](#regular-daily-usage)


One time configurations for all repositories
============================================

Open a free account at GitHub
-----------------------------

-   If you do not already have one, you will need to open a free account at <https://github.com>. In order to become a contributor / developer for REM3D, you need to create an account on our [Web site](http://avni.org/login/register) and get access rights from us by requesting [here](http://avni.org/join-us/github).

Configure ssh-keys for GitHub
-----------------------------

-   <https://help.github.com/articles/generating-ssh-keys>

Make sure your version of Git is fine
-------------------------------------

You need to have a version of Git greater or equal to 1.8, thus type:

        $ git --version


Install hub
-----------

-   Follow the instructions at https://github.com/github/hub under https://github.com/github/hub#standalone. Basically you just need to download the binary and put it somewhere in your PATH.

-   Permanently alias `git` to `hub` by putting this in your `.bashrc`:

        $ alias git='hub'

Download our configuration script `config_repo`
-----------------------------------------------

-   Copy/paste the following script in `~/bin/config_repo`

        #!/bin/bash
        # Set configuration for github.
        # Configuration is saved in .git/config

        ARGS=1
        if [ $# -ne $ARGS ] || [ $1 == "-help" ]
        then
          echo "usage: config_repo your_github_user_name "
          exit 0;
        fi

        GITHUB_USER=$1

        GITVERSIONMAJOR=`git --version | head -1 | cut -f 3 -d " " | cut -c 1`
        GITVERSIONMINOR=`git --version | head -1 | cut -f 3 -d " " | cut -c 3`

        if [[ "$GITVERSIONMAJOR" -eq 0 || ("$GITVERSIONMAJOR" -eq 1 && "$GITVERSIONMINOR" -le 7) ]]; then
          echo "git version >= 1.8 is needed for the scripts to work, but you have"
          git --version
          echo "please install git version >= 1.8 and try again"
          exit -1
        fi

        git config --global url."https://".insteadOf git://
        git config branch.devel.remote origin
        git config branch.devel.merge devel
        git config branch.devel.pushremote $GITHUB_USER
        git config push.default current
        git remote set-head origin devel

-   Make it executable

        $ chmod u+x ~/bin/config_repo

One-time configuration for each repository
==========================================

-   Clone the repository on your machine

        $ git clone -b devel https://username:password@github.com/geodynamics/avni.git

-   Checkout the `devel` branch

        $ cd avni
        $ git checkout -b devel origin/devel

-   Confirm that you are on the `devel` branch (* next to the name)

        $ git branch

-   Call our configuration script (replace "your_github_name" with your GitHub name, i.e., with your login on the GitHub Web site)

        $ config_repo your_github_name

-   create or update a fork, that is, a copy of the repository on your GitHub account. If you are using several computers (a desktop, your laptop etc.), each with a copy of the code, you need to type this on each machine; the first one will create a copy on GitHub, and all the others on other machines will make your local version aware of the existing copy on GitHub

        $ git fork

Regular daily usage
===================

-   update your copy of the repository

        $ git pull

-   if you get conflicts when doing so (i.e. if local changes you have made conflict with changes made by others on the same line of the same file of the source code), a powerful way of resolving them is to type this: (_meld_ needs to be installed on your system; if it is not, you can install it with _apt-get install meld_ or similar)

        $ git mergetool --tool=meld

-   make some changes to any file you want using your favorite editor (in the line below we use `vi` as an example)

        $ vim some_file.f90

-   commit your changes locally, adding a very short message (one line) explaining what you have changed; it is recommended to do a git pull right before that in order to make sure that your local copy is up-to-date

        $ git pull ; git commit -a -m "Explain your commit"

-   if you get conflicts when committing your changes (i.e. if your changes conflict with changes made by others on the same line of the same file of the source code), a powerful way of resolving them is to type this: (_meld_<http://meldmerge.org> needs to be installed on your system; if it is not, you can install it with _yum install meld_ in Linux, download MacOS version from <https://yousseb.github.io/meld>. )

        $ git mergetool --tool=meld

-   (optional) if you want to check what has changed (and thus what will be committed) before typing the `git commit` above, you can type one or both of these two commands:

        $ git status -s
        $ git diff

-   push your changes to your GitHub fork; it is recommended to do a git pull right before that in order to make sure that your local copy is up-to-date

        $ git pull ; git push

-   Create a pull-request to get your changes into the main repository (this is needed only once for each change; if you are fixing an existing change after receiving an error message from our BuildBot code-consistency checking system, you need the "git push" above again but you do NOT need to create a pull request a second time); it is recommended to do a git pull right before that in order to make sure that your local copy is up-to-date

        $ git pull ; git pull-request

_Note (optional):_ It is not strictly necessary to create a pull request for *every* commit you make if you do not want to, you can safely submit pull requests after making a few commits instead if you prefer. However, it also does not hurt to do so.
