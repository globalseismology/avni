Public Release Considerations
-----------------------------

There are several considerations for a public release:

* We want to maintain public/private demarcation  of code at all levels - subdirectory, file and specific subroutines. In other words, we should be able to specify which parts of the codebase we want users to have public access to at a granular level. We accomplish this by using a [public file](.public) that lists the codebase components that need to be retained for a public release. Within each file, we also demarcate subroutines that are private by listing them within private code blocks
```
######### private start #########
subroutine xyz():
######### private end #########
```

* We want to squash all commits made in the private repository to retain privacy of the development process. At the same time, we wish to retain SHA-256 hash of the steps to allow investigation of repository commit history. Both may be accomplished by isolating the squashing to a downstream branch (`public`) and using merge commits. Each merge commit makes a new hash that is common to the two branches.

* We wish to encourage user involvement on the public repository through Github [issues](https://github.com/globalseismology/avni/issues) and pull requests of features developed on top of the releases. These improvements can also be leveraged in projects that use our private repositories. This is done by syncing any new changes from public repository through a shared `public` branch.

Layout of branches
------------------

We maintain 3 major branches for our client libraries that are relevant to public releases. Read/write access to these are restricted to main administrators and developers:
* `main` — Active development occurs on this branch or in branches that spin off from this.
* `public` — Development for bug fixes happens here. We also bump versions and update the changelog on this branch. 
* `release` — We squash commits from the public branch into single release commits on this branch as well as tagging releases. 

Separate `public` and `release` branches are needed so that all the testing, removing private files and merging from `main` can be removed within the isolated `'public` branch.  These branches were created as `git checkout -b public 6aa31bf` and `git checkout -b release 6aa31bf`, based on the first SHA commit of the main branch (`6aa31bf`), which was found using `git log --abbrev-commit --pretty=oneline`. In order for merging to happen seamlessly, this procedure of accounting for commit history was needed to preserve the record.

New branches may be created for individual projects. Please clone the `main` branch to build upon the latest codes
`git checkout -b new_branch main`
You can push this locally created branch to the remote `globalseismology/avni` with
`git push -u origin new_branch`

Repository Public Release Workflow
-----------------------

The `public` branch on [avni-private repository](https://github.com/globalseismology/avni-private) tracks the publicly open [avni repository](https://github.com/globalseismology/avni). The workflow for every public release involves getting all edits done in `public` based on a subset of codes selected from `main`. We need to write a shell script that removes the files not in [public file](.public) and parts of the files that are marked as private.

Choose a version  number in the file `avni/version.py`. Once `public` is tested and ready, an admin squashes all comments to the `release` branch and tags it as the same version as in the `avni/version.py`:
`git checkout release`  
`git merge --squash public`  
`git commit -m "v1.0.0"`  
`git tag v1.0.0 -m "v1.0.0"`  
`git push origin v1.0.0`  

Then the admin pulls the latest changes from and syncs with the [public repository](https://github.com/globalseismology/avni) using a series of commands such as below:
`git remote add cig git@github.com:globalseismology/avni.git`  
`git pull cig main`  
`git push cig HEAD:main`  
`git push cig v1.0.0`  

We also need to push these changes to the branches on origin and merge the squashed commit back to `main`. You may suspect that git would be confused merging a squashed commit back into branches containing the non-collapsed commits, but it all works just as expected. Git is smart enough to realize no changes need to be made when merging in the squashed commit, but we should still merge to keep our branches in sync.
`git push origin release`  
`git checkout public`  
`git merge release`  
`git push origin public`  

Publishing on PyPi
-----------------------

Now that the code for the python package is almost complete, you can start building the distribution archives. Archives are compressed files that help your package to be deployed across multiple platforms and also make it platform independent. To get started, navigate to PyPI at `https://pypi.org/` to register for an account; we’ll use TestPyPI to test that our release looks okay and works as expected, so also navigate to TestPyPI at `https://test.pypi.org/` to register for an account on TestPyPI (TestPyPI is entirely separate from PyPI, so your PyPI account is separate from your TestPyPI account). Interacting with PyPI and TestPyPI is done, of course, using a Python package called twine, which you can install with

`pip install twine`  

At first, we will simply package the source of our package and upload it to PyPI as a first release. To create a source distribution, do

`python setup.py sdist`


Running this command creates a source distribution for the current released version in a `dist/` directory, with a filename of something like `dist/PACKAGENAME-VERSION.tar.gz`, e.g., `dist/avni-0.1.tar.gz` for the first release of this package. If this is not your first release, you’ll want to remove old releases from the `dist/` directory before continuing (but you should be working with a fresh clone, so an empty dist/ directory).

The first thing you want to check is that the `long_description` that you have specified in your setup.py’s setuptools.setup command can be correctly displayed by PyPI; to check this, run

`twine check dist/*`

which performs this check for all source distributions in the `dist/` directory (you can also manually specify the latest one as the filename). It’s useful to run this check before creating the final release, so you can fix any issues before tagging the release. 

Next, upload the source distribution to TestPyPI. TestPyPI is a clone of PyPI and you can therefore check that your package release upload and the way it appears on the website look okay before publishing the final release. Once you upload a version to PyPI you can no longer easily change it, so for every release you should check with TestPyPI to make sure you are not making any mistakes (e.g., formatting errors in the long_description that will form the webpage for your package on PyPI. You can upload your source distribution by doing

`twine upload --repository-url https://test.pypi.org/legacy/ dist/*`

which uploads every release file in the `dist/` directory

Then navigate to your project’s page on [TestPyPI](https://test.pypi.org/project/avni/) and check that all looks well. You’ll want to check that (a) your `long_description`, which is typically the README.md or README.rst that you have on GitHub, is rendered correctly as the webpage of your release, and (b) that the source distribution was properly uploaded by checking the `Download files` tab. 

This now looks as we expect it to look, similar to the GitHub page. TestPyPI is a fully functional Python package index, so you can install with pip from TestPyPI using, e.g.

`pip install -i https://test.pypi.org/simple/ avni`
which is useful to check your package installs without a hitch using pip.

Once you’re happy with the way your package release looks on TestPyPI, it’s time to upload your release to PyPI itself. This you do with

`twine upload dist/*`

This works the same way as uploading to TestPyPI. Once the upload is finished, you can navigate to your package’s PyPI site, e.g. `https://pypi.org/project/avni/`.
