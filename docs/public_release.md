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
* `public` — Development for bug fixes happens here. We also bump versions and update the changelog on this branch. We squash commits from the release branch into single release commits on this branch as well as tagging releases.

New branches may be created for individual projects. Please clone the `main` branch to build upon the latest codes
`git checkout -b new_branch main`
You can push this locally created branch to the remote `globalseismology/avni` with
`git push -u origin new_branch`

Public Release Workflow
-----------------------

The `public` branch on [avni-private repository](https://github.com/globalseismology/avni-private) tracks the publicly open [avni repository](https://github.com/globalseismology/avni). The workflow for every public release involves getting all edits done in `public` based on a subset of codes selected from `main`. We need to write a shell script that removes the files not in [public file](.public) and parts of the files that are marked as private.

Once `public` is tested and ready, an admin squashes all comments to a `'release-v1.0.0` branch and tags it:
`git checkout release-v1.0.0`
`git merge --squash public`
`git commit -m "1.0.0"`
`git tag 1.0.0 -m "1.0.0"`

Then the admin pulls the latest changes from and syncs with the [public repository](https://github.com/globalseismology/avni) using a series of commands such as below:
`git remote add cig git@github.com:globalseismology/avni.git`
`git pull cig master`
`git push cig HEAD:master`

We also need to push these changes to the branches on origin and merge the squashed commit back to `public` and `devel`. You may suspect that git would be confused merging a squashed commit back into branches containing the non-collapsed commits, but it all works just as expected. Git is smart enough to realize no changes need to be made when merging in the squashed commit, but we should still merge to keep our branches in sync.
`git push origin master`
`git checkout public`
`git merge master`
`git push origin public`
`git checkout devel`
`git merge public`
`git push origin devel`