Public Release Workflow
-----------------------

We maintain 3 major branches for our client libraries and these are relevant to public releases. Read/write access to these are restricted to main administrators and developers:  
* `devel` — Active development occurs on this branch or in branches that spin off from this.  
* `public` — Development for bug fixes happens here. We also bump versions and update the changelog on this branch.  
* `master` — We squash commits from the release branch into single release commits on this branch as well as tagging releases.  

New branches may be created for individual projects. Please clone the `devel` branch to build upon the latest codes  
`git checkout -b new_branch devel`  
You can push this locally created branch to the remote `globalseismology/avni` with  
`git push -u origin new_branch`  

The `master` branch on `https://github.com/globalseismology/avni` tracks the publicly open repository `https://github.com/geodynamics/avni`. The workflow for every public release involves getting all edits done in `public` based on a subset of codes selected from `devel`. Once `public` is tested and ready, an admin pushes it downstream to `master`, squashes all comments and tags it:  
`git checkout master`  
`git merge --squash public`  
`git commit -m "1.0.0"`  
`git tag 1.0.0 -m "1.0.0"`  

Then the admin pulls the latest changes and syncs with the public CIG repository using a series of commands such as below:  
`git remote add cig git@github.com:geodynamics/avni.git`  
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