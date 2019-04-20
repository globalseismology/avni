Public Release Workflow
-----------------------

The `master` branch tracks the publicly open repository `https://github.com/geodynamics/rem3d`. The workflow for every public release involves getting all edits done in `public` based on a subset of codes selected from `devel`. Once `public` is tested and ready, an admin pushes it downstream to `master`, squashes all comments and tags it:
`git checkout master`  
`git merge --squash public`  
`git commit -m "1.0.0"`  
`git tag 1.0.0 -m "1.0.0"`  
Then the admin pulls the latest changes and syncs with the public CIG repository using a series of commands such as below:
`git remote add cig git@github.com:geodynamics/rem3d.git`  
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