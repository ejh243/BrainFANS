---
sidebar_position: 6
title: "Step 6: Finalise the pull request"
description: Complete a code review, fix merge conflicts and merge
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Step 6: Finalise the pull request

## Code review

Now you have created a pull request, someone else (the reviewer that was assigned to the pull request) should complete a code review of your changes. To learn more about code reviews and some best practices, head over to [these pages](/category/introduction-to-code-reviews).

:::tip[Size of pull request]
Ideally your pull request should be in a state where a reviewer will take no more than about 1 hour to complete. In most cases this means that less than about 400 lines will have changed. If more than 400 lines have been changed, consider splitting your pull request into multiple pull requests (it is likely that the pull request is too complex). If a pull request has too many lines to review, it will likely be pushed back (and therefore less likely to be integrated into the main branch).
:::

## Merge branches

Once the code review has been completed and approved, the branches can be merged. At this stage, any merge conflicts should be resolved. This is usually a simple case of choosing which lines to keep from each file. For an in-depth look at merge conflicts, consult [these pages](/Developer-information/merge-conflicts.md)

If this proves difficult, add someone else onto the assignee list of the pull request.

### Types of merges

Depending on the branch that is being merged, there are different types of merges that we recommend. Simply using 'merge' can lead to a lot of *'commit noise'* in the repository over time (making it difficult to see which commits are related to one another). 

If you plan on deleting your development branch immediately after merging, use 'sqaush and merge'. However, if you suspect that the branch will be required again in the future (which will be the case with branches like *web-development* or *atac-develop*), then use a simple merge. There is also a third option available in 'rebase and merge'. To keep things simple, we currently do not recommend using rebase and merge (though you can still learn about it in the tabs below). For the majority of cases, if you are working on a feature branch, squash and merge is recommended. To learn more about these types of merges, you can use the tabs below.

:::danger[Squash merges on long standing branches]
The most glaring problem with squash merges is when a developer continues to use their feature branch after merging. If this is you, save yourself the hassle and **stop now**. 

Why?
\
When you use squash and merge, the timeline of the git repository is pushed into a single point (a single commit). As such, when you next go to merge your changes you will be greeted with a merge conflict **for every single file you have changed** (which is rather unpleasant). The reason for this is that git does not know which file (the one on the base branch or the one on the comparison branch) preceeds the other. As a result, it is up to the developer to decide which lines from each file should be kept (which can be an incredibly lengthy process).
:::

<Tabs>
  <TabItem value="Merge" label="Merge">
    This is the default option for merging branches with git. When you click this button, every commit you make gets effectively moved over to the base branch. This is recommended for branches that you do not want to delete. 
    
    If there are a lot of small, inconsequential commits in the pull request, then this type of merge can lead to a cluttered repository. If another developer inspects commits that have been simply merged into the base branch, it can be difficult to gather all of the context behind said commit. This is because it could be surrounded by hundreds (or thousands) of other commits that aren't related to it in anyway. 
    Use this merge option sparingly to avoid disorganisation in the repository.  
  </TabItem>
  <TabItem value="Squash" label="Squash and merge" default>
    As mentioned above, this is our recommendation for feature branches (branches you expect to delete upon merge) as it removes a lot of *'commit noise'* in the repository. 
    
    Ideally your pull request will cover one overarching feature/enhancement. The feature/enhancement that is brought with the pull request should be used as the name of the squashed commit (and the individual commits will be automatically added to the commit description). An added bonus of squashing commits is that if another developer is searching through the repository they will be able to 'zoom and enhance' on any particular pull request (from this they should be able to gather more context behind certain file changes). This is particularly useful when an individual commit does not make much sense without an accompanying commit. 
    
    To ensure that the squashed commit can be traced back easily to the pull request, ensure that the reference to the pull request (#[pull request number]) is in the title of the commit (this is done for you by default).

    To do a squash merge, use this button:

    ![Screenshot of squash and merge](/development-pipeline/squash-commit.png)
  </TabItem>
  <TabItem value="Rebase" label="Rebase and merge">
    This moves the comparison branch to a new base commit, which basically rewrites the commit history of the base branch. This can be useful for maintaining a more linear commit history, which in turn can be very useful when working with integrating changes from a main branch into a long standing feature branch.

    The reason we don't recommend this is due to the potential of losing commits (which can lead to confusion of other developers later down the line). Also, if any merge conflicts do indeed arise due to this type of merge, it can be very difficult to resolve them, especially if there are a large number of commits involved in such merge conflicts.
    
    Only use rebase and merge if you know what you are doing and other developers are aware of the possible side effects.
  </TabItem>
</Tabs>



## Delete your branch

In most cases your branch is no longer required, GitHub should prompt the user to delete the branch safely after the merge has been completed.

![Screenshot of delete branch button](/development-pipeline/delete-branch-button.png)

If you created a new development branch in [step 3](./Creating-new-branch.md), this will not apply to you. If unsure, notify an admin.

Clicking the 'delete branch' button will only delete the branch on GitHub, not the branch on your copy of the repository. If you want to delete your branch locally, please complete the following:

```bash
# Ensure that there are no uncommited changes to the branch
git checkout <your-branch-name>
git status

# If the last line of the output of git status is 
# 'nothing to commit, working tree clean', you can continue safely

# Delete the branch
git checkout master
git branch -d <your-branch-name>

# Verify the deletion (your branch should not appear in the output)
git branch
```

## Mark related issues as complete

Go back to the issues that were resolved with this pull request and mark them as completed/closed. This is necessary as it cleans up the issues page, making it easier to navigate.

:::info[Automation]
GitHub allows issues to be autocompleted if the issues are linked to the pull request (located in the developments section of the pull request).
:::

![Screenshot of closing an issue](/development-pipeline/close-issue.png)

## Celebrate the completion of your contribution

Congratulations! You have completed a contribution. As you contribute more, this process will become second nature.

:::tip[Future contributions]
If this was your first contribution, we recommend you continue to use this guide for your next few contributions as well.
:::