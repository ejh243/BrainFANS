---
sidebar_position: 6
title: "Step 6: Finalise the pull request"
description: Complete a code review, fix merge conflicts and merge
---

# Step 6: Finalise the pull request

## Code review

Now you have created a pull request, someone else (the reviewer that was assigned to the pull request) should complete a code review of your changes. To learn more about code reviews and some best practices, head over to [these pages](/category/introduction-to-code-reviews).

:::tip[Size of pull request]
Ideally your pull request should be in a state where a reviewer will take no more than about 1 hour to complete. In most cases this means that less than about 400 lines will have changed. If more than 400 lines have been changed, consider splitting your pull request into multiple pull requests (it is likely that the pull request is too complex). If a pull request has too many lines to review, it will likely be pushed back (and therefore less likely to be integrated into the main branch).
:::

## Merge branches

Once the code review has been completed and approved, the branches can be merged. At this stage, any merge conflicts should be resolved. This is usually a simple case of choosing which lines to keep from each file.

If this proves difficult, add someone else onto the assignee list of the pull request.

### Squash and merge

Currently we recommend using the 'squash and merge' option when closing a pull request. 

![Screenshot of squash and merge](/development-pipeline/squash-commit.png)

The reason for this is that it makes the commit history on the master branch have less 'commit noise'. Again, ideally your pull request will cover one overarching feature. The feature that is being implemented should be used as the name of the squashed commit (and the individual commits will be automatically added to the commit description). An added bonus of squashing commits is that if another developer is searching through the repository they will be able to 'zoom and enhance' on any particular pull request and from this they will gather more context behind certain file changes. This is particularly useful when an individual commit does not make much sense without an accompanying commit. To ensure that the squashed commit can be traced back easily to the pull request, ensure that the reference to the pull request (#[pull request number]) is in the title of the commit (this is done for you by default).  

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