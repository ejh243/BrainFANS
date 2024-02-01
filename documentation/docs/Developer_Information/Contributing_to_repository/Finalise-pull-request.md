---
sidebar_position: 6
title: "Step Six: Finalise the pull request"
description: Complete a code review, fix merge conflicts and merge
---

# Step Six: Finalise the pull request

## Code review

Now you have created a pull request, someone else (the reviewer that was assigned to the pull request) should complete a code review of your changes. A separate document will be available in the future detailing how this process is expected to play out. For now, please consult [Google's code review guide](https://google.github.io/eng-practices/review/reviewer/).

## Merge branches

Once the code review has been completed, the branches can be merged. At this stage, any merge conflicts should be resolved. This is usually a simple case of choosing which lines to keep from each file.

If this proves difficult, add someone else onto the assignee list of the pull request.

## Delete your branch

In most cases your branch is no longer required, GitHub should prompt the user to delete the branch safely after the merge has been completed.

![Screenshot of delete branch button](/development-pipeline/delete-branch-button.png)

If you created a new development branch in [step 3](./Creating-new-branch.md), this will not apply to you. If unsure, notify an admin.

Clicking the 'delete branch' button will only delete the branch on GitHub, not the branch on your copy of the repository. In order to delete your branch locally, please complete the following:

```console
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

:::info[Clean up]
Go back to the issues that were resolved with this pull request and mark them as completed/closed. This is necessary as it cleans up the issues page, making it easier to navigate.
:::

![Screenshot of closing an issue](/development-pipeline/close-issue.png)

## Celebrate the completion of your contribution

Congratulations! You have completed a contribution. As you contribute more, this process will become second nature.

:::tip[Future contributions]
If this was your first contribution, we recommend you continue to use this guide for your next few contributions as well
:::