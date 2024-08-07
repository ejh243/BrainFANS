---
sidebar_position: 5
title: "Step 5: Creating a pull request"
description: Create a pull request on GitHub
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Step 5: Creating a pull request

Now that you have finished resolving the issue, it is time to create a pull request in order to create a formal request to merge your contributions into the main/development branch of the repository.

:::tip[Don't worry]
*If this feels scary, don't worry. Because you made a separate branch in [step three](./Creating-new-branch.md) you are not in danger of messing up the repository.*
:::

## Pushing to the remote repository

To begin, run these commands in the terminal:

```bash
# Get the latest changes to the repository
git fetch

# Merge these changes into your branch
git checkout your-branch-name
git merge origin/branch # Where branch is the one you are working off of
```

The reason we do this is to ensure that merge conflicts can be adressed before a pull request is made. Please visit [this page](/Developer-information/merge-conflicts.md) for a thorough explanation and example of a merge conflict. If you are working with files no one else is touching, the chances of a merge conflict will be low, but this is still a good habit to get into.

Once any merge conflicts that arise from the above have been resolved, you can publish your branch (push your branch to GitHub) with the following command:

```bash
# Push your changes to the remote repository
git push -u origin your-branch-name
```

:::warning[multiple people working on the same branch]
If you are working on the same branch with other people and the branch has already been published (pushed for the first time), please run the following periodically through the development process:

```bash
# Pull the latest changes from the remote repository
git pull --rebase origin your-branch-name

# Push your changes to the remote repository
git push origin your-branch-name
```

Note that sometimes you will get this error:

```text
CONFLICT (content): Merge conflict in filename.extension
error: Failed to merge in the changes.
Patch failed at commit commit-name
```

If you do get this error, enter the following:

```bash
# Discard the rebase
git rebase --abort

# Complete a normal pull
git pull origin your-branch-name

# Work through the merge conflicts that likely caused --rebase to fail

# Push your changes to the remote repository
git push origin your-branch-name
```

This will ensure that you are not causing any merge conflicts for others when you go to push your changes. `git pull` will merge the version online with your local version of the branch. Warnings of merge conflicts will appear in the terminal if they are apparent at this stage. Fix these before pushing.
:::

## Return to github to finalise the pull request
<Tabs>
  <TabItem value="contributor" label="You are a contributor" default>
    Move to the pull request tab on GitHub:

    ![Screenshot of the location of pull requests](/development-pipeline/pull-request-location.png)

    Click on 'New pull request':

    ![Screenshot of the location of pull request button](/development-pipeline/pull-request-button.png)

    Select the development branch you are working off of for the **base field** and your newly created branch for the **compare field**:

    ![Screenshot of the compare and base](/development-pipeline/pull-request-branch-selection.png)

    Click on 'Create pull request':

    ![Screenshot of the create pull request button](/development-pipeline/create-pull-request-button.png)
  </TabItem>
  <TabItem value="Non-contibutor" label="You are not a contributor">
    If you are not a contributor you just need to be a little more careful
    that your pull request is indeed going to be merged into the master branch
    of ejh243/BrainFANS (and not *your* fork's master branch).

    Move to the pull request tab of your fork on GitHub:

    ![Screenshot of the location of pull requests](/development-pipeline/pull-request-location.png)

    Click on 'New pull request'

    ![Screenshot of the location of pull request button](/development-pipeline/pull-request-button.png)

    Select the development branch you are working off of for the **base field** and your newly created branch for the **compare field**.
    Additionallly, ensure that 'base repository' is ejh243/BrainFANS and 'head repository' is pointing to your own fork.

    ![Screenshot of the compare and base](/development-pipeline/pull-request-branch-selection-fork.png)

    Click on 'Create pull request':

    ![Screenshot of the create pull request button](/development-pipeline/create-pull-request-button.png)
  </TabItem>
</Tabs>

## Fill in the pull request template

![Screenshot of the pull request template](/development-pipeline/pull-request-template-boxes.png)

**PLEASE COMPLETE THE FOLLOWING:**

- Give your pull request a title (What is the main feature this pull request brings)
- Fill in the description with what you have done (be specific)
- Replace [number] with the numbers of any issues that the pull request resolves
- Convert `[ ]`to `[x]` for any checkboxes that are relevant to the pull request

## Assign people and labels to the request

Reviewers -> Assigned people will be notified to review your changes and provide feedback (code review). Your pull request cannot be accepted until a reviewer has approved your commits.
\
Assignees -> Assigned people will be responsible for overseeing the pull request and the merge process, we recommend at least assigning yourself
\
Labels -> Relevant tags for the pull request (best to pick the same ones as used in the issues being addressed)

![Screenshot of applying labels](/development-pipeline/assign-labels-pull-request.png)

## Click 'Create pull request'

![Screenshot of create pull request button](/development-pipeline/create-pull-request.png)
