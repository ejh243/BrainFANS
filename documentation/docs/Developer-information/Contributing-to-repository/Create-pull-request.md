---
sidebar_position: 5
title: "Step 5: Creating a pull request"
description: Create a pull request on GitHub
---

# Step 5: Creating a pull request

Now that you have finished resolving the issue, it is time to create a pull request in order to create a formal request to merge your contributions into the main/development branch of the repository.

:::tip[Don't worry]
*If this feels scary, don't worry. Because you made a separate branch in [step three](./Creating-new-branch.md) you are not in danger of messing up the repository.*
:::

## Pushing to the remote repository

To begin, run this command in the terminal:

```bash
# Push your changes to the remote repository
git push origin <your-branch-name>
```

:::warning[multiple people working on the same branch]
If you are working on the same branch with other people and the branch has already been published (pushed for the first time) to the online repository, please run the following instead:

```bash
# Pull the latest changes from the remote repository
git pull origin <your-branch-name>

# Push your changes to the remote repository
git push origin <your-branch-name>
```

This will ensure that you are not causing any merge conflicts when you go to push your changes. `git pull` will merge the version online with your local version of the branch. Warnings of merge conflicts will appear in the terminal if they are apparent at this stage. Fix these before pushing.
:::

## Return to github to finalise the pull request

Move to the pull request tab on GitHub:

![Screenshot of the location of pull requests](/development-pipeline/pull-request-location.png)

Click on 'New pull request':

![Screenshot of the location of pull request button](/development-pipeline/pull-request-button.png)

Select the development branch you are working off of for the **base field** and your newly created branch for the **compare field**:

![Screenshot of the compare and base](/development-pipeline/pull-request-branch-selection.png)

Click on 'Create pull request':

![Screenshot of the create pull request button](/development-pipeline/create-pull-request-button.png)

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