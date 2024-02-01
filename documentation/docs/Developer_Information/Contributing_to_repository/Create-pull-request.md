---
sidebar_position: 4
title: "Step Five: Creating a pull request"
description: Create a pull request on GitHub
---

# Step Five: Creating a pull request

Now that you have finished resolving the issue, it is time to create a pull request in order to create a formal request to merge your contributions into the main/development branch of the repository.

*If this feels scary, don't worry. Because you made a separate branch in [step three](./Creating-new-branch.md) you are not in danger of messing up the repository.*

## Pull and push

To begin, run these two commands sequentially:

```console
# Pull any changes from the remote repository that may have changed since step one
# This is recommeded to protect you from finding merge conflicts later on
git pull origin <your-branch-name>

# Push your changes to the remote repository
git push origin <your-branch-name>
```

Using `git pull` and then `git push` might be unecessary in some circumstances. But it is a good habit to get into, especially as more people work on the repository (it will save you time later, trust me).

## Return to github to finalise the pull request

Move to the pull request tab on GitHub

![Screenshot of the location of pull requests](/development-pipeline/pull-request-location.png)

Click on 'New pull request'

![Screenshot of the location of pull request button](/development-pipeline/pull-request-button.png)

Select 'master' for the **base field** and [your development branch] for the **compare field**

![Screenshot of the compare and base](/development-pipeline/pull-request-branch-selection.png)

Click on 'Create pull request'

![Screenshot of the create pull request button](/development-pipeline/create-pull-request-button.png)

## Fill in the pull request template

![Screenshot of the pull request template](/development-pipeline/pull-request-template-boxes.png)

**PLEASE COMPLETE THE FOLLOWING:**

- Give your pull request a title
- Fill in the description with what you have done
- Replace [number] with the number of the issue that the pull request resolves
- Convert `[ ]`to `[x]` for any checkboxes that are relevant to the pull request

## Assign people and labels to the request

Reviewers -> Assigned people will be notified to review your changes and provide feedback (code review)
\
Assignees -> Assigned people will be responsible for overseeing the pull request and the merge process, we recommend at least assigning yourself
\
Labels -> Relevant tags for the pull request (best to pick the same ones as used in the issue being addressed)

![Screenshot of applying labels](/development-pipeline/assign-labels-pull-request.png)

## Click 'Create pull request'

![Screenshot of create pull request button](/development-pipeline/create-pull-request.png)