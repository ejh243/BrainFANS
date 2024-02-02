---
sidebar_position: 2
title: "Step 2: Creating a new branch"
description: Creating a new branch off of the web-development branch
---

# Step 2: Creating a new branch

:::info[Cloning the repository]
If you do not already have the BrainFANS repository on your local machine, please refer to the [existing documentation](/Developer-information/Contributing-to-repository/Getting-latest-version-of-BrainFANS.md#cloning-the-repository) on this before continuing
:::

In the vast majority of cases, documentation changes will be independent of one another. This helps out with merge conflicts and so the importance on branch creation is not as high as with contributions to the pipelines in the repository.

Even with this in mind, we still recommend creating a new branch off of the web-development branch so that changes to the documentation (especially if the changes are [larger in scope](./Getting-started.md#larger-scope)) are easy to track and follow.

To create a new branch off of the web-development branch, head over to the terminal and use the following unix commands:

```console
# Change your directory to the position of your local BrainFANS repository
cd /path/to/BrainFANS

# Pull the latest changes to the web-development branch
git pull origin web-development

# Create a new branch off of web-development and move to it
git checkout -b web-development <your-branch-name>
```

Now you are on the new branch you can start making commits.
