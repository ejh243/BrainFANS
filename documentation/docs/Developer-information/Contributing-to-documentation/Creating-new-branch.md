---
sidebar_position: 2
title: "Step 2: Creating a new branch"
description: Creating a new branch off of the web-development branch
---

# Step 2: Creating a new branch

:::info[Cloning the repository]
If you do not already have the BrainFANS repository on your local machine, please refer to the [existing documentation](/Developer-information/Contributing-to-repository/Getting-latest-version-of-BrainFANS.md) on this before continuing.
:::

Start by creating a new branch off of the **web-development** (not master) branch so that changes to the documentation are easy to track and follow.

To create a new branch off of the web-development branch, head over to the terminal and use the following git commands:

```bash
# Change your directory to the position of your local BrainFANS repository
cd /path/to/BrainFANS

# Fetch the latest changes to the repository
git fetch

# Create a new branch off of web-development and move to it
git checkout -b <your-branch-name> origin/web-development
```

Now you are on the new branch you can start making commits.
