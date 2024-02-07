---
sidebar_position: 2
title: "Step 2: Get the latest version of BrainFANS"
description: Update or clone the repository to your local machine
---

# Step 2: Get the latest version of BrainFANS

## First time cloning the repository

Most contributions will be a part of an existing development branch. For example, if the issue you are resolving relates to atac pipeline, then you should have the latest version of the `atac-develop` pipeline. Please use the following unix commands to clone a specific branch of the repository to a personal directory of yours:

```bash
# Move to the directory you want the repository to be cloned to
cd path/to/directory

# Create a copy/clone of the online repository found on GitHub
git clone -b <branch-name> https://github.com/ejh243/BrainFANS.git
```

### In the event you already have the branch cloned

If you have already cloned the branch you are working off of in the past, there is no need to clone the repository again, instead you can just move to the branch in question and pull the latest changes.

```bash
# Move to the location of the BrainFANS repository
cd path/to/BrainFANS

# Change current branch to the branch you wish to work off of
git checkout <branch-name>

# Pull the latest changes
git pull origin <branch-name>
```