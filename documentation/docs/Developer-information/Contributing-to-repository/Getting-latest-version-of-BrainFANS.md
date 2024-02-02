---
sidebar_position: 2
title: "Step 2: Get the latest version of BrainFANS"
description: Update or clone the repository to your local machine
---

# Step 2: Get the latest version of BrainFANS

## Cloning the repository

If this is your first contribution, please use the following unix commands to clone the repository to a personal directory of yours:

```console
# Move to the directory you want the repository to be cloned to
cd path/to/directory

# Create a copy/clone of the online repository found on GitHub
git clone https://github.com/ejh243/BrainFANS.git
```

## Getting the latest version of the repository

For subsequent contributions, navigate to the location of your **personal** BrainFANS repository and pull the latest changes.

Most contributions will be a part of an existing development branch. For example, if the issue you are resolving relates to atac pipeline, then you should have the latest version of the `atac-develop` pipeline. In order to avoid merge conflicts, pull the latest changes for the branch that you plan to work off of with these unix commands:

```console
# Change directory to location of your local BrainFANS repository
cd path/to/BrainFANS

# Updates your local BrainFANS repository development branch
# to the latest version that exists on GitHub
git pull origin <development-branch-name>
```