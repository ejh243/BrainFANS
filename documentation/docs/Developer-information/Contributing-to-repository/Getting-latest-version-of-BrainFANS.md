---
sidebar_position: 2
title: "Step 2: Get the latest version of BrainFANS"
description: Update or clone the repository to your local machine
---

# Step 2: Get the latest version of BrainFANS

If this is your first contribution, please use the following unix commands to clone the repository to a personal directory of yours:

```console
# Move to the directory you want the repository to be cloned to
cd path/to/directory

# Create a copy/clone of the online repository found on GitHub
git clone https://github.com/ejh243/BrainFANS.git
```

For subsequent contributions, navigate to the location of your **personal** BrainFANS repository and pull the latest changes:

```console
cd path/to/BrainFANS

# Updates your local BrainFANS to the latest version in the GitHub repository
git pull origin master
```

Most contributions will be a part of an existing development branch. For example, if the issue you are resolving relates to atac pipeline, then you should have the latest version of the `atac-develop` pipeline.

If you are working on a specific branch of BrainFANS, you will need to pull that branch to your local repository:

```console
cd path/to/BrainFANS

# Updates your local BrainFANS to the latest version in the GitHub repository
git pull origin <branch-name>
```