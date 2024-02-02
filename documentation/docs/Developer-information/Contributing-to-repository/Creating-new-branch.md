---
sidebar_position: 3
title: "Step 3: Create a new branch"
description: Create a new sub branch off of existing development branches
---

# Step 3: Create a new branch

Create a new development branch in the terminal with a name that is relevant to the issue that you are resolving. In general, you should be making a branch off of a development branch. If the issue you are resolving is very general, small or unrelated to development branches, you may create your branch off of the master branch. 

:::danger[Important]
If you do not follow this step you will not be able to push any of your changes or run the risk of creating merge conflicts (which can be a real headache)
:::


```console
# Creates new branch off of the provided development branch
# and automatically moves you onto this branch
git checkout -b <development-branch-name> <your-branch-name>
```