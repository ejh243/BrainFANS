---
sidebar_position: 3
title: "Step 3: Create a new branch"
description: Create a new sub branch off of existing development branches
---

# Step 3: Create a new branch

Create a new feature branch in the terminal with a name that is relevant 
to the issue that you are resolving. For most changes, you'll likely
be creating a new branch that stems from the master branch. However, if you
plan on working on an existing development branch (perhaps some pipeline
that has yet to be merged into master), you'll need to stem off of that
branch instead. 

The important thing to note here is that your branch will start off in the
same state as the branch that you stem from. If you create the branch off of
a really old branch that hasn't been touched in years, your branch will also
be out of date.

:::danger[Important]
If you do not follow this step you will not be able to push any of your 
changes due to the branch protection rules on master branche.
:::


```bash
# Creates new branch off of the master branch
# and automatically moves you onto your newly created branch
git checkout -b your-branch-name origin/master 

# Does the same, but creates the new branch off of a branch of your choosing 
git checkout -b your-branch-name origin/development-branch
```
