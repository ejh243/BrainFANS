---
sidebar_position: 4
title: "Step 4: Resolving the issue"
description: Making changes to the code to resolve the issue
---

# Step 4: Resolving the issue

## Setup

:::warning[Check this]
Before making any commits, run `git status` to check you are indeed in your development branch.
:::

![Screenshot of the output of git status](/development-pipeline/branch-check.png)

If this is not the branch you created in [step three](./Creating-new-branch.md), please run the following:

```console
# Switch to the development branch you created in step three
git checkout <your-branch-name>

# Check that you are indeed on your branch (!! this is very important !!)
git status
```

:::danger[Important] 
If you are on the wrong branch when making commits you run the risk of creating merge conflicts that can be difficult to resolve
:::

If you have skipped ahead and already made commits, please type the following into the terminal:

```console
# Clears any commits you have made (this does not delete your hard work)
git reset --soft HEAD
```

## Making commits

Edit the relevant files in your favourite text editor one at a time.

Once a **small** and complete change has been made use:

```console
# Stage the relevant file to be committed
git add path/to/file

# Commit the changes
git commit -m "Enter your commit message here"
```

To make the above easier for you, get into the habit of using `git status` often (think of it as the git version of `ls`). Using `git status` will show you which files have been changed and the path to said files (which you can then copy and paste to the terminal).

:::tip[Good rules of thumb for commit messages]
1) The message should complete the sentence: "If applied, this commit will [commit-message]"
2) Don't use the word **and**, be succinct. If you only commit small changes this will be much easier
:::

Sometimes it makes sense to make changes to multiple files at the same time (especially if they share the same issue). Don't commit loads of files at the same time out of laziness however, that defeats the point of version control.