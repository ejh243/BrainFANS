---
sidebar_position: 1
---
# The Development Pipeline

The following document details how to make a change to the [BrainFANS](https://github.com/ejh243/BrainFANS) repository through the command line (and your favourite text editor). Some text editors like VScode and Sublime Text can eliminate the need to remember git commands. We recommend that you use the terminal for your first couple of contributions in order to become more familiar with the version control process.

If you are ever unsure about one of the steps described in this document, please ask Sam Fletcher (s.o.fletcher@exeter.ac.uk) for further clarification.

## Overview of Process

1) [**Get the latest version of BrainFANS**](#step-one-get-the-latest-version-of-brainfans)
2) [**Submit an issue**](#step-two-submit-an-issue)
3) [**Create a new branch**](#step-three-create-a-new-branch)
4) [**Resolve the issue**](#step-four-resolve-the-issue)
5) [**Create the pull request**](#step-five-creating-a-pull-request)
6) [**Finalise the pull request**](#step-six-finalise-the-pull-request)

## Step One: Get the latest version of BrainFANS

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

## Step Two: Submit an issue

Before you make any changes to the repository, we recommend that you create an issue on GitHub. This will make it easier for others to know who is working on what. It is also a great way of documenting  the extent of the changes that you are planning to make and that you are already working on a solution. If you are solving an issue someone else has already logged, this step can be skipped.

### Open the issues tab on GitHub

![Screenshot of location of issues tab](/development-pipeline/issues-location.png)

### Click 'New issue'

![Screenshot of new issue button](/development-pipeline/new-issue-button.png)

### Select the relevant template for your report

- Bug: Make a report for any bugs you have encountered in the repository (*i.e.* Unexpected/unwanted behaviour)
- Code refactor: Make suggestions for how code can be optimised or improved
- Documentation: Report any faults in the documentation
- Enhancement: Suggest any features that should be added/altered

### Fill in the issue report

![Screenshot of filling in issue](/development-pipeline/fill-in-issues.png)

### Assign relevant tags (if applicable)

Add relevant personnel to assignees section, assignees will be notified of the issue. If you are planning on address the issue yourself, remember to assign yourself.  You can also add further labels to the issue if they are relevant (obvious labels related to the type of issue template you have selected will be added for you).

![Screenshot of applying labels](/development-pipeline/assign-labels.png)

### Submit the issue report

![Screenshot of submission button](/development-pipeline/submit-issue-button.png)

## Step Three: Create a new branch

Switch back over to the terminal and create a new development branch with a name that is relevant to the issue that you are resolving:

```console
# Creates new branch and automatically moves you onto this branch
git checkout -b <your-branch-name>
```

## Step Four: Resolve the issue

### Setup



Before making any commits, run `git status` to check you are indeed in your development branch.

![Screenshot of the output of git status](/development-pipeline/branch-check.png)

If this is not the branch you created in [step three](#step-three-create-a-new-branch), please run the following:

```console
# Switch to the development branch you created in step three
git checkout <your-branch-name>

# Check that you are indeed on your branch (!! this is very important !!)
git status
```

### Making commits

Edit the relevant files in your favourite text editor one at a time.

Once a **small** and complete change has been made use:

```console
# Stage the relevant file to be committed
git add path/to/file

# Commit the changes
git commit -m "Enter your commit message here"
```

To make the above easier for you, get into the habit of using `git status` often (think of it as the git version of `ls`). Using `git status` will show you which files have been changed and the path to said files (which you can then copy and paste to the terminal).


There are two good rules of thumb for your commit messages:

1) The message should complete the sentence: "If applied, this commit will [commit-message]"
2) Don't use the word **and**, be succinct. If you only commit small changes this will be much easier

Sometimes it makes sense to make changes to multiple files at the same time (especially if they share the same issue). Don't commit loads of files at the same time out of laziness however, that defeats the point of version control.

## Step Five: Creating a pull request

Once you have finished resolving the issue, it is time to create a pull request in order to merge your contributions into the main branch of the repository.

*If this feels scary, don't worry. Because you made a separate branch in [step three](#step-three-create-a-new-branch) you are not in danger of messing up the repository.*

### Pull and push

To begin, run these two commands sequentially:

```console
# Pull any changes from the remote repository that may have changed since step one
git pull origin <your-branch-name>

# Push your changes to the remote repository
git push origin <your-branch-name>
```

Using `git pull` and then `git push` might be unecessary in some circumstances. But it is a good habit to get into, especially as more people work on the repository (it will save you time later, trust me).

### Return to github to finalise the pull request

Move to the pull request tab on GitHub

![Screenshot of the location of pull requests](/development-pipeline/pull-request-location.png)

Click on 'New pull request'

![Screenshot of the location of pull request button](/development-pipeline/pull-request-button.png)

Select 'master' for the **base field** and [your development branch] for the **compare field**

![Screenshot of the compare and base](/development-pipeline/pull-request-branch-selection.png)

Click on 'Create pull request'

![Screenshot of the create pull request button](/development-pipeline/create-pull-request-button.png)

### Fill in the pull request template

![Screenshot of the pull request template](/development-pipeline/pull-request-template-boxes.png)

**PLEASE COMPLETE THE FOLLOWING:**

- Give your pull request a title
- Fill in the description with what you have done
- Replace [number] with the number of the issue that the pull request resolves
- Convert `[ ]`to `[x]` for any checkboxes that are relevant to the pull request

### Assign people and labels to the request

Reviewers -> Assigned people will be notified to review your changes and provide feedback (code review)
\
Assignees -> Assigned people will be responsible for overseeing the pull request and the merge process, we recommend at least assigning yourself
\
Labels -> Relevant tags for the pull request (best to pick the same ones as used in the issue being addressed)

![Screenshot of applying labels](/development-pipeline/assign-labels-pull-request.png)

### Click 'Create pull request'

![Screenshot of create pull request button](/development-pipeline/create-pull-request.png)

## Step Six: Finalise the pull request

### Code review

From here, someone else (the reviewer) should complete a code review of your changes. A separate document will be available in the future detailing how this process is expected to play out. For now, please consult [Google's code review guide](https://google.github.io/eng-practices/review/reviewer/).

### Merge branches

Once the code review has been completed, the branches can be merged. At this stage, any merge conflicts should be resolved. This is usually a simple case of choosing which lines to keep from each file.

If this proves difficult, add someone else onto the assignee list of the pull request.

### Delete the development branch

The development branch is no longer required, GitHub should prompt the user to delete the branch safely after the merge has been completed.

![Screenshot of delete branch button](/development-pipeline/delete-branch-button.png)

However, this will only delete the branch on GitHub, not the branch on your copy of the repository. In order to delete your branch locally, please complete the following:

```console
# Ensure that there are no uncommited changes to the branch
git checkout <your-branch-name>
git status

# If the last line of the output of git status is 
# 'nothing to commit, working tree clean', you can continue safely

# Delete the branch
git checkout master
git branch -d <your-branch-name>

# Verify the deletion (your branch should not appear in the output)
git branch
```

### Mark related issues as complete

Go back to the issues that were resolved with this pull request and mark them as completed/closed.

![Screenshot of closing an issue](/development-pipeline/close-issue.png)

## Step Seven

Celebrate your first successful contribution to the pipeline!
