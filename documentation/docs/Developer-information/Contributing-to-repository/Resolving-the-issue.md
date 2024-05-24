---
sidebar_position: 4
title: "Step 4: Resolving the issue"
description: Making changes to the code to resolve the issue
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Step 4: Resolving the issue

## Setup

:::warning[Check this]
Before making any commits, run `git status` to check you are indeed in your 
development branch.
:::

![Screenshot of the output of git status](/development-pipeline/branch-check.png)

If this is not the branch you created in [step three](./Creating-new-branch.md), 
please run the following:

```bash
# Switch to the development branch you created in step three
git checkout your-branch-name

# Check that you are indeed on your branch (!! this is very important !!)
git status
```

:::danger[Important] 
If you are on the wrong branch when making commits you run the risk of 
creating merge conflicts that can be difficult to resolve. If you make commits 
to a development (or master) branch, you will not be able to push these 
commits at all (due to branch protection rules).
:::

If you have skipped ahead and already made commits to the master branch 
accidently, please type the following into the terminal:

```bash
# Check you are indeed on the master branch
git branch # You should see a '*' next to 'master'

# Clears any commits you have made (this does not delete your hard work)
git reset --soft origin/master

# Move over to the branch you made in step 3
git checkout your-branch-name
```

The above command will not work if your changed files will be overwritten due 
to the checkout. In such cases, git will throw the following error message:

```text
error: Your local changes to the following files would be overwritten by checkout:
    somedir/somefile.extension
Please, commit your changes or stash them before you can switch branches.
Aborting
```

In such cases, you will need to stash your changes before switching branches.
To do this, enter the following commands into the terminal.

```bash
# Stash your changes
git stash push

# Move over to your development branch without needing to commit your changes
git checkout your-branch-name

# Bring your changes back to your workspace from the stash
git stash pop
```

## Making commits

Edit the relevant files in your favourite text editor one at a time.

Once you have completed a 'small' change, you can create a commit. In order to 
create a commit, you need to add the changed file(s) to the index with 
`git add`, then commit them with `git commit`.

<Tabs>
  <TabItem value="inline" label="Writing a commit in the terminal" default>
    ```bash
    # Stage the relevant file to be committed
    git add path/to/file

    # Commit the changes
    git commit -m "Enter your commit message here"
    ```
  </TabItem>
  <TabItem value="text editor" label="Writing a commit in a text editor">
    ```bash
    # Stage the relevant file to be committed
    git add path/to/file

    # Open up your text editor for commit message
    git commit
    ```

    If you ommit the `-m` flag when using `git commit`. Git will open up a file
    called `COMMIT_MSG` which you can write your commit message into. Git will 
    open your configured text editor (if you haven't set this manually, this 
    is likely vim). To change the text editor that opens up when you type
    `git commit`, run:

    ```bash
    git config --global core.editor "your-editor-here"
    ```

    Once you are happy with your commit message, save and close the file.
  </TabItem>
</Tabs>

To make the above easier for you, get into the habit of using `git status` 
often (think of it as the git version of `ls`). Using `git status` will show 
you which files have been changed and the path to said files (which you can 
then copy and paste to the terminal).

Above we mentioned that you should commit 'small' changes. By 'small' we really 
mean [atomic](https://en.wikipedia.org/wiki/Atomic_commit). Atomic
commits are great as they can be very helpful when identifying when a bug was
first introduced into the pull request (or the master branch if the bug gets
past the pull request). An atomic commit (for our purposes) is a commit that
does one thing. If you find yourself putting 'and' or 'also' (*etc.*) into the 
commit message, you probably aren't creating an atomic commit.

You *may* feel like you want to state the intentions of your commit or perhaps
give caveats (things that might break, parts that don't quite work *etc.*) with
your commit messages. To do this, we recommend putting a separate paragraph
below the commit separated by a blank line. If you are writing the commit
in the terminal (with the `-m` flag), you can do this by:

```bash
git commit -m "Short summary of commit" -m "Longer description of commit"
```

Please don't use an extra paragraph as a way of getting around commiting large
multi-feature commits. Continue to keep your summaries short and concise and
your commits as atomic as possile.

Sometimes it makes sense to commit multiple files at the same time 
(especially when the changes are the same across each file). Don't commit loads
of files at the same time out of laziness however. Remember to keep commits
atomic wherever possible.

### Example commit messages

It can be difficult to know how to structure your commit messages. We aren't
fussy here, as long as it is clear what the commit does you can word it however
you like.

Online, you may see specific structures for commit messages. These are not
required, but might help you if you are struggling. Below we consider an
example change to a small python function and some associated commit messages.

```python
# Old code 
 def factorial(n):
     if n < 0:
         return None
     result = 1
     for i in range(1, n + 1):
         result = result * i
     return result

# New code 
 def factorial(n):
~    if not isinstance(n, int) or n < 0:
        return None
     result = 1
     for i in range(1, n + 1):
         result = result * i
     return result
```

#### Conventional commit message

The conventional commit message style is of the form:

```text title="COMMIT_MSG"
type-of-change: description

Optional longer description
```

So for our change above we would put

```text title="COMMIT_MSG"
fix: factorial now handles non-integer inputs

factorial previously threw an error if a non-integer was inputted which caused
problems in filexyz.py.
```

#### Complete the sentence commits

For some developers, it feels awkward to start the commit message. To get
around this, it can be beneficial to write your commit message by attempting
to complete the following sentence:

> If applied, this commit will...

So in our example we would write:

```text title="COMMIT_MSG"
allow the factorial function to handle non-integer inputs
```

So that the full sentence would read:
> If applied, this commit will allow the factorial function to handle 
non-integer inputs

#### Your own style

Provided that the message clearly states what you have done and possibly why
(if applicable), then the style really doesn't matter. The above examples are
just here to help developers that are struggling with this.
