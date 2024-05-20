---
sidebar_class_name: hidden
---

# Merge conflicts

This page details how to resolve a simple merge conflict with an example exercise. This excerise can be carried out in the terminal and with your favourite text editor. 

## Create the git repository

To start this tutorial, go into your terminal and initialise a new repository:

```bash
# Change the name of this folder to be anything you want it to be
mkdir merge-conflict-example
cd merge-conflict-example

git init
```

You should see the output:

```text
Initialized empty Git repository in /path/to/merge-conflict-example/.git/
```

We are now going to create a simple text file on the master branch

```bash
echo "line 1" > file.txt

# Commit this file to the repository
git add file.txt
git commit -m "initial commit"
```

You should see the output:

```text
[master (root-commit) 65135f9] initial commit
 1 file changed, 1 insertion(+)
 create mode 100644 file.txt
```

## A functioning merge example

Before we create a merge conflict, we need a better understanding of what is happening when you merge two branches. To do this we are going to work through a very simple merge using a feature branch.

First, we are going to create a feature branch and change `file.txt` on this feature branch.

```bash
git branch new-branch
git checkout new-branch
# >> switched to branch 'new-branch'

echo "line 2" >> file.txt
```

We now have `file.txt` with the content:

```text title="file.txt"
line 1
line 2
```

Riveting.

Let's commit our important changes to `file.txt`

```bash
git add file.txt
git commit -m "add second line"
```

Now lets merge this branch into the master branch.

```bash
git checkout master
git merge new-branch
```

Our merge (provided the steps were followed correctly) should have been successful with the following terminal output:

```text
Updating ef44a5f..f57f4b0
Fast-forward
 file.txt | 1 +
 1 file changed, 1 insertion(+)
```

This should all be very familiar with you already. But there is one question here. How does git know which files to update and merge during this process? Seriously, think about it. Try these commands next:

```bash
git checkout new-branch

echo "line 3" >> file.txt
git add file.txt
git commit -m "add third line"

git checkout master

echo "This is the first line" > new-file.txt
git add new-file.txt
git commit -m "add second file"

git merge new-branch
# This will open up your configured editor. Just close out with the default
# merge commit message. 
```

Once you save and close the file `MERGE_MSG`, you will see this in the terminal:

```text
Merge made by the 'recursive' strategy.
 file.txt | 1 +
 1 file changed, 1 insertion(+)
```

Let's break down what happened here. We added one more line to `file.txt`, moved back to the master branch and created `new-file.txt`. However, after merging git knew that we wanted to update `file.txt` but keep `new-file.txt`. How did git know we wanted this? 


## Why the merge example worked

When me made a change to `file.txt` on `new-branch`, git knew 'when' this change happened. As a result, when we merged `new-branch` into the `master` branch, git knew that the changed file on the feature branch is the one the user wants to carry forward. The merge results in the old `file.txt` on the master branch being effectively substituted for the one on the feature branch. The only reason this works is because git knows which file is deemed 'old' and which one is deemed 'new'. The new file (the one from the feature branch) is the one to keep. In the same way, git knew that `new-file.txt` was deemed the 'new state' and the absence of `new-file.txt` on `new-branch` was deemed the 'old state'.

Let's dive deeper into this example. When we carried out the merge, git saw that `file.txt` on the master branch had the same commit history as the one on the feature branch. Here is the logs of each version of the file before we merged the branches

```text title="master version of file.txt commit history"
[user]$ git log --oneline file.txt
f57f4b0 (HEAD -> master) add second line
ef44a5f initial commit

# Note that your commit id strings (like f57f4b0) will likely be different
# to the above. These values are just unique identifiers given to each commit.
```

```text title="new-branch version of file.txt commit history"
[user]$ git log --oneline file.txt
a655f63 (HEAD -> new-branch) add third line
f57f4b0 (master) add second line
ef44a5f initial commit
```

The commit histories for `file.txt` on both branches are compatible as the 'stem' of the history is the same on both files. Both files have commits `ef44a5f` and `f57f4b0` and so they are deemed compatible. To merge the versions we just need to add the commit `a655f63` onto `file.txt`on the master branch. This results in both versions of `file.txt` to have the same commit log (and therefore matching contents).

In the same way, `new-file.txt` on the master branch is seen as the newest version of the file. Thus, despite the fact that we are merging from `new-branch` into `master` (where `new-file.txt` doesn't exist), `new-file.txt` is kept after the merge. 

## A non-functioning merge example

What happens when git *doesn't* know which file to keep. Earlier, we said that git knows 'when' the file was changed. The reason for the inverted commas here is to bring attention to the fact that 'when' isn't derived from a timestamp (or something similar). Instead, all git knows is the contents of the commit history for the repository. Git doesn't actually know the exact time a file was modified (and even if it did, that wouldn't help things), but instead knows the history of commits for each file. 

A merge conflict arrises when the commit history for two versions of the same file no longer line up with each other. Let's create one such conflict now:

```bash
# Change the file on the master branch
echo "This line was entered from the master branch" >> file.txt
git add file.txt
git commit -m "master commit"

# Now we change the file differently on the new-branch
git checkout new-branch
echo "This line was entered from new-branch" >> file.txt
git add file.txt
git commit -m "new-branch commit"
```

We just brewed up a merge conflict. The two versions of `file.txt` now differ on line 4. This usually isn't a problem, look back on our last merge and `file.txt` differed between the two branches as well (master branch didn't have "line 3" on line 3 yet). Why is this case different then? Well let's look at the commit logs for both versions of `file.txt`.

```text title="master version of file.txt commit history"
[user]$ git checkout master
[user]$ git log --oneline file.txt
5495b54 (HEAD -> master) master commit
a655f63 add third line
f57f4b0 add second line
ef44a5f initial commit
```

```text title="new-branch version of file.txt commit history"
[user]$ git checkout new-branch
[user]$ git log --oneline file.txt
3381285 (HEAD -> new-branch) new-branch commit
a655f63 add third line
f57f4b0 add second line
ef44a5f initial commit
```

Our commit histories no longer line up. On the `master` branch, our fourth commit to `file.txt` is `5495b54`, whilst on `new-branch` it is `3381285`. Not only this, but both commits change the same line (line 3) in the file. Git no longer knows which version of `file.txt` is 'newer' and will ask the user for help. Note that if we changed different lines in each branch, this merge conflict would not arise. Let's see the conflict in action by attempting a merge like before.

```bash
git checkout master
git merge new-branch
```

You should see the following error:

```text
Auto-merging file.txt
CONFLICT (content): Merge conflict in file.txt
Automatic merge failed; fix conflicts and then commit the result.
```

Just as expected, git doesn't know which version of the file is 'newer' and is asking us for help. Usually, the user knows which bits of each file to keep.

- Should we keep the 4th line from the `master` version? 
- Should we keep the 4th line from the `new-branch` version?
- Should we keep both lines and place them sequentially after each other?

## Resolving a merge conflict

Although our merge just failed, our repository is not in the same state it was before the merge. Try running `git status` and you will see the following output:

```text
[user]$ git status
On branch master
You have unmerged paths.
  (fix conflicts and run "git commit")
  (use "git merge --abort" to abort the merge)

Unmerged paths:
  (use "git add <file>..." to mark resolution)
        both modified:   file.txt

no changes added to commit (use "git add" and/or "git commit -a")
```

Our repository is in an intermediatory state where it is expecting us to fix the unmerged paths. To resolve these paths, we need to open `file.txt` in our favourite text editor. Upon opening `file.txt`, you'll see the following:

```text title="file.txt"
 line 1
 line 2
 line 3
 <<<<<<< HEAD
 This line was entered from the master branch
 =======
 This line was entered from new-branch
 >>>>>>> new-branch
```

This is an intermediate form of `file.txt` that git created during the merge process. The merge conflict can be seen between the `<<<<<<<`, `=======` and `>>>>>>>` markers. The first three lines are the same (commit wise) between the two versions of the file, the only conflicts occur at line 4. 

### What do these markers mean? 

- The first line of the conflict: `<<<<<<< HEAD`, marks the start of the content that is seen on HEAD. 
  - Recall that we saw HEAD in the git log for `file.txt`: `5495b54 (HEAD -> master) master commit`
  - HEAD means the top of the commit tree for the branch we are currently on (in our case, `master`).
  - The content below  `<<<<<<< HEAD` is just the content of `file.txt` from the perspective of the `master` branch
- The third line of our conflict:  `=======`, is a separator. It signifies where the end of the HEAD (`master`) section is and where the start of the incoming change (from `new-branch`) is.
- The final line of our conflict: `>>>>>>> new-branch` tells us where the merge conflict ends 
  - In this example, the merge conflict is on the last line of `file.txt` so this isn't really necessary. However, in lots of other scenarios, merge conflicts are in the middle of a file, so this line is very useful.
  - The content between `=======` and `>>>>>>> new-branch` is just the content of `file.txt` from the perspective of `new-branch`.

### Resolution

To resolve the merge conflict, we just need to remove these markers and keep the lines we want to keep. So if we wanted to keep the line we added on the master branch, we would remove the markers and the line from `new-branch` to get:

```text title="file.txt"
line 1
line 2
line 3
This line was entered from the master branch
```

:::note
You may have realised that in this resolution step, all you needed to do to resolve the merge conflict was remove the markers. The content of the file can be changed in any way you wish, git will not complain if the file looks entirely different after the merge conflict (the user knows best). Despite not being recommended (as it isn't very transparent), you *could* change the other lines (unrelated to the merge conflict) here as well. For example:

```text title="file.txt"
This is a completely different file now. We removed all of the text and

replaced it with this text instead. The merge is still technically resolved.
```
:::

Now that we have resolved our merge conflict, we can sucessfully complete the merge. Use the following commands in your terminal:

```bash
git add file.txt
git commit
```

This will open your text editor with the file: `COMMIT_EDITMSG`. In this file you can see a summary of the merge conflict:

```text title="COMMIT_EDITMSG"
# Conflicts:
#	file.txt
#
# It looks like you may be committing a merge.
# If this is not correct, please remove the file
#	.git/MERGE_HEAD
# and try again.


# Please enter the commit message for your changes. Lines starting
# with '#' will be ignored, and an empty message aborts the commit.
#
# On branch master
# All conflicts fixed but you are still merging.
#
```

Add a commit message to this file, then save and close it. Your commit message should be like any other commit message. However, we recommend that you add how you resolved the merge conflict. For our example, we should put something along the lines of:

```text title="COMMIT_EDITMSG"
resolve merge conflict with file.txt by keeping the line from master and 

removing the line from new-branch. This was because...
```

That's it. The merge conflict is resolved. 

## Further reading

Merge conflicts can become a bit of a headache when there are mutliple files that interact with one another. It is important to consider the implications of the changes you make when resolving a merge conflict. Sometimes you will need to keep bits of each commit tree for the file to be functional with the rest of the codebase. If at all possible, work through any merge conflicts with other developers so that you can be sure that there aren't any problems with the resolution.

Our example here worked purely in the command line, but sometimes you can resolve conflicts on GitHub as well (provided there aren't too many affected files). If you ever end up creating a merge conflict with your pull request, you will see this at the bottom of the pull request.

![merge conflict box on GitHub](/merge-conflicts/github-merge-conflicts.png)

Clicking on 'Resolve conflicts' will bring you to a online text editor on GitHub. From here you can follow the same steps as [above](#resolving-a-merge-conflict).
