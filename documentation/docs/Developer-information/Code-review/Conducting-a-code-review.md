---
sidebar_position: 1
title: "Conducting a code review"
description: "Learn how to complete a code review"
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Conducting a code review

When another contributor creates a pull request they will assign others to review their commits. This is a required step in the contribution process and new feature branches cannot be merged without an approved code review. If you have not conducted a code review (or submitted a pull request) before, this page will help get you up to speed.

:::warning[Good headspace]
If you, the reviewer, are having a rough day (or simply have a lot of other work to do). We recommend you do not conduct a code review whilst in this headspace. Submitting a pull request can be a scary prospect for developers (especially newer developers). A bad mood will be reflected in your review and is likely to discourage future contributions, a lack of contributions will lead to a lack of codebase improvements. **A late code review is better than an overly negative one.**
:::

## Starting the code review

To start your code review navigate over to the pull request that you have been assigned to. At the bottom of the pull request (after all listed commits) you will see these boxes:

![Screenshot of code review requirement](/code-reviews/add-review-button.png)

This will bring you to a page on GitHub that displays every change that has been made to each of the files affected by the pull request. Alternatively you could switch over to the commits tab to view changes by commit instead of by file.

![Screenshot of commits tab in code review](/code-reviews/switch-to-commits.png)

![Screenshot of commits list](/code-reviews/commit-list.png)

:::info[Large number of changed files]
Sometimes, the number of files that have changed for a single pull request can be very large. As a reviewer, if you feel like too many files have been changed for a single pull request, communicate this to the developer that has submitted the pull request. It can be very difficult to keep track of a large number of changes, making the associated code review difficult. Generally, a large number of files being changed signfies that the pull request is doing too many things at once. Pull requests should ideally tackle a single feature, not multiple.
:::

## Reviewing changes

From here, the reviewer should check through the changes that have been made to the code base. GitHub provides the reviewer with a side-by-side view for each of the changes made. Sometimes however, GitHub is unable to show the difference between two files (for example, a change to a .png file or .pdf *etc.*). In these cases, it is best to pull these changes to your local machine to inspect them properly. Some reviewers may prefer to pull the changes made to their local machine regardless (due to a preference with their text editor for example).

As of 2020 you also have access to GitHub codespaces, which is just vscode ran inside your web browser. However, this is a paid service if you use it for longer than 120 hours a month. Provided you don't use the service for longer than 120 hours a month, we recommend using this for code reviews for simplicity. Use the tabs below to find out how to use codespaces or view changes locally on your own machine. For the official documentation for codespaces, please consult [this page](https://docs.github.com/en/codespaces/developing-in-a-codespace/using-github-codespaces-for-pull-requests).

<Tabs>
  <TabItem value="codespaces" label="GitHub Codespaces" default>
    If you wish to use GitHub codespaces, simply go to the top right of the pull request page and click on the '\<\> Code' button. Then select 'Create codespace...'.

    ![Screenshot of codespaces button](/code-reviews/create-new-codespace.png)

    If you already have a code space (for a different pull request), you can use the + icon to create a new one specifically for this pull request.

    ![Screenshot of new codespaces button](/code-reviews/new-codespace.png)

    This will load up VScode in your browser with the majority of its features enabled (cannot use WSL, Dev Containers, SSH *etc.*). In this window you will have access to the terminal, which can be useful for conducting small unit tests. You can also comment any line of code (in the exact same way you can on GitHub). To find out more about commenting, see [the next section](#leaving-comments).
  </TabItem>
  <TabItem value="New-user" label="If you do not have BrainFANS locally">
    If you do not already have BrainFANS on your local machine, you will need to clone the development branch, use the follow git commands to achieve this:

    ```bash
    # Move to the directory you want the repository to be cloned to
    cd path/to/directory

    # Create a copy/clone of the online repository found on GitHub
    git clone -b development-branch-name https://github.com/ejh243/BrainFANS.git
    ```
  </TabItem>
  <TabItem value="Existing-user" label="If you already have BrainFANS locally">
    If you have a version of BrainFANS already on your local machine, use these commands in the terminal:

    ```bash
    # Move to the location of the BrainFANS repository
    cd path/to/BrainFANS

    # Fetches latest changes to the repository
    git fetch

    # Switches you over to the latest version of the specific development branch
    git checkout -b development-branch-name origin/development-branch-name
    ```
  </TabItem>
</Tabs>


## Leaving comments

As the reviewer works through the changes made, they can add comments to their review. These are a great way to communicate to the reviewee about specific lines of code. Comments could be questions about the code, improvements that could be made or words of affirmation (code reviews don't have to be entirely negative comments). To create a comment the reviewer can click on the blue '+' button for any line of code in a changed file.

:::tip[multiple line comment]
Sometimes you may want to make a comment for a selection of lines in a file. You can do this by clicking on the blue '+' button and dragging up or down.
:::

![Screenshot of blue plus button](/code-reviews/add-comment-button.png)

Some converse is likely to occur during this stage between the reviewer and the reviewee. These conversations are incredibly important and we recommend that both parties look at our list of [best practices](./Best-practices.md).

## What to look for in a code review

Code reviews can end up being very subjective, it can be difficult to decide whether or not a change to the codebase is actually an improvement. To help with this we have created an ordered list of what makes code better. The list starts with the most important aspects the reviewer should focus on and ends with the least important. For code examples that emphasise the below points, click on the sub-headings.

### [Pivotal: Design](./Code-examples.md#design)
Above all else, the code should be written in a way such that it will integrate well into the existing codebase. 

Is the feature even required? 

If the pull request had been accepted without a review would anything else break? 

If the any of the above points fail, then this is a good reason to reject the pull request. If scripts interact in anyway, a change to one script could impact another. If the changes and the existing repository do not mesh together well, no amount of refactoring will fix this issue. We do not want to accept pull requests that would break the existing architecture of the repository.

### [Critical: Functionality](./Code-examples.md#functionality)

The code should be functional. Alongside design, this is the most important aspect to assess in a code review. The functionality of code can be checked using test data on the reviewer's local machine (see [above](#reviewing-changes) for how to get the changes onto your local machine). If a change causes code to become non-functional (or if specific edge cases are not accounted for), the code **should not** be approved until it is fixed.

If you, as the reviewer, can identify any non-functional code and provide suggestions for how to fix it in your comments, please do this (it would be highly appreciated). Please note however that the reviewer is not *required* to provide solutions and that the reviewee is not entilted to this treatment.

### [Important: Complexitity](./Code-examples.md#complexity)

Are indivdual lines/files/functions too complex? 

If code is too complex, it is unlikely that these lines/files/functions will be touched in the future by other developers. Code that no one wants to touch is bad code and should be avoided.

### [Beneficial: Readability/Understandability/Maintainability](./Code-examples.md#readabilityunderstandabilitymaintainability)

Is the code easy to read, understand and maintain?

These usually go hand in hand. There is some subtlety between them, but for the most part they boil down to the same question:

Do you, the reviewer, feel like you would be comfortable (though not necessarily happy) working with this code in the future (ignoring domain knowledge gaps)?

If you, the reviewer, are struggling to understand some sections of code or find yourself getting lost: **Do not** skip over it. Leave a comment and ask. Leaving a comment is the first step towards resolving problems of this type.

### [Beneficial: Comments](./Code-examples.md#comments)

Are there any comments written in plain English? 

Are these comments actually necessary? 

A good rule of thumb is that comments should explain "*why, not what*". If a comment is required to explain the *what*, then this usually signifies that the code is not written well or the idea behind the code is too complex.

### [Relevant: Scalability/Expandability](./Code-examples.md#scalabilityexpandability)

Can the changed files be expanded in scope? 

Ideally, a script (or source code) can be easily repurposed if required. A tool might work for a specific dataset, but a better tool could work with other datasets with a few small changes. The script (or source code) does not need to be multi-purpose itself, but the ability to be such can be highly beneficial.

## Finalising the code review

Once the reviewer has worked through all of the changes made to the code base, they can then make their decision for the pull request:

1) Approve the pull request, then merge the changes (working through any merge conflicts first)
2) Start a conversation with the reviewee regarding the changes required to get approval (this is the most common decision)
3) Reject the pull request. We hope this doesn't happen, but sometimes this is required (perhaps the issue has already been resolved by another pull request that just did it better)

![Screenshot of finishing a review](/code-reviews/finish-review.png)

Add a summarising comment to the text box in the 'Finish your review' window and hit 'Submit review' to finish. 

Thank you for reviewing a pull request. We know the code review process can take a long time. We hope you learned something along the way and that the BrainFANS repository has improved as a result.
