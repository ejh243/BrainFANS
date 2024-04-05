---
sidebar_position: 2
title: "Step 2: Get the latest version of BrainFANS"
description: Update or clone the repository to your local machine
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Step 2: Get the latest version of BrainFANS
Most contributions will be a part of an existing development branch. For example, if the issue you are resolving relates to atac pipeline, then you should have the latest version of the `atac-develop` pipeline. The commands you will need to use will change depending on the user, please pick the option below that best describes you.


<Tabs>
  <TabItem value="New-user" label="First time cloning BrainFANS" default>
     Please use the following git commands to clone the specific development branch of the repository to a personal directory of yours:

    ```bash
    # Move to the directory you want the repository to be cloned to
    cd path/to/directory

    # Create a copy/clone of the online repository found on GitHub
    git clone -b <branch-name> https://github.com/ejh243/BrainFANS.git
    ```
  </TabItem>
  <TabItem value="Existing-user" label="You do not have the development branch locally">
    If you have already cloned the repository in the past, there is no need to clone it again. Instead, use these commands to get the latest version of the development branch you plan to work off of.

    ```bash
    # Move to the location of the BrainFANS repository
    cd path/to/BrainFANS

    # Fetches latest changes to the repository
    git fetch

    # Switches you over to the latest version of the specific development branch 
    # you plan to work off of
    git checkout -b <development-branch-name> origin/<development-branch-name>
    ```
  </TabItem>
  <TabItem value="Existing-branch" label="You already have the development branch locally">
    If you already have the development branch locally from the past, it might not be up to date with the online repository. You will need to pull the latest changes first.

    ```bash
    # Move to the location of the BrainFANS repository
    cd path/to/BrainFANS

    # Move over to the development branch
    git checkout <development-branch-name>

    # Pulls (and merges) the latest changes to the development branch
    git pull origin <development-branch-name> 
    ```
  </TabItem>
  <TabItem value="Non contributor" label="You are not a contributor">
    If you are not a BrainFANS contributor, you will not be able to submit/push
    any of your changes directly to the repository. Instead, you will need to
    create a [fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo) 
    of the repository. To create a fork, you will need to head over to the
    GitHub page for [BrainFANS](https://github.com/ejh243/BrainFANS) and click
    on the fork button at the top right of the page.

    ![fork button](/Forks/Create-fork-button.png)

    Once you have forked the repository, you should clone this to your local
    machine by typing the following into the terminal:

    ```bash
    # Move to the directory you want the repository to be cloned to
    cd path/to/directory

    # 'Username' should be substituted with your GitHub username
    git clone https://github.com/username/BrainFANS.git
    ```

    If you are working with a fork, you will need to periodically sync your
    fork with the changes made in the official BrainFANS repository (and
    subsequently use `git fetch` on your local version). You can do this by
    clicking on the 'Sync fork' button on the fork's GitHub page.

    ![Sync fork button](/Forks/Sync-fork-button.png)

    To bring these changes to your local branch, type the following into the
    terminal:

    ```bash
    cd path/to/BrainFANS
    git fetch
    ```
  </TabItem>
</Tabs>
