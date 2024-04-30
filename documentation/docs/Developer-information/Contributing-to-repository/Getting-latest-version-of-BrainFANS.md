---
sidebar_position: 2
title: "Step 2: Get the latest version of BrainFANS"
description: Update or clone the repository to your local machine
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Step 2: Get the latest version of BrainFANS
Most contributions will be a part of an existing development branch. For example, if the issue you are resolving relates to atac pipeline, then you should have the latest version of the `atac-develop` pipeline. The commands you will need to use will change depending on the user, please pick the option below that best describes you.

If you do not have a GitHub account already, please go to 
[this page](./Git-account-creation.md) in order to set one up. At some point 
during this guide you will no longer be able to continue, it is best you set 
one up now.

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
</Tabs>
