---
sidebar_position: 2
title: "Step 2: Get the latest version of BrainFANS"
description: Update or clone the repository to your local machine
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Step 2: Get the latest version of BrainFANS
Whenever you are working with the BrainFANS repository, please ensure that you
have the latest version of the repository. The commands you will need to use 
will change depending on the user, please pick the option below that best 
describes you.


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
    ```
  </TabItem>
</Tabs>
