---
sidebar_position: 5
title: Running the website locally 
description: Learn how to build and serve the website locally
---

import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Running the website locally

If you are making more significant changes to the website (or just want to know how your changes will look on these pages), we recommend that you build and serve the website off of your local machine. The benefits of doing this are two fold. First, you can see how your changes will actually look on these pages. Second, you can check for any build errors. Build errors can be quite common when working on building new pages, linking images or cross-linking to other pages. These are very difficult to spot normally and if such breaking changes are merged into the master branch, the website will not deploy until these errors are fixed. 

## Prerequisites

:::tip[ISCA users]
If you are an ISCA user, you do not need to install Node.js. A module is already available to use: `nodejs/18.12.1-GCCcore-12.2.0`. Load this module and ignore the installation process below. Make sure to be in a **compute node** if you are using ISCA for the next steps.
:::

In order to run the website locally, you will need Node.js (version 18.0 or greater). You can get the latest version of node in a variety of file formats via the [official website](https://nodejs.org/en/download/). After downloading Node.js, you will need to install it. If you need any assistance in doing this on linux-based systems (or if you plan on installing via the source code) follow [this guide](https://medium.com/@tgmarinho/how-to-install-node-js-via-binary-archive-on-linux-ab9bbe1dd0c2).

After installing Node.js, check that the version is indeed 18.0 or greater by running this command in the terminal:

```bash
# Checks the version of node that is installed into a directory on PATH
node -v
```

If this command returns `node: command not found` or a version of node lower than 18.0, you'll need to restart and troubleshoot, the following steps will not work without this binary.

## Building and serving

There are two steps involved with running the site locally: building, where markdown and config files are converted to html files; and serving, where the website is ran on a local server for you to access within a web browser. To complete these steps, follow the below commands in your terminal:

```bash
# Change to the 'documentation' directory in BrainFANS
cd path/to/local/BrainFANS
cd documentation

# Perform a clean install of docusaurus and its dependencies
# Ignore the warning on Array#sort(), this doesn't cause any problems
npm ci

# Build the website
# This will be slow the first time, but faster with subsequent runs
npm run build

# Run the server
npm run serve
```

:::info[build warnings]
During the `npm run build` stage, any warnings (or errors) will be written to the terminal. Take your time reading these errors/warnings. In the majority of cases, the error messages are very informative and will point you in the right direction in order to fix them.

<centre>
![Screenshot of build errors](/SSH-tunnelling/build-errors.png)
</centre>
:::

## Viewing the website

Provided `npm run build` did not return any errors you will now be able to view the updated website via the server created by `npm run serve`. If you have ran these commands on your local machine, this website should open up automatically in your default web browser. If it doesn't then go to the link that is displayed in the terminal (which is likely: http://localhost:3000/BrainFANS/).

![Screenshot of npm run serve output](/SSH-tunnelling/npm-run-serve-output.png)

## Viewing the website over ssh connection

If you are running the website on an ssh connection (ISCA users), you will not be able to view the website immediately like you would if running locally. One more step is required in ssh tunnelling. This process creates a secret pathway between your machine and the server that is running the website, allowing you to connect to the web address safely and securely.

If you do not know how to set up an SSH tunnel, click the relevant tab below to learn how.

<Tabs>
  <TabItem value="Mac/Linux" label="Mac/Linux" default>
    You can create an ssh tunnel to the server via the following command:

    ```bash
    ssh -L 3000:remote-server:3000 Username@SSH-server
    ```

    The above assumes the forward port is 3000, change this if the output of `npm run serve` says otherwise. For the other fields consult the below list:

    * Remote server -> This is the name of the node you ran the server on (for ISCA users: look at your command history, you'll see 'xxxyyy@mrc-compzzz' or 'xxxyyy@compzzz', the text after the '@' is the remote server)
    * Username -> This is the username you use to log into the SSH server
    * SSH server -> This is the name of the SSH server that you are logged into
  
    You should now be able to go to connect to the server: http://localhost:3000/BrainFANS/ 
    (Change the port number here if yours is not 3000)
  </TabItem>
  <TabItem value="MobaXterm" label="MobaXterm">
    At the top of your window you should see a 'Tunneling' (American English spelling) button. Click this button.

    ![Screenshot of MobaXterm navbar](/SSH-tunnelling/MobaXterm-navbar.png)

    A window will now pop up, click 'New SSH tunnel' at the bottom of this window.

    ![Screenshot of MobaXterm new SSH tunnel button](/SSH-tunnelling/MobaXterm-new-tunnel-button.png)

    You'll now need to fill in each of the boxes in the new window. The list of terminology below should help with this:

    * Forwarded port -> This is likely 3000. It is the number seen after 'localhost:' displayed in the output of `npm run serve`
    * SSH server -> This is the name of the SSH server that you are logged into
    * SSH login -> This is the user name you use to log into the SSH server
    * SSH port -> The port number of the SSH server (default port number is 22)
    * Remote server -> This is the name of the node you ran the server on (for ISCA users: look at your command history, you'll see 'xxxyyy@mrc-compzzz' or 'xxxyyy@compzzz', paste the text after the '@' in this box)
    * Remote port -> This will be the same as the forwarded port (likely 3000)

    Click save, name the tunnel and hit start.

    ![Screenshot of tunnel](/SSH-tunnelling/MobaXterm-tunnel.png)

    You will be prompted for your password (the one you use to connect to the SSH server), enter this and hit OK.

    ![Screenshot of password box](/SSH-tunnelling/MobaXterm-password-box.png)

    You should now be able to go to connect to the server: http://localhost:3000/BrainFANS/ 
    (Change the port number here if yours is not 3000)
  </TabItem>
</Tabs>
