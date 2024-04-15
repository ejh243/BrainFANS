---
sidebar_class_name: hidden
---
import Tabs from '@theme/Tabs';
import TabItem from '@theme/TabItem';

# Step 0: Setting up a GitHub account

The steps in the contribution guides assume that the user has already created
and set up their github account and has git installed on their system. The
following page details how to set this up.

## Installing git

<Tabs>
  <TabItem value="Windows" label="Windows Users" default>
    If you are working on windows, I'm so sorry to hear that. Ideally you
    are working on a WSL client (and then you can follow the steps on linux).
    If this is not possible (or is too much of a hassle), then you'll need to
    get [Git for Windows](https://gitforwindows.org).

    Follow the instructions from the wizard and you should be done. Access
    git via git bash, powershell or the terminal.
  </TabItem>
  <TabItem value="Mac" label="Mac Users">
    You likely already have git installed if you are a mac user. Type
    `git version` into the terminal to check for this. 

    In the event that you don't have git, you can install the program via
    homebrew with:

    `brew install git`
  </TabItem>
  <TabItem value="Linux" label="Linux Users">
    You can rather easily install git on linux using the package manager for
    your specific distribution. Examples include:

    ```bash
    # Debian/Ubuntu
    sudo apt-get install git-all

    # Fedora
    sudo yum install git
    sudo dnf install git

    # Arch (btw)
    sudo pacman -S git
    ```

    Obviously there are a lot of distros, look 
    [here](https://git-scm.com/download/linux) for a comprehensive list.
  </TabItem>
</Tabs>

## Creating a GitHub account

With git, you can immediately start version control. But in order to contribute
to the BrainFANS repository, you will need a GitHub account. This is a simple
matter of going to [GitHub](https://github.com/), then clicking on "sign up 
for GitHub", then following the rest of the instructions.

Once you have signed up, please proceed to log into your account.

## Setting up your account

:::warning[This is actually important]
So far everything has been very obvious. This is the step that isn't. GitHub
no longer allows users to use their account's password on the command line
(for security reasons). Using your github account's password will not work.
:::

<Tabs>
  <TabItem value="https" label="https" default>
    For most users, it is easier to set up your Git-GitHub connection using
    https. As stated above, you cannot use your normal account password when
    using commands like `git push`. If you don't believe this, try it out
    yourself, try pushing some changes to a repository. You'll be greeted with:

    ```text
    fatal: Authentication failed for 'https://github.com/user/repo.git'
    ```

    No, you didn't type your password incorrectly, this simply just doesn't
    work anymore. Instead you'll need to generate a personal access token
    (PAT) to use specifically for these kinds of operations.

    First, go to [GitHub](https://github.com) and log into your account.

    Second, go to [this page](https://github.com/settings/tokens) and follow
    these steps:

    ![PAT step1](/git-PATS/PAT-step1.png)
    ![PAT step2](/git-PATS/PAT-step2.png)
    ![PAT step3](/git-PATS/PAT-step3.png)
    ![PAT step4](/git-PATS/PAT-step4.png)
    ![PAT step5](/git-PATS/PAT-step5.png)

    Keep a note of this access token securely. If you lose it, you'll just
    have to make another token.
    ![PAT step6](/git-PATS/PAT-step6.png)

    Note that these tokens have a limited life span, so you'll need to repeat
    this process often. There is also an option to set the PAT to last
    indefinitely. DO NOT DO THIS (it is a big security risk).
  </TabItem>
  <TabItem value="ssh" label="ssh">
    Alternatively, you can use an ssh connection in order to push your changes
    to GitHub. This is preferred by some developers as you no longer need to
    enter a password whenever you use `git push` (or similar) and it is 
    generally safer.

    In order to use an ssh connection you will need to set up an ssh key pair
    (you don't have to do this if you already have an ssh keypair that you
    want to reuse). You can accomplish this by heading over to the terminal 
    and typing:

    ```bash
    # Generates a new ssh key pair. You can replace the default name if you
    # want to organise your ssh keys, otherwise just press enter on the text
    # prompt. Note: ed25519 is a type of encryption key

    ssh-keygen -t ed25519 -C "your_email@domain.com" 
    ```

    Next go to [this page on GitHub](https://github.com/settings/keys) and
    click on "New SSH Key".

    ![New ssh key button](/SSH-token/step1.png)

    Head back over to the terminal and navigate to wherever the ssh key was
    generated to (by default this will be ~/.ssh/). Then extract the contents 
    of the public key to your clipboard. 

    ```bash title="Example of extracting public key"
    # You don't have to use cat, anything that prints the file is sufficient
    cat ~/.ssh/Your-SSH-Key.pub
    ```
    Be sure to get the public key and not the private key (the public key has 
    the file extension `.pub`). 

    Now go back to the page on GitHub and paste the public key into the box
    labelled "key". Give the key a name (the name doesn't matter) and click
    add SSH key.

    ![Add ssh key button](/SSH-token/step2.png)

    Now, your terminal will not yet recognise this key as one to use when
    connecting to GitHub via ssh. To make things easier for you, you can add
    the following commands to your `~/.bashrc` (or the rc for whatever shell
    you use) file.

    ```bash
    # .rc file
    # Replace -s with -c if using a C-type shell (like csh)
    eval "$(ssh-agent -s)" 
    ssh-add "path/to/your/private/key"
    ```

    For the first time doing this, you will need to restart your shell or
    `source` the rc file before this will work. After however, all you will
    need to do is input your key's passphrase (if you set one up) whenever you
    open a new shell.

  </TabItem>
</Tabs>

:::warning[incorrect information for ssh]
It is assumed on these pages that the developer is using a https connection
and not an ssh connection to GitHub. As such, if you plan on using an ssh 
connection, you will need to change the remote url. Please change instances of:
> https://github.com/ejh243/BrainFANS.git

to:
> git@github.com:ejh243/BrainFANS.git

If you have already cloned the repoistory. You will need to change the
origin location. You can do this with:

```bash
git remote set-url origin git@github.com:ejh243/BrainFANS.git
```

Please check that this has indeed worked with:

```bash
# Displays the remote urls
git remote -v
```
:::
