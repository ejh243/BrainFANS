"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[5379],{4965:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>a,contentTitle:()=>s,default:()=>d,frontMatter:()=>o,metadata:()=>l,toc:()=>h});var i=n(4848),r=n(8453);const o={sidebar_class_name:"hidden"},s="Merge conflicts",l={id:"Developer-information/merge-conflicts",title:"Merge conflicts",description:"This page details how to resolve a simple merge conflict with an example exercise. This excerise can be carried out in the terminal and with your favourite text editor.",source:"@site/docs/Developer-information/merge-conflicts.md",sourceDirName:"Developer-information",slug:"/Developer-information/merge-conflicts",permalink:"/BrainFANS/Developer-information/merge-conflicts",draft:!1,unlisted:!1,tags:[],version:"current",frontMatter:{sidebar_class_name:"hidden"},sidebar:"developerSidebar",previous:{title:"Examples of good and bad code",permalink:"/BrainFANS/Developer-information/Code-review/Code-examples"}},a={},h=[{value:"Create the git repository",id:"create-the-git-repository",level:2},{value:"A functioning merge example",id:"a-functioning-merge-example",level:2},{value:"Why the merge example worked",id:"why-the-merge-example-worked",level:2},{value:"A non-functioning merge example",id:"a-non-functioning-merge-example",level:2},{value:"Resolving a merge conflict",id:"resolving-a-merge-conflict",level:2},{value:"What do these markers mean?",id:"what-do-these-markers-mean",level:3},{value:"Resolution",id:"resolution",level:3},{value:"Further reading",id:"further-reading",level:2}];function c(e){const t={a:"a",admonition:"admonition",code:"code",em:"em",h1:"h1",h2:"h2",h3:"h3",header:"header",img:"img",li:"li",p:"p",pre:"pre",ul:"ul",...(0,r.R)(),...e.components};return(0,i.jsxs)(i.Fragment,{children:[(0,i.jsx)(t.header,{children:(0,i.jsx)(t.h1,{id:"merge-conflicts",children:"Merge conflicts"})}),"\n",(0,i.jsx)(t.p,{children:"This page details how to resolve a simple merge conflict with an example exercise. This excerise can be carried out in the terminal and with your favourite text editor."}),"\n",(0,i.jsx)(t.h2,{id:"create-the-git-repository",children:"Create the git repository"}),"\n",(0,i.jsx)(t.p,{children:"To start this tutorial, go into your terminal and initialise a new repository:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:"# Change the name of this folder to be anything you want it to be\nmkdir merge-conflict-example\ncd merge-conflict-example\n\ngit init\n"})}),"\n",(0,i.jsx)(t.p,{children:"You should see the output:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",children:"Initialized empty Git repository in /path/to/merge-conflict-example/.git/\n"})}),"\n",(0,i.jsx)(t.p,{children:"We are now going to create a simple text file on the master branch"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:'echo "line 1" > file.txt\n\n# Commit this file to the repository\ngit add file.txt\ngit commit -m "initial commit"\n'})}),"\n",(0,i.jsx)(t.p,{children:"You should see the output:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",children:"[master (root-commit) 65135f9] initial commit\n 1 file changed, 1 insertion(+)\n create mode 100644 file.txt\n"})}),"\n",(0,i.jsx)(t.h2,{id:"a-functioning-merge-example",children:"A functioning merge example"}),"\n",(0,i.jsx)(t.p,{children:"Before we create a merge conflict, we need a better understanding of what is happening when you merge two branches. To do this we are going to work through a very simple merge using a feature branch."}),"\n",(0,i.jsxs)(t.p,{children:["First, we are going to create a feature branch and change ",(0,i.jsx)(t.code,{children:"file.txt"})," on this feature branch."]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:"git branch new-branch\ngit checkout new-branch\n# >> switched to branch 'new-branch'\n\necho \"line 2\" >> file.txt\n"})}),"\n",(0,i.jsxs)(t.p,{children:["We now have ",(0,i.jsx)(t.code,{children:"file.txt"})," with the content:"]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="file.txt"',children:"line 1\nline 2\n"})}),"\n",(0,i.jsx)(t.p,{children:"Riveting."}),"\n",(0,i.jsxs)(t.p,{children:["Let's commit our important changes to ",(0,i.jsx)(t.code,{children:"file.txt"})]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:'git add file.txt\ngit commit -m "add second line"\n'})}),"\n",(0,i.jsx)(t.p,{children:"Now lets merge this branch into the master branch."}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:"git checkout master\ngit merge new-branch\n"})}),"\n",(0,i.jsx)(t.p,{children:"Our merge (provided the steps were followed correctly) should have been successful with the following terminal output:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",children:"Updating ef44a5f..f57f4b0\nFast-forward\n file.txt | 1 +\n 1 file changed, 1 insertion(+)\n"})}),"\n",(0,i.jsx)(t.p,{children:"This should all be very familiar with you already. But there is one question here. How does git know which files to update and merge during this process? Seriously, think about it. Try these commands next:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:'git checkout new-branch\n\necho "line 3" >> file.txt\ngit add file.txt\ngit commit -m "add third line"\n\ngit checkout master\n\necho "This is the first line" > new-file.txt\ngit add new-file.txt\ngit commit -m "add second file"\n\ngit merge new-branch\n# This will open up your configured editor. Just close out with the default\n# merge commit message. \n'})}),"\n",(0,i.jsxs)(t.p,{children:["Once you save and close the file ",(0,i.jsx)(t.code,{children:"MERGE_MSG"}),", you will see this in the terminal:"]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",children:"Merge made by the 'recursive' strategy.\n file.txt | 1 +\n 1 file changed, 1 insertion(+)\n"})}),"\n",(0,i.jsxs)(t.p,{children:["Let's break down what happened here. We added one more line to ",(0,i.jsx)(t.code,{children:"file.txt"}),", moved back to the master branch and created ",(0,i.jsx)(t.code,{children:"new-file.txt"}),". However, after merging git knew that we wanted to update ",(0,i.jsx)(t.code,{children:"file.txt"})," but keep ",(0,i.jsx)(t.code,{children:"new-file.txt"}),". How did git know we wanted this?"]}),"\n",(0,i.jsx)(t.h2,{id:"why-the-merge-example-worked",children:"Why the merge example worked"}),"\n",(0,i.jsxs)(t.p,{children:["When me made a change to ",(0,i.jsx)(t.code,{children:"file.txt"})," on ",(0,i.jsx)(t.code,{children:"new-branch"}),", git knew 'when' this change happened. As a result, when we merged ",(0,i.jsx)(t.code,{children:"new-branch"})," into the ",(0,i.jsx)(t.code,{children:"master"})," branch, git knew that the changed file on the feature branch is the one the user wants to carry forward. The merge results in the old ",(0,i.jsx)(t.code,{children:"file.txt"})," on the master branch being effectively substituted for the one on the feature branch. The only reason this works is because git knows which file is deemed 'old' and which one is deemed 'new'. The new file (the one from the feature branch) is the one to keep. In the same way, git knew that ",(0,i.jsx)(t.code,{children:"new-file.txt"})," was deemed the 'new state' and the absence of ",(0,i.jsx)(t.code,{children:"new-file.txt"})," on ",(0,i.jsx)(t.code,{children:"new-branch"})," was deemed the 'old state'."]}),"\n",(0,i.jsxs)(t.p,{children:["Let's dive deeper into this example. When we carried out the merge, git saw that ",(0,i.jsx)(t.code,{children:"file.txt"})," on the master branch had the same commit history as the one on the feature branch. Here is the logs of each version of the file before we merged the branches"]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="master version of file.txt commit history"',children:"[user]$ git log --oneline file.txt\nf57f4b0 (HEAD -> master) add second line\nef44a5f initial commit\n\n# Note that your commit id strings (like f57f4b0) will likely be different\n# to the above. These values are just unique identifiers given to each commit.\n"})}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="new-branch version of file.txt commit history"',children:"[user]$ git log --oneline file.txt\na655f63 (HEAD -> new-branch) add third line\nf57f4b0 (master) add second line\nef44a5f initial commit\n"})}),"\n",(0,i.jsxs)(t.p,{children:["The commit histories for ",(0,i.jsx)(t.code,{children:"file.txt"})," on both branches are compatible as the 'stem' of the history is the same on both files. Both files have commits ",(0,i.jsx)(t.code,{children:"ef44a5f"})," and ",(0,i.jsx)(t.code,{children:"f57f4b0"})," and so they are deemed compatible. To merge the versions we just need to add the commit ",(0,i.jsx)(t.code,{children:"a655f63"})," onto ",(0,i.jsx)(t.code,{children:"file.txt"}),"on the master branch. This results in both versions of ",(0,i.jsx)(t.code,{children:"file.txt"})," to have the same commit log (and therefore matching contents)."]}),"\n",(0,i.jsxs)(t.p,{children:["In the same way, ",(0,i.jsx)(t.code,{children:"new-file.txt"})," on the master branch is seen as the newest version of the file. Thus, despite the fact that we are merging from ",(0,i.jsx)(t.code,{children:"new-branch"})," into ",(0,i.jsx)(t.code,{children:"master"})," (where ",(0,i.jsx)(t.code,{children:"new-file.txt"})," doesn't exist), ",(0,i.jsx)(t.code,{children:"new-file.txt"})," is kept after the merge."]}),"\n",(0,i.jsx)(t.h2,{id:"a-non-functioning-merge-example",children:"A non-functioning merge example"}),"\n",(0,i.jsxs)(t.p,{children:["What happens when git ",(0,i.jsx)(t.em,{children:"doesn't"})," know which file to keep. Earlier, we said that git knows 'when' the file was changed. The reason for the inverted commas here is to bring attention to the fact that 'when' isn't derived from a timestamp (or something similar). Instead, all git knows is the contents of the commit history for the repository. Git doesn't actually know the exact time a file was modified (and even if it did, that wouldn't help things), but instead knows the history of commits for each file."]}),"\n",(0,i.jsx)(t.p,{children:"A merge conflict arrises when the commit history for two versions of the same file no longer line up with each other. Let's create one such conflict now:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:'# Change the file on the master branch\necho "This line was entered from the master branch" >> file.txt\ngit add file.txt\ngit commit -m "master commit"\n\n# Now we change the file differently on the new-branch\ngit checkout new-branch\necho "This line was entered from new-branch" >> file.txt\ngit add file.txt\ngit commit -m "new-branch commit"\n'})}),"\n",(0,i.jsxs)(t.p,{children:["We just brewed up a merge conflict. The two versions of ",(0,i.jsx)(t.code,{children:"file.txt"})," now differ on line 4. This usually isn't a problem, look back on our last merge and ",(0,i.jsx)(t.code,{children:"file.txt"})," differed between the two branches as well (master branch didn't have \"line 3\" on line 3 yet). Why is this case different then? Well let's look at the commit logs for both versions of ",(0,i.jsx)(t.code,{children:"file.txt"}),"."]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="master version of file.txt commit history"',children:"[user]$ git checkout master\n[user]$ git log --oneline file.txt\n5495b54 (HEAD -> master) master commit\na655f63 add third line\nf57f4b0 add second line\nef44a5f initial commit\n"})}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="new-branch version of file.txt commit history"',children:"[user]$ git checkout new-branch\n[user]$ git log --oneline file.txt\n3381285 (HEAD -> new-branch) new-branch commit\na655f63 add third line\nf57f4b0 add second line\nef44a5f initial commit\n"})}),"\n",(0,i.jsxs)(t.p,{children:["Our commit histories no longer line up. On the ",(0,i.jsx)(t.code,{children:"master"})," branch, our fourth commit to ",(0,i.jsx)(t.code,{children:"file.txt"})," is ",(0,i.jsx)(t.code,{children:"5495b54"}),", whilst on ",(0,i.jsx)(t.code,{children:"new-branch"})," it is ",(0,i.jsx)(t.code,{children:"3381285"}),". Not only this, but both commits change the same line (line 3) in the file. Git no longer knows which version of ",(0,i.jsx)(t.code,{children:"file.txt"})," is 'newer' and will ask the user for help. Note that if we changed different lines in each branch, this merge conflict would not arise. Let's see the conflict in action by attempting a merge like before."]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:"git checkout master\ngit merge new-branch\n"})}),"\n",(0,i.jsx)(t.p,{children:"You should see the following error:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",children:"Auto-merging file.txt\nCONFLICT (content): Merge conflict in file.txt\nAutomatic merge failed; fix conflicts and then commit the result.\n"})}),"\n",(0,i.jsx)(t.p,{children:"Just as expected, git doesn't know which version of the file is 'newer' and is asking us for help. Usually, the user knows which bits of each file to keep."}),"\n",(0,i.jsxs)(t.ul,{children:["\n",(0,i.jsxs)(t.li,{children:["Should we keep the 4th line from the ",(0,i.jsx)(t.code,{children:"master"})," version?"]}),"\n",(0,i.jsxs)(t.li,{children:["Should we keep the 4th line from the ",(0,i.jsx)(t.code,{children:"new-branch"})," version?"]}),"\n",(0,i.jsx)(t.li,{children:"Should we keep both lines and place them sequentially after each other?"}),"\n"]}),"\n",(0,i.jsx)(t.h2,{id:"resolving-a-merge-conflict",children:"Resolving a merge conflict"}),"\n",(0,i.jsxs)(t.p,{children:["Although our merge just failed, our repository is not in the same state it was before the merge. Try running ",(0,i.jsx)(t.code,{children:"git status"})," and you will see the following output:"]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",children:'[user]$ git status\nOn branch master\nYou have unmerged paths.\n  (fix conflicts and run "git commit")\n  (use "git merge --abort" to abort the merge)\n\nUnmerged paths:\n  (use "git add <file>..." to mark resolution)\n        both modified:   file.txt\n\nno changes added to commit (use "git add" and/or "git commit -a")\n'})}),"\n",(0,i.jsxs)(t.p,{children:["Our repository is in an intermediatory state where it is expecting us to fix the unmerged paths. To resolve these paths, we need to open ",(0,i.jsx)(t.code,{children:"file.txt"})," in our favourite text editor. Upon opening ",(0,i.jsx)(t.code,{children:"file.txt"}),", you'll see the following:"]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="file.txt"',children:" line 1\n line 2\n line 3\n <<<<<<< HEAD\n This line was entered from the master branch\n =======\n This line was entered from new-branch\n >>>>>>> new-branch\n"})}),"\n",(0,i.jsxs)(t.p,{children:["This is an intermediate form of ",(0,i.jsx)(t.code,{children:"file.txt"})," that git created during the merge process. The merge conflict can be seen between the ",(0,i.jsx)(t.code,{children:"<<<<<<<"}),", ",(0,i.jsx)(t.code,{children:"======="})," and ",(0,i.jsx)(t.code,{children:">>>>>>>"})," markers. The first three lines are the same (commit wise) between the two versions of the file, the only conflicts occur at line 4."]}),"\n",(0,i.jsx)(t.h3,{id:"what-do-these-markers-mean",children:"What do these markers mean?"}),"\n",(0,i.jsxs)(t.ul,{children:["\n",(0,i.jsxs)(t.li,{children:["The first line of the conflict: ",(0,i.jsx)(t.code,{children:"<<<<<<< HEAD"}),", marks the start of the content that is seen on HEAD.","\n",(0,i.jsxs)(t.ul,{children:["\n",(0,i.jsxs)(t.li,{children:["Recall that we saw HEAD in the git log for ",(0,i.jsx)(t.code,{children:"file.txt"}),": ",(0,i.jsx)(t.code,{children:"5495b54 (HEAD -> master) master commit"})]}),"\n",(0,i.jsxs)(t.li,{children:["HEAD means the top of the commit tree for the branch we are currently on (in our case, ",(0,i.jsx)(t.code,{children:"master"}),")."]}),"\n",(0,i.jsxs)(t.li,{children:["The content below  ",(0,i.jsx)(t.code,{children:"<<<<<<< HEAD"})," is just the content of ",(0,i.jsx)(t.code,{children:"file.txt"})," from the perspective of the ",(0,i.jsx)(t.code,{children:"master"})," branch"]}),"\n"]}),"\n"]}),"\n",(0,i.jsxs)(t.li,{children:["The third line of our conflict:  ",(0,i.jsx)(t.code,{children:"======="}),", is a separator. It signifies where the end of the HEAD (",(0,i.jsx)(t.code,{children:"master"}),") section is and where the start of the incoming change (from ",(0,i.jsx)(t.code,{children:"new-branch"}),") is."]}),"\n",(0,i.jsxs)(t.li,{children:["The final line of our conflict: ",(0,i.jsx)(t.code,{children:">>>>>>> new-branch"})," tells us where the merge conflict ends","\n",(0,i.jsxs)(t.ul,{children:["\n",(0,i.jsxs)(t.li,{children:["In this example, the merge conflict is on the last line of ",(0,i.jsx)(t.code,{children:"file.txt"})," so this isn't really necessary. However, in lots of other scenarios, merge conflicts are in the middle of a file, so this line is very useful."]}),"\n",(0,i.jsxs)(t.li,{children:["The content between ",(0,i.jsx)(t.code,{children:"======="})," and ",(0,i.jsx)(t.code,{children:">>>>>>> new-branch"})," is just the content of ",(0,i.jsx)(t.code,{children:"file.txt"})," from the perspective of ",(0,i.jsx)(t.code,{children:"new-branch"}),"."]}),"\n"]}),"\n"]}),"\n"]}),"\n",(0,i.jsx)(t.h3,{id:"resolution",children:"Resolution"}),"\n",(0,i.jsxs)(t.p,{children:["To resolve the merge conflict, we just need to remove these markers and keep the lines we want to keep. So if we wanted to keep the line we added on the master branch, we would remove the markers and the line from ",(0,i.jsx)(t.code,{children:"new-branch"})," to get:"]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="file.txt"',children:"line 1\nline 2\nline 3\nThis line was entered from the master branch\n"})}),"\n",(0,i.jsxs)(t.admonition,{type:"note",children:[(0,i.jsxs)(t.p,{children:["You may have realised that in this resolution step, all you needed to do to resolve the merge conflict was remove the markers. The content of the file can be changed in any way you wish, git will not complain if the file looks entirely different after the merge conflict (the user knows best). Despite not being recommended (as it isn't very transparent), you ",(0,i.jsx)(t.em,{children:"could"})," change the other lines (unrelated to the merge conflict) here as well. For example:"]}),(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="file.txt"',children:"This is a completely different file now. We removed all of the text and\n\nreplaced it with this text instead. The merge is still technically resolved.\n"})})]}),"\n",(0,i.jsx)(t.p,{children:"Now that we have resolved our merge conflict, we can sucessfully complete the merge. Use the following commands in your terminal:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:"git add file.txt\ngit commit\n"})}),"\n",(0,i.jsxs)(t.p,{children:["This will open your text editor with the file: ",(0,i.jsx)(t.code,{children:"COMMIT_EDITMSG"}),". In this file you can see a summary of the merge conflict:"]}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="COMMIT_EDITMSG"',children:"# Conflicts:\n#\tfile.txt\n#\n# It looks like you may be committing a merge.\n# If this is not correct, please remove the file\n#\t.git/MERGE_HEAD\n# and try again.\n\n\n# Please enter the commit message for your changes. Lines starting\n# with '#' will be ignored, and an empty message aborts the commit.\n#\n# On branch master\n# All conflicts fixed but you are still merging.\n#\n"})}),"\n",(0,i.jsx)(t.p,{children:"Add a commit message to this file, then save and close it. Your commit message should be like any other commit message. However, we recommend that you add how you resolved the merge conflict. For our example, we should put something along the lines of:"}),"\n",(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-text",metastring:'title="COMMIT_EDITMSG"',children:"resolve merge conflict with file.txt by keeping the line from master and \n\nremoving the line from new-branch. This was because...\n"})}),"\n",(0,i.jsx)(t.p,{children:"That's it. The merge conflict is resolved."}),"\n",(0,i.jsx)(t.h2,{id:"further-reading",children:"Further reading"}),"\n",(0,i.jsx)(t.p,{children:"Merge conflicts can become a bit of a headache when there are mutliple files that interact with one another. It is important to consider the implications of the changes you make when resolving a merge conflict. Sometimes you will need to keep bits of each commit tree for the file to be functional with the rest of the codebase. If at all possible, work through any merge conflicts with other developers so that you can be sure that there aren't any problems with the resolution."}),"\n",(0,i.jsx)(t.p,{children:"Our example here worked purely in the command line, but sometimes you can resolve conflicts on GitHub as well (provided there aren't too many affected files). If you ever end up creating a merge conflict with your pull request, you will see this at the bottom of the pull request."}),"\n",(0,i.jsx)(t.p,{children:(0,i.jsx)(t.img,{alt:"merge conflict box on GitHub",src:n(953).A+"",width:"890",height:"114"})}),"\n",(0,i.jsxs)(t.p,{children:["Clicking on 'Resolve conflicts' will bring you to a online text editor on GitHub. From here you can follow the same steps as ",(0,i.jsx)(t.a,{href:"#resolving-a-merge-conflict",children:"above"}),"."]})]})}function d(e={}){const{wrapper:t}={...(0,r.R)(),...e.components};return t?(0,i.jsx)(t,{...e,children:(0,i.jsx)(c,{...e})}):c(e)}},953:(e,t,n)=>{n.d(t,{A:()=>i});const i=n.p+"assets/images/github-merge-conflicts-2b9eb9c3b5923b7907ce327d16335181.png"},8453:(e,t,n)=>{n.d(t,{R:()=>s,x:()=>l});var i=n(6540);const r={},o=i.createContext(r);function s(e){const t=i.useContext(o);return i.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function l(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(r):e.components||r:s(e.components),i.createElement(o.Provider,{value:t},e.children)}}}]);