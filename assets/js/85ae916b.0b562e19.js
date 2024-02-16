"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[347],{6975:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>c,contentTitle:()=>l,default:()=>m,frontMatter:()=>a,metadata:()=>u,toc:()=>h});var r=n(5893),o=n(1151),i=n(4866),s=n(5162);const a={sidebar_position:6,title:"Step 6: Finalise the pull request",description:"Complete a code review, fix merge conflicts and merge"},l="Step 6: Finalise the pull request",u={id:"Developer-information/Contributing-to-repository/Finalise-pull-request",title:"Step 6: Finalise the pull request",description:"Complete a code review, fix merge conflicts and merge",source:"@site/docs/Developer-information/Contributing-to-repository/Finalise-pull-request.md",sourceDirName:"Developer-information/Contributing-to-repository",slug:"/Developer-information/Contributing-to-repository/Finalise-pull-request",permalink:"/BrainFANS/Developer-information/Contributing-to-repository/Finalise-pull-request",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:6,frontMatter:{sidebar_position:6,title:"Step 6: Finalise the pull request",description:"Complete a code review, fix merge conflicts and merge"},sidebar:"developerSidebar",previous:{title:"Step 5: Creating a pull request",permalink:"/BrainFANS/Developer-information/Contributing-to-repository/Create-pull-request"},next:{title:"Contributing to Documentation - A Quick Start Guide",permalink:"/BrainFANS/category/contributing-to-documentation---a-quick-start-guide"}},c={},h=[{value:"Code review",id:"code-review",level:2},{value:"Merge branches",id:"merge-branches",level:2},{value:"Types of merges",id:"types-of-merges",level:3},{value:"Delete your branch",id:"delete-your-branch",level:2},{value:"Mark related issues as complete",id:"mark-related-issues-as-complete",level:2},{value:"Celebrate the completion of your contribution",id:"celebrate-the-completion-of-your-contribution",level:2}];function d(e){const t={a:"a",admonition:"admonition",br:"br",code:"code",em:"em",h1:"h1",h2:"h2",h3:"h3",img:"img",p:"p",pre:"pre",strong:"strong",...(0,o.a)(),...e.components};return(0,r.jsxs)(r.Fragment,{children:[(0,r.jsx)(t.h1,{id:"step-6-finalise-the-pull-request",children:"Step 6: Finalise the pull request"}),"\n",(0,r.jsx)(t.h2,{id:"code-review",children:"Code review"}),"\n",(0,r.jsxs)(t.p,{children:["Now you have created a pull request, someone else (the reviewer that was assigned to the pull request) should complete a code review of your changes. To learn more about code reviews and some best practices, head over to ",(0,r.jsx)(t.a,{href:"/category/introduction-to-code-reviews",children:"these pages"}),"."]}),"\n",(0,r.jsx)(t.admonition,{title:"Size of pull request",type:"tip",children:(0,r.jsx)(t.p,{children:"Ideally your pull request should be in a state where a reviewer will take no more than about 1 hour to complete. In most cases this means that less than about 400 lines will have changed. If more than 400 lines have been changed, consider splitting your pull request into multiple pull requests (it is likely that the pull request is too complex). If a pull request has too many lines to review, it will likely be pushed back (and therefore less likely to be integrated into the main branch)."})}),"\n",(0,r.jsx)(t.h2,{id:"merge-branches",children:"Merge branches"}),"\n",(0,r.jsx)(t.p,{children:"Once the code review has been completed and approved, the branches can be merged. At this stage, any merge conflicts should be resolved. This is usually a simple case of choosing which lines to keep from each file."}),"\n",(0,r.jsx)(t.p,{children:"If this proves difficult, add someone else onto the assignee list of the pull request."}),"\n",(0,r.jsx)(t.h3,{id:"types-of-merges",children:"Types of merges"}),"\n",(0,r.jsxs)(t.p,{children:["Depending on the branch that is being merged, there are different types of merges that we recommend. Simply using 'merge' can lead to a lot of ",(0,r.jsx)(t.em,{children:"'commit noise'"})," in the repository over time (making it difficult to see which commits are related to one another)."]}),"\n",(0,r.jsxs)(t.p,{children:["If you plan on deleting your development branch immediately after merging, use 'sqaush and merge'. However, if you suspect that the branch will be required again in the future (which will be the case with branches like ",(0,r.jsx)(t.em,{children:"web-development"})," or ",(0,r.jsx)(t.em,{children:"atac-develop"}),"), then use a simple merge. There is also a third option available in 'rebase and merge'. To keep things simple, we currently do not recommend using rebase and merge (though you can still learn about it in the tabs below). For the majority of cases, if you are working on a feature branch, squash and merge is recommended. To learn more about these types of merges, you can use the tabs below."]}),"\n",(0,r.jsxs)(t.admonition,{title:"Squash merges on long standing branches",type:"danger",children:[(0,r.jsxs)(t.p,{children:["The most glaring problem with squash merges is when a developer continues to use their feature branch after merging. If this is you, save yourself the hassle and ",(0,r.jsx)(t.strong,{children:"stop now"}),"."]}),(0,r.jsxs)(t.p,{children:["Why?\n",(0,r.jsx)(t.br,{}),"\n","When you use squash and merge, the timeline of the git repository is pushed into a single point (a single commit). As such, when you next go to merge your changes you will be greeted with a merge conflict ",(0,r.jsx)(t.strong,{children:"for every single file you have changed"})," (which is rather unpleasant). The reason for this is that git does not know which file (the one on the base branch or the one on the comparison branch) preceeds the other. As a result, it is up to the developer to decide which lines from each file should be kept (which can be an incredibly lengthy process)."]})]}),"\n",(0,r.jsxs)(i.Z,{children:[(0,r.jsxs)(s.Z,{value:"Merge",label:"Merge",children:[(0,r.jsx)(t.p,{children:"This is the default option for merging branches with git. When you click this button, every commit you make gets effectively moved over to the base branch. This is recommended for branches that you do not want to delete."}),(0,r.jsx)(t.p,{children:"If there are a lot of small, inconsequential commits in the pull request, then this type of merge can lead to a cluttered repository. If another developer inspects commits that have been simply merged into the base branch, it can be difficult to gather all of the context behind said commit. This is because it could be surrounded by hundreds (or thousands) of other commits that aren't related to it in anyway.\nUse this merge option sparingly to avoid disorganisation in the repository."})]}),(0,r.jsxs)(s.Z,{value:"Squash",label:"Squash and merge",default:!0,children:[(0,r.jsxs)(t.p,{children:["As mentioned above, this is our recommendation for feature branches (branches you expect to delete upon merge) as it removes a lot of ",(0,r.jsx)(t.em,{children:"'commit noise'"})," in the repository."]}),(0,r.jsx)(t.p,{children:"Ideally your pull request will cover one overarching feature/enhancement. The feature/enhancement that is brought with the pull request should be used as the name of the squashed commit (and the individual commits will be automatically added to the commit description). An added bonus of squashing commits is that if another developer is searching through the repository they will be able to 'zoom and enhance' on any particular pull request (from this they should be able to gather more context behind certain file changes). This is particularly useful when an individual commit does not make much sense without an accompanying commit."}),(0,r.jsx)(t.p,{children:"To ensure that the squashed commit can be traced back easily to the pull request, ensure that the reference to the pull request (#[pull request number]) is in the title of the commit (this is done for you by default)."}),(0,r.jsx)(t.p,{children:"To do a squash merge, use this button:"}),(0,r.jsx)(t.p,{children:(0,r.jsx)(t.img,{alt:"Screenshot of squash and merge",src:n(891).Z+"",width:"499",height:"391"})})]}),(0,r.jsxs)(s.Z,{value:"Rebase",label:"Rebase and merge",children:[(0,r.jsx)(t.p,{children:"This moves the comparison branch to a new base commit, which basically rewrites the commit history of the base branch. This can be useful for maintaining a more linear commit history, which in turn can be very useful when working with integrating changes from a main branch into a long standing feature branch."}),(0,r.jsx)(t.p,{children:"The reason we don't recommend this is due to the potential of losing commits (which can lead to confusion of other developers later down the line). Also, if any merge conflicts do indeed arise due to this type of merge, it can be very difficult to resolve them, especially if there are a large number of commits involved in such merge conflicts."}),(0,r.jsx)(t.p,{children:"Only use rebase and merge if you know what you are doing and other developers are aware of the possible side effects."})]})]}),"\n",(0,r.jsx)(t.h2,{id:"delete-your-branch",children:"Delete your branch"}),"\n",(0,r.jsx)(t.p,{children:"In most cases your branch is no longer required, GitHub should prompt the user to delete the branch safely after the merge has been completed."}),"\n",(0,r.jsx)(t.p,{children:(0,r.jsx)(t.img,{alt:"Screenshot of delete branch button",src:n(8897).Z+"",width:"1794",height:"116"})}),"\n",(0,r.jsxs)(t.p,{children:["If you created a new development branch in ",(0,r.jsx)(t.a,{href:"/BrainFANS/Developer-information/Contributing-to-repository/Creating-new-branch",children:"step 3"}),", this will not apply to you. If unsure, notify an admin."]}),"\n",(0,r.jsx)(t.p,{children:"Clicking the 'delete branch' button will only delete the branch on GitHub, not the branch on your copy of the repository. If you want to delete your branch locally, please complete the following:"}),"\n",(0,r.jsx)(t.pre,{children:(0,r.jsx)(t.code,{className:"language-bash",children:"# Ensure that there are no uncommited changes to the branch\ngit checkout <your-branch-name>\ngit status\n\n# If the last line of the output of git status is \n# 'nothing to commit, working tree clean', you can continue safely\n\n# Delete the branch\ngit checkout master\ngit branch -d <your-branch-name>\n\n# Verify the deletion (your branch should not appear in the output)\ngit branch\n"})}),"\n",(0,r.jsx)(t.h2,{id:"mark-related-issues-as-complete",children:"Mark related issues as complete"}),"\n",(0,r.jsx)(t.p,{children:"Go back to the issues that were resolved with this pull request and mark them as completed/closed. This is necessary as it cleans up the issues page, making it easier to navigate."}),"\n",(0,r.jsx)(t.admonition,{title:"Automation",type:"info",children:(0,r.jsx)(t.p,{children:"GitHub allows issues to be autocompleted if the issues are linked to the pull request (located in the developments section of the pull request)."})}),"\n",(0,r.jsx)(t.p,{children:(0,r.jsx)(t.img,{alt:"Screenshot of closing an issue",src:n(309).Z+"",width:"942",height:"308"})}),"\n",(0,r.jsx)(t.h2,{id:"celebrate-the-completion-of-your-contribution",children:"Celebrate the completion of your contribution"}),"\n",(0,r.jsx)(t.p,{children:"Congratulations! You have completed a contribution. As you contribute more, this process will become second nature."}),"\n",(0,r.jsx)(t.admonition,{title:"Future contributions",type:"tip",children:(0,r.jsx)(t.p,{children:"If this was your first contribution, we recommend you continue to use this guide for your next few contributions as well."})})]})}function m(e={}){const{wrapper:t}={...(0,o.a)(),...e.components};return t?(0,r.jsx)(t,{...e,children:(0,r.jsx)(d,{...e})}):d(e)}},5162:(e,t,n)=>{n.d(t,{Z:()=>s});n(7294);var r=n(512);const o={tabItem:"tabItem_Ymn6"};var i=n(5893);function s(e){let{children:t,hidden:n,className:s}=e;return(0,i.jsx)("div",{role:"tabpanel",className:(0,r.Z)(o.tabItem,s),hidden:n,children:t})}},4866:(e,t,n)=>{n.d(t,{Z:()=>j});var r=n(7294),o=n(512),i=n(2466),s=n(6550),a=n(469),l=n(1980),u=n(7392),c=n(12);function h(e){return r.Children.toArray(e).filter((e=>"\n"!==e)).map((e=>{if(!e||(0,r.isValidElement)(e)&&function(e){const{props:t}=e;return!!t&&"object"==typeof t&&"value"in t}(e))return e;throw new Error(`Docusaurus error: Bad <Tabs> child <${"string"==typeof e.type?e.type:e.type.name}>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.`)}))?.filter(Boolean)??[]}function d(e){const{values:t,children:n}=e;return(0,r.useMemo)((()=>{const e=t??function(e){return h(e).map((e=>{let{props:{value:t,label:n,attributes:r,default:o}}=e;return{value:t,label:n,attributes:r,default:o}}))}(n);return function(e){const t=(0,u.l)(e,((e,t)=>e.value===t.value));if(t.length>0)throw new Error(`Docusaurus error: Duplicate values "${t.map((e=>e.value)).join(", ")}" found in <Tabs>. Every value needs to be unique.`)}(e),e}),[t,n])}function m(e){let{value:t,tabValues:n}=e;return n.some((e=>e.value===t))}function p(e){let{queryString:t=!1,groupId:n}=e;const o=(0,s.k6)(),i=function(e){let{queryString:t=!1,groupId:n}=e;if("string"==typeof t)return t;if(!1===t)return null;if(!0===t&&!n)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return n??null}({queryString:t,groupId:n});return[(0,l._X)(i),(0,r.useCallback)((e=>{if(!i)return;const t=new URLSearchParams(o.location.search);t.set(i,e),o.replace({...o.location,search:t.toString()})}),[i,o])]}function b(e){const{defaultValue:t,queryString:n=!1,groupId:o}=e,i=d(e),[s,l]=(0,r.useState)((()=>function(e){let{defaultValue:t,tabValues:n}=e;if(0===n.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(t){if(!m({value:t,tabValues:n}))throw new Error(`Docusaurus error: The <Tabs> has a defaultValue "${t}" but none of its children has the corresponding value. Available values are: ${n.map((e=>e.value)).join(", ")}. If you intend to show no default tab, use defaultValue={null} instead.`);return t}const r=n.find((e=>e.default))??n[0];if(!r)throw new Error("Unexpected error: 0 tabValues");return r.value}({defaultValue:t,tabValues:i}))),[u,h]=p({queryString:n,groupId:o}),[b,f]=function(e){let{groupId:t}=e;const n=function(e){return e?`docusaurus.tab.${e}`:null}(t),[o,i]=(0,c.Nk)(n);return[o,(0,r.useCallback)((e=>{n&&i.set(e)}),[n,i])]}({groupId:o}),g=(()=>{const e=u??b;return m({value:e,tabValues:i})?e:null})();(0,a.Z)((()=>{g&&l(g)}),[g]);return{selectedValue:s,selectValue:(0,r.useCallback)((e=>{if(!m({value:e,tabValues:i}))throw new Error(`Can't select invalid tab value=${e}`);l(e),h(e),f(e)}),[h,f,i]),tabValues:i}}var f=n(2389);const g={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};var y=n(5893);function v(e){let{className:t,block:n,selectedValue:r,selectValue:s,tabValues:a}=e;const l=[],{blockElementScrollPositionUntilNextRender:u}=(0,i.o5)(),c=e=>{const t=e.currentTarget,n=l.indexOf(t),o=a[n].value;o!==r&&(u(t),s(o))},h=e=>{let t=null;switch(e.key){case"Enter":c(e);break;case"ArrowRight":{const n=l.indexOf(e.currentTarget)+1;t=l[n]??l[0];break}case"ArrowLeft":{const n=l.indexOf(e.currentTarget)-1;t=l[n]??l[l.length-1];break}}t?.focus()};return(0,y.jsx)("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,o.Z)("tabs",{"tabs--block":n},t),children:a.map((e=>{let{value:t,label:n,attributes:i}=e;return(0,y.jsx)("li",{role:"tab",tabIndex:r===t?0:-1,"aria-selected":r===t,ref:e=>l.push(e),onKeyDown:h,onClick:c,...i,className:(0,o.Z)("tabs__item",g.tabItem,i?.className,{"tabs__item--active":r===t}),children:n??t},t)}))})}function w(e){let{lazy:t,children:n,selectedValue:o}=e;const i=(Array.isArray(n)?n:[n]).filter(Boolean);if(t){const e=i.find((e=>e.props.value===o));return e?(0,r.cloneElement)(e,{className:"margin-top--md"}):null}return(0,y.jsx)("div",{className:"margin-top--md",children:i.map(((e,t)=>(0,r.cloneElement)(e,{key:t,hidden:e.props.value!==o})))})}function x(e){const t=b(e);return(0,y.jsxs)("div",{className:(0,o.Z)("tabs-container",g.tabList),children:[(0,y.jsx)(v,{...e,...t}),(0,y.jsx)(w,{...e,...t})]})}function j(e){const t=(0,f.Z)();return(0,y.jsx)(x,{...e,children:h(e.children)},String(t))}},309:(e,t,n)=>{n.d(t,{Z:()=>r});const r=n.p+"assets/images/close-issue-673b9fb2d2aee1dffe3497171b93c2ce.png"},8897:(e,t,n)=>{n.d(t,{Z:()=>r});const r=n.p+"assets/images/delete-branch-button-c693a89d7456e976d561374a05c00b7e.png"},891:(e,t,n)=>{n.d(t,{Z:()=>r});const r=n.p+"assets/images/squash-commit-7cd009e44e481df4632112e89ef9294b.png"},1151:(e,t,n)=>{n.d(t,{Z:()=>a,a:()=>s});var r=n(7294);const o={},i=r.createContext(o);function s(e){const t=r.useContext(i);return r.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function a(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(o):e.components||o:s(e.components),r.createElement(i.Provider,{value:t},e.children)}}}]);