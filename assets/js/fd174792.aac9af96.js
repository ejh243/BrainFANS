"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[3017],{7846:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>h,contentTitle:()=>l,default:()=>p,frontMatter:()=>s,metadata:()=>c,toc:()=>d});var i=n(5893),o=n(1151),a=n(4866),r=n(5162);const s={sidebar_position:1,title:"Conducting a code review",description:"Learn how to complete a code review"},l="Conducting a code review",c={id:"Developer-information/Code-review/Conducting-a-code-review",title:"Conducting a code review",description:"Learn how to complete a code review",source:"@site/docs/Developer-information/Code-review/Conducting-a-code-review.md",sourceDirName:"Developer-information/Code-review",slug:"/Developer-information/Code-review/Conducting-a-code-review",permalink:"/BrainFANS/Developer-information/Code-review/Conducting-a-code-review",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:1,frontMatter:{sidebar_position:1,title:"Conducting a code review",description:"Learn how to complete a code review"},sidebar:"developerSidebar",previous:{title:"Introduction to Code Reviews",permalink:"/BrainFANS/category/introduction-to-code-reviews"},next:{title:"Code Review - Best Practices",permalink:"/BrainFANS/Developer-information/Code-review/Best-practices"}},h={},d=[{value:"Starting the code review",id:"starting-the-code-review",level:2},{value:"Reviewing changes",id:"reviewing-changes",level:2},{value:"Leaving comments",id:"leaving-comments",level:2},{value:"What to look for in a code review",id:"what-to-look-for-in-a-code-review",level:2},{value:"Pivotal: Design",id:"pivotal-design",level:3},{value:"Critical: Functionality",id:"critical-functionality",level:3},{value:"Important: Complexitity",id:"important-complexitity",level:3},{value:"Beneficial: Readability/Understandability/Maintainability",id:"beneficial-readabilityunderstandabilitymaintainability",level:3},{value:"Beneficial: Comments",id:"beneficial-comments",level:3},{value:"Relevant: Scalability/Expandability",id:"relevant-scalabilityexpandability",level:3},{value:"Finalising the code review",id:"finalising-the-code-review",level:2}];function u(e){const t={a:"a",admonition:"admonition",code:"code",em:"em",h1:"h1",h2:"h2",h3:"h3",img:"img",li:"li",ol:"ol",p:"p",pre:"pre",strong:"strong",...(0,o.a)(),...e.components};return(0,i.jsxs)(i.Fragment,{children:[(0,i.jsx)(t.h1,{id:"conducting-a-code-review",children:"Conducting a code review"}),"\n",(0,i.jsx)(t.p,{children:"When another contributor creates a pull request they will assign others to review their commits. This is a required step in the contribution process and new feature branches cannot be merged without an approved code review. If you have not conducted a code review (or submitted a pull request) before, this page will help get you up to speed."}),"\n",(0,i.jsx)(t.admonition,{title:"Good headspace",type:"warning",children:(0,i.jsxs)(t.p,{children:["If you, the reviewer, are having a rough day (or simply have a lot of other work to do). We recommend you do not conduct a code review whilst in this headspace. Submitting a pull request can be a scary prospect for developers (especially newer developers). A bad mood will be reflected in your review and is likely to discourage future contributions, a lack of contributions will lead to a lack of codebase improvements. ",(0,i.jsx)(t.strong,{children:"A late code review is better than an overly negative one."})]})}),"\n",(0,i.jsx)(t.h2,{id:"starting-the-code-review",children:"Starting the code review"}),"\n",(0,i.jsx)(t.p,{children:"To start your code review navigate over to the pull request that you have been assigned to. At the bottom of the pull request (after all listed commits) you will see these boxes:"}),"\n",(0,i.jsx)(t.p,{children:(0,i.jsx)(t.img,{alt:"Screenshot of code review requirement",src:n(1733).Z+"",width:"763",height:"200"})}),"\n",(0,i.jsx)(t.p,{children:"This will bring you to a page on GitHub that displays every change that has been made to each of the files affected by the pull request. Alternatively you could switch over to the commits tab to view changes by commit instead of by file."}),"\n",(0,i.jsx)(t.p,{children:(0,i.jsx)(t.img,{alt:"Screenshot of commits tab in code review",src:n(718).Z+"",width:"647",height:"51"})}),"\n",(0,i.jsx)(t.p,{children:(0,i.jsx)(t.img,{alt:"Screenshot of commits list",src:n(2290).Z+"",width:"403",height:"212"})}),"\n",(0,i.jsx)(t.admonition,{title:"Large number of changed files",type:"info",children:(0,i.jsx)(t.p,{children:"Sometimes, the number of files that have changed for a single pull request can be very large. As a reviewer, if you feel like too many files have been changed for a single pull request, communicate this to the developer that has submitted the pull request. It can be very difficult to keep track of a large number of changes, making the associated code review difficult. Generally, a large number of files being changed signfies that the pull request is doing too many things at once. Pull requests should ideally tackle a single feature, not multiple."})}),"\n",(0,i.jsx)(t.h2,{id:"reviewing-changes",children:"Reviewing changes"}),"\n",(0,i.jsxs)(t.p,{children:["From here, the reviewer should check through the changes that have been made to the code base. GitHub provides the reviewer with a side-by-side view for each of the changes made. Sometimes however, GitHub is unable to show the difference between two files (for example, a change to a .png file or .pdf ",(0,i.jsx)(t.em,{children:"etc."}),"). In these cases, it is best to pull these changes to your local machine to inspect them properly. Some reviewers may prefer to pull the changes made to their local machine regardless (due to a preference with their text editor for example)."]}),"\n",(0,i.jsxs)(t.p,{children:["As of 2020 you also have access to GitHub codespaces, which is just vscode ran inside your web browser. However, this is a paid service if you use it for longer than 120 hours a month. Provided you don't use the service for longer than 120 hours a month, we recommend using this for code reviews for simplicity. Use the tabs below to find out how to use codespaces or view changes locally on your own machine. For the official documentation for codespaces, please consult ",(0,i.jsx)(t.a,{href:"https://docs.github.com/en/codespaces/developing-in-a-codespace/using-github-codespaces-for-pull-requests",children:"this page"}),"."]}),"\n",(0,i.jsxs)(a.Z,{children:[(0,i.jsxs)(r.Z,{value:"codespaces",label:"GitHub Codespaces",default:!0,children:[(0,i.jsx)(t.p,{children:"If you wish to use GitHub codespaces, simply go to the top right of the pull request page and click on the '<> Code' button. Then select 'Create codespace...'."}),(0,i.jsx)(t.p,{children:(0,i.jsx)(t.img,{alt:"Screenshot of codespaces button",src:n(960).Z+"",width:"364",height:"316"})}),(0,i.jsx)(t.p,{children:"If you already have a code space (for a different pull request), you can use the + icon to create a new one specifically for this pull request."}),(0,i.jsx)(t.p,{children:(0,i.jsx)(t.img,{alt:"Screenshot of new codespaces button",src:n(3404).Z+"",width:"353",height:"53"})}),(0,i.jsxs)(t.p,{children:["This will load up VScode in your browser with the majority of its features enabled (cannot use WSL, Dev Containers, SSH ",(0,i.jsx)(t.em,{children:"etc."}),"). In this window you will have access to the terminal, which can be useful for conducting small unit tests. You can also comment any line of code (in the exact same way you can on GitHub). To find out more about commenting, see ",(0,i.jsx)(t.a,{href:"#leaving-comments",children:"the next section"}),"."]})]}),(0,i.jsxs)(r.Z,{value:"New-user",label:"If you do not have BrainFANS locally",children:[(0,i.jsx)(t.p,{children:"If you do not already have BrainFANS on your local machine, you will need to clone the development branch, use the follow git commands to achieve this:"}),(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:"# Move to the directory you want the repository to be cloned to\ncd path/to/directory\n\n# Create a copy/clone of the online repository found on GitHub\ngit clone -b development-branch-name https://github.com/ejh243/BrainFANS.git\n"})})]}),(0,i.jsxs)(r.Z,{value:"Existing-user",label:"If you already have BrainFANS locally",children:[(0,i.jsx)(t.p,{children:"If you have a version of BrainFANS already on your local machine, use these commands in the terminal:"}),(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:"# Move to the location of the BrainFANS repository\ncd path/to/BrainFANS\n\n# Fetches latest changes to the repository\ngit fetch\n\n# Switches you over to the latest version of the specific development branch\ngit checkout -b development-branch-name origin/development-branch-name\n"})}),(0,i.jsxs)(t.p,{children:["Alternatively, if you have the ",(0,i.jsx)(t.a,{href:"https://github.com/cli/cli",children:"GitHub cli tool"}),"\ninstalled you can directly checkout to the latest version of the pull\nrequest with:"]}),(0,i.jsx)(t.pre,{children:(0,i.jsx)(t.code,{className:"language-bash",children:"# Move to the location of the BrainFANS repository\ncd path/to/BrainFANS\n\n# Gets latest version of the code seen in the pull request and places you\n# on that branch.\ngh pr checkout pull-request-number\n"})}),(0,i.jsxs)(t.p,{children:["The pull request number can be found in multiple areas on GitHub (for\nexample, next to the title of the pull request). However, if you are using\nthis cli tool, you may as well use ",(0,i.jsx)(t.code,{children:"gh pr list"})," to list all pull requests.\nAlternatively, you could use ",(0,i.jsx)(t.code,{children:"gh pr status"})," to see how the pull requests\nfor the repository are coming along. This has the added bonus of showing\nyou which pull requests you are assigned as a reviewer for."]})]})]}),"\n",(0,i.jsx)(t.h2,{id:"leaving-comments",children:"Leaving comments"}),"\n",(0,i.jsx)(t.p,{children:"As the reviewer works through the changes made, they can add comments to their review. These are a great way to communicate to the reviewee about specific lines of code. Comments could be questions about the code, improvements that could be made or words of affirmation (code reviews don't have to be entirely negative comments). To create a comment the reviewer can click on the blue '+' button for any line of code in a changed file."}),"\n",(0,i.jsx)(t.admonition,{title:"multiple line comment",type:"tip",children:(0,i.jsx)(t.p,{children:"Sometimes you may want to make a comment for a selection of lines in a file. You can do this by clicking on the blue '+' button and dragging up or down."})}),"\n",(0,i.jsx)(t.p,{children:(0,i.jsx)(t.img,{alt:"Screenshot of blue plus button",src:n(5430).Z+"",width:"708",height:"254"})}),"\n",(0,i.jsxs)(t.p,{children:["Some converse is likely to occur during this stage between the reviewer and the reviewee. These conversations are incredibly important and we recommend that both parties look at our list of ",(0,i.jsx)(t.a,{href:"/BrainFANS/Developer-information/Code-review/Best-practices",children:"best practices"}),"."]}),"\n",(0,i.jsx)(t.h2,{id:"what-to-look-for-in-a-code-review",children:"What to look for in a code review"}),"\n",(0,i.jsx)(t.p,{children:"Code reviews can end up being very subjective, it can be difficult to decide whether or not a change to the codebase is actually an improvement. To help with this we have created an ordered list of what makes code better. The list starts with the most important aspects the reviewer should focus on and ends with the least important. For code examples that emphasise the below points, click on the sub-headings."}),"\n",(0,i.jsx)(t.h3,{id:"pivotal-design",children:(0,i.jsx)(t.a,{href:"/BrainFANS/Developer-information/Code-review/Code-examples#design",children:"Pivotal: Design"})}),"\n",(0,i.jsx)(t.p,{children:"Above all else, the code should be written in a way such that it will integrate well into the existing codebase."}),"\n",(0,i.jsx)(t.p,{children:"Is the feature even required?"}),"\n",(0,i.jsx)(t.p,{children:"If the pull request had been accepted without a review would anything else break?"}),"\n",(0,i.jsx)(t.p,{children:"If the any of the above points fail, then this is a good reason to reject the pull request. If scripts interact in anyway, a change to one script could impact another. If the changes and the existing repository do not mesh together well, no amount of refactoring will fix this issue. We do not want to accept pull requests that would break the existing architecture of the repository."}),"\n",(0,i.jsx)(t.h3,{id:"critical-functionality",children:(0,i.jsx)(t.a,{href:"/BrainFANS/Developer-information/Code-review/Code-examples#functionality",children:"Critical: Functionality"})}),"\n",(0,i.jsxs)(t.p,{children:["The code should be functional. Alongside design, this is the most important aspect to assess in a code review. The functionality of code can be checked using test data on the reviewer's local machine (see ",(0,i.jsx)(t.a,{href:"#reviewing-changes",children:"above"})," for how to get the changes onto your local machine). If a change causes code to become non-functional (or if specific edge cases are not accounted for), the code ",(0,i.jsx)(t.strong,{children:"should not"})," be approved until it is fixed."]}),"\n",(0,i.jsxs)(t.p,{children:["If you, as the reviewer, can identify any non-functional code and provide suggestions for how to fix it in your comments, please do this (it would be highly appreciated). Please note however that the reviewer is not ",(0,i.jsx)(t.em,{children:"required"})," to provide solutions and that the reviewee is not entilted to this treatment."]}),"\n",(0,i.jsx)(t.h3,{id:"important-complexitity",children:(0,i.jsx)(t.a,{href:"/BrainFANS/Developer-information/Code-review/Code-examples#complexity",children:"Important: Complexitity"})}),"\n",(0,i.jsx)(t.p,{children:"Are indivdual lines/files/functions too complex?"}),"\n",(0,i.jsx)(t.p,{children:"If code is too complex, it is unlikely that these lines/files/functions will be touched in the future by other developers. Code that no one wants to touch is bad code and should be avoided."}),"\n",(0,i.jsx)(t.h3,{id:"beneficial-readabilityunderstandabilitymaintainability",children:(0,i.jsx)(t.a,{href:"/BrainFANS/Developer-information/Code-review/Code-examples#readabilityunderstandabilitymaintainability",children:"Beneficial: Readability/Understandability/Maintainability"})}),"\n",(0,i.jsx)(t.p,{children:"Is the code easy to read, understand and maintain?"}),"\n",(0,i.jsx)(t.p,{children:"These usually go hand in hand. There is some subtlety between them, but for the most part they boil down to the same question:"}),"\n",(0,i.jsx)(t.p,{children:"Do you, the reviewer, feel like you would be comfortable (though not necessarily happy) working with this code in the future (ignoring domain knowledge gaps)?"}),"\n",(0,i.jsxs)(t.p,{children:["If you, the reviewer, are struggling to understand some sections of code or find yourself getting lost: ",(0,i.jsx)(t.strong,{children:"Do not"})," skip over it. Leave a comment and ask. Leaving a comment is the first step towards resolving problems of this type."]}),"\n",(0,i.jsx)(t.h3,{id:"beneficial-comments",children:(0,i.jsx)(t.a,{href:"/BrainFANS/Developer-information/Code-review/Code-examples#comments",children:"Beneficial: Comments"})}),"\n",(0,i.jsx)(t.p,{children:"Are there any comments written in plain English?"}),"\n",(0,i.jsx)(t.p,{children:"Are these comments actually necessary?"}),"\n",(0,i.jsxs)(t.p,{children:['A good rule of thumb is that comments should explain "',(0,i.jsx)(t.em,{children:"why, not what"}),'". If a comment is required to explain the ',(0,i.jsx)(t.em,{children:"what"}),", then this usually signifies that the code is not written well or the idea behind the code is too complex."]}),"\n",(0,i.jsx)(t.h3,{id:"relevant-scalabilityexpandability",children:(0,i.jsx)(t.a,{href:"/BrainFANS/Developer-information/Code-review/Code-examples#scalabilityexpandability",children:"Relevant: Scalability/Expandability"})}),"\n",(0,i.jsx)(t.p,{children:"Can the changed files be expanded in scope?"}),"\n",(0,i.jsx)(t.p,{children:"Ideally, a script (or source code) can be easily repurposed if required. A tool might work for a specific dataset, but a better tool could work with other datasets with a few small changes. The script (or source code) does not need to be multi-purpose itself, but the ability to be such can be highly beneficial."}),"\n",(0,i.jsx)(t.h2,{id:"finalising-the-code-review",children:"Finalising the code review"}),"\n",(0,i.jsx)(t.p,{children:"Once the reviewer has worked through all of the changes made to the code base, they can then make their decision for the pull request:"}),"\n",(0,i.jsxs)(t.ol,{children:["\n",(0,i.jsx)(t.li,{children:"Approve the pull request, then merge the changes (working through any merge conflicts first)"}),"\n",(0,i.jsx)(t.li,{children:"Start a conversation with the reviewee regarding the changes required to get approval (this is the most common decision)"}),"\n",(0,i.jsx)(t.li,{children:"Reject the pull request. We hope this doesn't happen, but sometimes this is required (perhaps the issue has already been resolved by another pull request that just did it better)"}),"\n"]}),"\n",(0,i.jsx)(t.p,{children:(0,i.jsx)(t.img,{alt:"Screenshot of finishing a review",src:n(6310).Z+"",width:"637",height:"535"})}),"\n",(0,i.jsx)(t.p,{children:"Add a summarising comment to the text box in the 'Finish your review' window and hit 'Submit review' to finish."}),"\n",(0,i.jsx)(t.p,{children:"Thank you for reviewing a pull request. We know the code review process can take a long time. We hope you learned something along the way and that the BrainFANS repository has improved as a result."})]})}function p(e={}){const{wrapper:t}={...(0,o.a)(),...e.components};return t?(0,i.jsx)(t,{...e,children:(0,i.jsx)(u,{...e})}):u(e)}},5162:(e,t,n)=>{n.d(t,{Z:()=>r});n(7294);var i=n(512);const o={tabItem:"tabItem_Ymn6"};var a=n(5893);function r(e){let{children:t,hidden:n,className:r}=e;return(0,a.jsx)("div",{role:"tabpanel",className:(0,i.Z)(o.tabItem,r),hidden:n,children:t})}},4866:(e,t,n)=>{n.d(t,{Z:()=>j});var i=n(7294),o=n(512),a=n(2466),r=n(6550),s=n(469),l=n(1980),c=n(7392),h=n(12);function d(e){return i.Children.toArray(e).filter((e=>"\n"!==e)).map((e=>{if(!e||(0,i.isValidElement)(e)&&function(e){const{props:t}=e;return!!t&&"object"==typeof t&&"value"in t}(e))return e;throw new Error(`Docusaurus error: Bad <Tabs> child <${"string"==typeof e.type?e.type:e.type.name}>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.`)}))?.filter(Boolean)??[]}function u(e){const{values:t,children:n}=e;return(0,i.useMemo)((()=>{const e=t??function(e){return d(e).map((e=>{let{props:{value:t,label:n,attributes:i,default:o}}=e;return{value:t,label:n,attributes:i,default:o}}))}(n);return function(e){const t=(0,c.l)(e,((e,t)=>e.value===t.value));if(t.length>0)throw new Error(`Docusaurus error: Duplicate values "${t.map((e=>e.value)).join(", ")}" found in <Tabs>. Every value needs to be unique.`)}(e),e}),[t,n])}function p(e){let{value:t,tabValues:n}=e;return n.some((e=>e.value===t))}function m(e){let{queryString:t=!1,groupId:n}=e;const o=(0,r.k6)(),a=function(e){let{queryString:t=!1,groupId:n}=e;if("string"==typeof t)return t;if(!1===t)return null;if(!0===t&&!n)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return n??null}({queryString:t,groupId:n});return[(0,l._X)(a),(0,i.useCallback)((e=>{if(!a)return;const t=new URLSearchParams(o.location.search);t.set(a,e),o.replace({...o.location,search:t.toString()})}),[a,o])]}function f(e){const{defaultValue:t,queryString:n=!1,groupId:o}=e,a=u(e),[r,l]=(0,i.useState)((()=>function(e){let{defaultValue:t,tabValues:n}=e;if(0===n.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(t){if(!p({value:t,tabValues:n}))throw new Error(`Docusaurus error: The <Tabs> has a defaultValue "${t}" but none of its children has the corresponding value. Available values are: ${n.map((e=>e.value)).join(", ")}. If you intend to show no default tab, use defaultValue={null} instead.`);return t}const i=n.find((e=>e.default))??n[0];if(!i)throw new Error("Unexpected error: 0 tabValues");return i.value}({defaultValue:t,tabValues:a}))),[c,d]=m({queryString:n,groupId:o}),[f,v]=function(e){let{groupId:t}=e;const n=function(e){return e?`docusaurus.tab.${e}`:null}(t),[o,a]=(0,h.Nk)(n);return[o,(0,i.useCallback)((e=>{n&&a.set(e)}),[n,a])]}({groupId:o}),b=(()=>{const e=c??f;return p({value:e,tabValues:a})?e:null})();(0,s.Z)((()=>{b&&l(b)}),[b]);return{selectedValue:r,selectValue:(0,i.useCallback)((e=>{if(!p({value:e,tabValues:a}))throw new Error(`Can't select invalid tab value=${e}`);l(e),d(e),v(e)}),[d,v,a]),tabValues:a}}var v=n(2389);const b={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};var g=n(5893);function w(e){let{className:t,block:n,selectedValue:i,selectValue:r,tabValues:s}=e;const l=[],{blockElementScrollPositionUntilNextRender:c}=(0,a.o5)(),h=e=>{const t=e.currentTarget,n=l.indexOf(t),o=s[n].value;o!==i&&(c(t),r(o))},d=e=>{let t=null;switch(e.key){case"Enter":h(e);break;case"ArrowRight":{const n=l.indexOf(e.currentTarget)+1;t=l[n]??l[0];break}case"ArrowLeft":{const n=l.indexOf(e.currentTarget)-1;t=l[n]??l[l.length-1];break}}t?.focus()};return(0,g.jsx)("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,o.Z)("tabs",{"tabs--block":n},t),children:s.map((e=>{let{value:t,label:n,attributes:a}=e;return(0,g.jsx)("li",{role:"tab",tabIndex:i===t?0:-1,"aria-selected":i===t,ref:e=>l.push(e),onKeyDown:d,onClick:h,...a,className:(0,o.Z)("tabs__item",b.tabItem,a?.className,{"tabs__item--active":i===t}),children:n??t},t)}))})}function y(e){let{lazy:t,children:n,selectedValue:o}=e;const a=(Array.isArray(n)?n:[n]).filter(Boolean);if(t){const e=a.find((e=>e.props.value===o));return e?(0,i.cloneElement)(e,{className:"margin-top--md"}):null}return(0,g.jsx)("div",{className:"margin-top--md",children:a.map(((e,t)=>(0,i.cloneElement)(e,{key:t,hidden:e.props.value!==o})))})}function x(e){const t=f(e);return(0,g.jsxs)("div",{className:(0,o.Z)("tabs-container",b.tabList),children:[(0,g.jsx)(w,{...e,...t}),(0,g.jsx)(y,{...e,...t})]})}function j(e){const t=(0,v.Z)();return(0,g.jsx)(x,{...e,children:d(e.children)},String(t))}},5430:(e,t,n)=>{n.d(t,{Z:()=>i});const i=n.p+"assets/images/add-comment-button-1433f9c50acbecaeaeba1bdce2c05638.png"},1733:(e,t,n)=>{n.d(t,{Z:()=>i});const i=n.p+"assets/images/add-review-button-030d55ee02833ef220e40627cc05005b.png"},2290:(e,t,n)=>{n.d(t,{Z:()=>i});const i=n.p+"assets/images/commit-list-7c75eb675c95eeaad38313cacd47f686.png"},960:(e,t,n)=>{n.d(t,{Z:()=>i});const i=n.p+"assets/images/create-new-codespace-7291fbb79cd7e81fbe4f938b7ecd6b82.png"},6310:(e,t,n)=>{n.d(t,{Z:()=>i});const i=n.p+"assets/images/finish-review-596987a4373f57248cc313d81ad4c6f2.png"},3404:(e,t,n)=>{n.d(t,{Z:()=>i});const i="data:image/png;base64,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"},718:(e,t,n)=>{n.d(t,{Z:()=>i});const i="data:image/png;base64,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"},1151:(e,t,n)=>{n.d(t,{Z:()=>s,a:()=>r});var i=n(7294);const o={},a=i.createContext(o);function r(e){const t=i.useContext(a);return i.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function s(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(o):e.components||o:r(e.components),i.createElement(a.Provider,{value:t},e.children)}}}]);