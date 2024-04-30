"use strict";(self.webpackChunkdocumentation=self.webpackChunkdocumentation||[]).push([[6030],{9065:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>c,contentTitle:()=>l,default:()=>p,frontMatter:()=>s,metadata:()=>u,toc:()=>d});var o=n(5893),r=n(1151),a=n(4866),i=n(5162);const s={sidebar_position:2,title:"Step 2: Get the latest version of BrainFANS",description:"Update or clone the repository to your local machine"},l="Step 2: Get the latest version of BrainFANS",u={id:"Developer-information/Contributing-to-repository/Getting-latest-version-of-BrainFANS",title:"Step 2: Get the latest version of BrainFANS",description:"Update or clone the repository to your local machine",source:"@site/docs/Developer-information/Contributing-to-repository/Getting-latest-version-of-BrainFANS.md",sourceDirName:"Developer-information/Contributing-to-repository",slug:"/Developer-information/Contributing-to-repository/Getting-latest-version-of-BrainFANS",permalink:"/BrainFANS/Developer-information/Contributing-to-repository/Getting-latest-version-of-BrainFANS",draft:!1,unlisted:!1,tags:[],version:"current",sidebarPosition:2,frontMatter:{sidebar_position:2,title:"Step 2: Get the latest version of BrainFANS",description:"Update or clone the repository to your local machine"},sidebar:"developerSidebar",previous:{title:"Step 1: Choosing an issue",permalink:"/BrainFANS/Developer-information/Contributing-to-repository/Choosing-an-issue"},next:{title:"Step 3: Create a new branch",permalink:"/BrainFANS/Developer-information/Contributing-to-repository/Creating-new-branch"}},c={},d=[];function h(e){const t={a:"a",code:"code",h1:"h1",img:"img",p:"p",pre:"pre",...(0,r.a)(),...e.components};return(0,o.jsxs)(o.Fragment,{children:[(0,o.jsx)(t.h1,{id:"step-2-get-the-latest-version-of-brainfans",children:"Step 2: Get the latest version of BrainFANS"}),"\n",(0,o.jsxs)(t.p,{children:["Most contributions will be a part of an existing development branch. For example, if the issue you are resolving relates to atac pipeline, then you should have the latest version of the ",(0,o.jsx)(t.code,{children:"atac-develop"})," pipeline. The commands you will need to use will change depending on the user, please pick the option below that best describes you."]}),"\n",(0,o.jsxs)(t.p,{children:["If you do not have a GitHub account already, please go to\n",(0,o.jsx)(t.a,{href:"/BrainFANS/Developer-information/Contributing-to-repository/Git-account-creation",children:"this page"})," in order to set one up. At some point\nduring this guide you will no longer be able to continue, it is best you set\none up now."]}),"\n",(0,o.jsxs)(a.Z,{children:[(0,o.jsxs)(i.Z,{value:"New-user",label:"First time cloning BrainFANS",default:!0,children:[(0,o.jsx)(t.p,{children:"Please use the following git commands to clone the specific development branch of the repository to a personal directory of yours:"}),(0,o.jsx)(t.pre,{children:(0,o.jsx)(t.code,{className:"language-bash",children:"# Move to the directory you want the repository to be cloned to\ncd path/to/directory\n\n# Create a copy/clone of the online repository found on GitHub\ngit clone -b branch-name https://github.com/ejh243/BrainFANS.git\n"})})]}),(0,o.jsxs)(i.Z,{value:"Existing-user",label:"You do not have the development branch locally",children:[(0,o.jsx)(t.p,{children:"If you have already cloned the repository in the past, there is no need to clone it again. Instead, use these commands to get the latest version of the development branch you plan to work off of."}),(0,o.jsx)(t.pre,{children:(0,o.jsx)(t.code,{className:"language-bash",children:"# Move to the location of the BrainFANS repository\ncd path/to/BrainFANS\n\n# Fetches latest changes to the repository\ngit fetch\n\n# Switches you over to the latest version of the specific development branch \n# you plan to work off of\ngit checkout -b development-branch-name origin/development-branch-name\n"})})]}),(0,o.jsxs)(i.Z,{value:"Existing-branch",label:"You already have the development branch locally",children:[(0,o.jsx)(t.p,{children:"If you already have the development branch locally from the past, it might not be up to date with the online repository. You will need to pull the latest changes first."}),(0,o.jsx)(t.pre,{children:(0,o.jsx)(t.code,{className:"language-bash",children:"# Move to the location of the BrainFANS repository\ncd path/to/BrainFANS\n\n# Move over to the development branch\ngit checkout development-branch-name\n\n# Pulls (and merges) the latest changes to the development branch\ngit pull origin development-branch-name \n"})})]}),(0,o.jsxs)(i.Z,{value:"Non contributor",label:"You are not a contributor",children:[(0,o.jsxs)(t.p,{children:["If you are not a BrainFANS contributor, you will not be able to submit/push\nany of your changes directly to the repository. Instead, you will need to\ncreate a ",(0,o.jsx)(t.a,{href:"https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo",children:"fork"}),"\nof the repository. To create a fork, you will need to head over to the\nGitHub page for ",(0,o.jsx)(t.a,{href:"https://github.com/ejh243/BrainFANS",children:"BrainFANS"})," and click\non the fork button at the top right of the page."]}),(0,o.jsx)(t.p,{children:(0,o.jsx)(t.img,{alt:"fork button",src:n(6295).Z+"",width:"486",height:"53"})}),(0,o.jsx)(t.p,{children:"Once you have forked the repository, you should clone this to your local\nmachine by typing the following into the terminal:"}),(0,o.jsx)(t.pre,{children:(0,o.jsx)(t.code,{className:"language-bash",children:"# Move to the directory you want the repository to be cloned to\ncd path/to/directory\n\n# 'Username' should be substituted with your GitHub username\ngit clone https://github.com/username/BrainFANS.git\n"})}),(0,o.jsxs)(t.p,{children:["If you are working with a fork, you will need to periodically sync your\nfork with the changes made in the official BrainFANS repository (and\nsubsequently use ",(0,o.jsx)(t.code,{children:"git fetch"})," on your local version). You can do this by\nclicking on the 'Sync fork' button on the fork's GitHub page."]}),(0,o.jsx)(t.p,{children:(0,o.jsx)(t.img,{alt:"Sync fork button",src:n(8083).Z+"",width:"923",height:"59"})}),(0,o.jsx)(t.p,{children:"To bring these changes to your local branch, type the following into the\nterminal:"}),(0,o.jsx)(t.pre,{children:(0,o.jsx)(t.code,{className:"language-bash",children:"cd path/to/BrainFANS\ngit fetch\n"})})]})]})]})}function p(e={}){const{wrapper:t}={...(0,r.a)(),...e.components};return t?(0,o.jsx)(t,{...e,children:(0,o.jsx)(h,{...e})}):h(e)}},5162:(e,t,n)=>{n.d(t,{Z:()=>i});n(7294);var o=n(512);const r={tabItem:"tabItem_Ymn6"};var a=n(5893);function i(e){let{children:t,hidden:n,className:i}=e;return(0,a.jsx)("div",{role:"tabpanel",className:(0,o.Z)(r.tabItem,i),hidden:n,children:t})}},4866:(e,t,n)=>{n.d(t,{Z:()=>S});var o=n(7294),r=n(512),a=n(2466),i=n(6550),s=n(469),l=n(1980),u=n(7392),c=n(12);function d(e){return o.Children.toArray(e).filter((e=>"\n"!==e)).map((e=>{if(!e||(0,o.isValidElement)(e)&&function(e){const{props:t}=e;return!!t&&"object"==typeof t&&"value"in t}(e))return e;throw new Error(`Docusaurus error: Bad <Tabs> child <${"string"==typeof e.type?e.type:e.type.name}>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.`)}))?.filter(Boolean)??[]}function h(e){const{values:t,children:n}=e;return(0,o.useMemo)((()=>{const e=t??function(e){return d(e).map((e=>{let{props:{value:t,label:n,attributes:o,default:r}}=e;return{value:t,label:n,attributes:o,default:r}}))}(n);return function(e){const t=(0,u.l)(e,((e,t)=>e.value===t.value));if(t.length>0)throw new Error(`Docusaurus error: Duplicate values "${t.map((e=>e.value)).join(", ")}" found in <Tabs>. Every value needs to be unique.`)}(e),e}),[t,n])}function p(e){let{value:t,tabValues:n}=e;return n.some((e=>e.value===t))}function f(e){let{queryString:t=!1,groupId:n}=e;const r=(0,i.k6)(),a=function(e){let{queryString:t=!1,groupId:n}=e;if("string"==typeof t)return t;if(!1===t)return null;if(!0===t&&!n)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return n??null}({queryString:t,groupId:n});return[(0,l._X)(a),(0,o.useCallback)((e=>{if(!a)return;const t=new URLSearchParams(r.location.search);t.set(a,e),r.replace({...r.location,search:t.toString()})}),[a,r])]}function b(e){const{defaultValue:t,queryString:n=!1,groupId:r}=e,a=h(e),[i,l]=(0,o.useState)((()=>function(e){let{defaultValue:t,tabValues:n}=e;if(0===n.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(t){if(!p({value:t,tabValues:n}))throw new Error(`Docusaurus error: The <Tabs> has a defaultValue "${t}" but none of its children has the corresponding value. Available values are: ${n.map((e=>e.value)).join(", ")}. If you intend to show no default tab, use defaultValue={null} instead.`);return t}const o=n.find((e=>e.default))??n[0];if(!o)throw new Error("Unexpected error: 0 tabValues");return o.value}({defaultValue:t,tabValues:a}))),[u,d]=f({queryString:n,groupId:r}),[b,m]=function(e){let{groupId:t}=e;const n=function(e){return e?`docusaurus.tab.${e}`:null}(t),[r,a]=(0,c.Nk)(n);return[r,(0,o.useCallback)((e=>{n&&a.set(e)}),[n,a])]}({groupId:r}),g=(()=>{const e=u??b;return p({value:e,tabValues:a})?e:null})();(0,s.Z)((()=>{g&&l(g)}),[g]);return{selectedValue:i,selectValue:(0,o.useCallback)((e=>{if(!p({value:e,tabValues:a}))throw new Error(`Can't select invalid tab value=${e}`);l(e),d(e),m(e)}),[d,m,a]),tabValues:a}}var m=n(2389);const g={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};var v=n(5893);function y(e){let{className:t,block:n,selectedValue:o,selectValue:i,tabValues:s}=e;const l=[],{blockElementScrollPositionUntilNextRender:u}=(0,a.o5)(),c=e=>{const t=e.currentTarget,n=l.indexOf(t),r=s[n].value;r!==o&&(u(t),i(r))},d=e=>{let t=null;switch(e.key){case"Enter":c(e);break;case"ArrowRight":{const n=l.indexOf(e.currentTarget)+1;t=l[n]??l[0];break}case"ArrowLeft":{const n=l.indexOf(e.currentTarget)-1;t=l[n]??l[l.length-1];break}}t?.focus()};return(0,v.jsx)("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,r.Z)("tabs",{"tabs--block":n},t),children:s.map((e=>{let{value:t,label:n,attributes:a}=e;return(0,v.jsx)("li",{role:"tab",tabIndex:o===t?0:-1,"aria-selected":o===t,ref:e=>l.push(e),onKeyDown:d,onClick:c,...a,className:(0,r.Z)("tabs__item",g.tabItem,a?.className,{"tabs__item--active":o===t}),children:n??t},t)}))})}function x(e){let{lazy:t,children:n,selectedValue:r}=e;const a=(Array.isArray(n)?n:[n]).filter(Boolean);if(t){const e=a.find((e=>e.props.value===r));return e?(0,o.cloneElement)(e,{className:"margin-top--md"}):null}return(0,v.jsx)("div",{className:"margin-top--md",children:a.map(((e,t)=>(0,o.cloneElement)(e,{key:t,hidden:e.props.value!==r})))})}function E(e){const t=b(e);return(0,v.jsxs)("div",{className:(0,r.Z)("tabs-container",g.tabList),children:[(0,v.jsx)(y,{...e,...t}),(0,v.jsx)(x,{...e,...t})]})}function S(e){const t=(0,m.Z)();return(0,v.jsx)(E,{...e,children:d(e.children)},String(t))}},6295:(e,t,n)=>{n.d(t,{Z:()=>o});const o="data:image/png;base64,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"},8083:(e,t,n)=>{n.d(t,{Z:()=>o});const o="data:image/png;base64,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"},1151:(e,t,n)=>{n.d(t,{Z:()=>s,a:()=>i});var o=n(7294);const r={},a=o.createContext(r);function i(e){const t=o.useContext(a);return o.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function s(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(r):e.components||r:i(e.components),o.createElement(a.Provider,{value:t},e.children)}}}]);