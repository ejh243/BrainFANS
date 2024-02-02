---
sidebar_position: 4
title: "Docusaurus features"
description: "Find out about the additional features that come with docusaurus"
---

# Docusaurus features

Thanks to mdx and react (javascript framework), the pages on this wiki are more customisable than a regular markdown document. This page provides some of the features that can be added to your pages and links to additional official documentation.


## Admonitions

Considering you are on this page, you have likely seen these special blocks already. Docusaurus enables you to create five different kinds of admonitions: note, info, tip, warning and danger. You can create these admonitions on your pages by using the following syntax:

```markdown
:::info[name of info box]
This is an info box that has **markdown** _capabilites_
:::
```

Which has output:

:::info[name of info box]
This is an info box that has **markdown** _capabilites_
:::

## Mathematics

Docusaurus allows you to use mathematics in a LaTeX format with KaTeX. You can implement mathematical notation into your pages using the $ sign:

```markdown
$E=mc^2$
```

Which will return: $E=mc^2$

### Mathematical blocks

You can use a double $ sign to make mathematical expression have their own box:

```markdown
$$
E=mc^2
$$
```

Which would return:

$$
E=mc^2
$$
