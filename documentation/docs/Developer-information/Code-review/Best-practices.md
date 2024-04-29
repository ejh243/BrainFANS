---
sidebar_position: 2
title: "Code Review - Best Practices"
description: "A list of best practices for code reviewers/reviewees"
---

# Code Review - Best Practices

Code reviews are challenging and can be very scary for both the reviewer and the reviewee. This page outlines some code review principles that we consider to be best practices. We encourage all developers to read through this document in order to nuture a good working culture.

## Code reviews are **mentoring** sessions

It is all too common that more junior developers feel as though they are not 'good enough' to conduct a code review. **This is nonsense**. Code reviews are inherently a mentoring session, not a lesson. It is highly likely that the reviewer and the reviewee will learn from a code review. Perhaps the reviewer, despite being more senior, notices a neat trick that could be useful in their own projects. Perhaps the reviewer, despite being more junior, notices a small flaw in the overall design of a feature.

The important thing to note here is that a code review goes both ways. All contributors are good programmers, regardless of their experience, their opinions and points should be respected. If you ever feel like you are not qualified to review someone else's code, ask questions. Do not simply skip over sections you don't understand.

## Review the **code** not the **person**

The aim of a code review is to assess the changes made. Sometimes, mistakes are made or code is poorly written. The reviewer **should not** critque the developer that wrote the poorly written code. Reviewing code in this manner does not help anyone. In most cases, critiquing the person results in the developer being less likely to make contributions in the future. We do not want this, everyone should feel motivated towards contributing to the repository. The hard line is: "it is the code that should be critiqued, not the person".

To illustrate this point, consider the below code blocks (comments are reviewer comments):

```python title="A bad review"
# YOU did not account for negative integer inputs with YOUR function
# YOU should have implemented error catching into YOUR function to avoid this
def factorial(n): 
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```


```python title="A good review"
# Currently, this function will return the incorrect result if the user inputs
# a negative number, say factorial(-1). To fix this, the function should
# have error catching implemented that warns the user if a negative value is
# inputted
def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```

## There is no **perfect** code

The full expression is:

> There is no such thing as "perfect" code - there is only better code -- Google

In general, if the pull request has changes that improve the existing code, you should accept it. Do not constantly push back a pull request just because of a small thing you do not like. All this ends up doing is discouraging developers from contributing in the future. For a definition of "better", we suggest the reader consults our [list of what to look for in a code review](./Conducting-a-code-review.md#what-to-look-for-in-a-code-review). 

The exception to this rule is when the pull request adds uneeded functionality. For example, someone could be adding a script that implements PCA (from scratch) to some pipeline. This might work well with the current codebase and even be functional. But there is no point in reinventing the wheel, other tools do this already (and likely do it better). Such a change violates our views on [Design](./Conducting-a-code-review.md#pivotal-design) and should be rejected.

## Justify your comments

As a reviewer, it is not very helpful to the reviewee if you give unjustified comments. Be a little more specific with your comments and explain *why* the line (or section) is a problem. If you are unsure of why some line is incorrect, try starting the comment with "I *suggest* changing this because..." or "I *think* this is wrong because...". This shows that you think there is a problem with the code without definitively stating such. 

To illustrate, here are some examples of good and bad reviews:

```python title="A bad review"
# This function doesn't work for me 
def factorial(n): 
    for i in range(1, n + 1):
        result = result * i 
    return result
```


```python title="A good review"
# I think this doesn't work because 'result' has not been defined yet.
# I suggest that the code should initialise 'result' to some number
# before this loop. I think this is because on the initial loop (i=1)
# 'result' doesn't have a value and so python doesn't know how to interpret
# [empty] * 1.
def factorial(n):
    for i in range(1, n + 1):
        result = result * i 
    return result
```


## On the matter of style

Sometimes a developer writes code in a different style to the reviewer, which can cause arguments about the *best* way to write your code. Maybe one developer enjoys the power of list comprehension in python, whilst the other prefers to use loops. Provided the style doesn't severly impact the consistency across the codebase, it is not a reasonable thing to comment on (or reject a pull request for). We encourage that contributors use linters to check their code for errors and style so this is unlikely to be a big problem.

We understand that it can be difficult to read someone else's code when their style heavily differs from yours. Please be respectful of other's coding styles, in general you are as unlikely to change someone else's coding style as they are to change your own. Again, all that matters is that the style across the codebase is consistent. Consider the following examples:

```python title="Inconsistent codebase"
# Python does not care about the number of spaces used for indentation, some
# languages (like bash) don't care about indentation at all. 
# But the lack of consistency here can make code less readable.
# (especially when there is heavier nesting).
def some_function(i, j, k):
     do_something()

def factorial(n):
  result = 1
  for i in range(1, n + 1):
    result = result * i
    return result

def some_other_function(l, m, n, o, p):
          do_something()
```

```python title="Consistent codebase"
# Consistent style generally makes code more readable. The number of spaces
# used for indentation might be hotly discussed, but we don't care about how 
# many spaces are just used. We only care about consistency across the codebase.
def some_function(i, j, k):
    do_something()

def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result

def some_other_function(l, m, n, o, p):
    do_something()
```

:::info[Nits]
If the reviewer really wants to comment on the code style, it is common practice that the comment begins with "*Nit:*". This way the reviewee is made aware of the changes they *could* make, but they are not pressured into changing their coding style if the change doesn't suit them.

Please be sparing with your useage of *Nit*. Flooding a code review with nits is generally a bad thing, it distracts the developer from the actual problems with the code. If your code review has a comparable number of nits to actually useful comments, please stop. 
:::

## Write something nice

Code reviews all too often focus on the negatives. If you are the reviewer, find something positive to say about the code. Not every comment needs to explain how to improve a function or a specific line. Try starting off a comment with "I like the way you..." or "This is a great way to ... I hope to start using this myself" (*etc.*).

This can make code reviews much more inviting to newer developers, which is of course a good thing.