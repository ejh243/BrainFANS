---
sidebar_position: 3
title: "Examples of good and bad code"
description: "Examples of good and bad code and what to comment on"
---

# Examples of good and bad code

When conducting a code review, there are many aspects the reviewer should be looking out for. These aspects have been laid out in [this list](./Conducting-a-code-review.md#what-to-look-for-in-a-code-review). This page aims to give good and bad examples of code for each of the points on the list. Code is written in python as, *generally*, python is easy to read for those new to programming. To keep consistency across code examples, all code will centre around the factorial function. If you have not heard of this function (or want a reminder), please follow [this link](https://en.wikipedia.org/wiki/Factorial).

## Design

This is a rather abstract concept and goes beyond the following example, the full description given in the [list](./Conducting-a-code-review.md#what-to-look-for-in-a-code-review) is more general than the example given.

Suppose you have the file: `main.py`. At some point in `main.py`, the factorial function is required (for some calculation). Suppose further that `main.py` is currently using an external package (like numpy) to source this function.

Now assume that a new pull request is made where the contributor has made their own version of the factorial function, claiming that their version is better. This is an example of bad design. If the external, original version of the factorial function is used in other files, accepting this pull request would result in inconsistencies. In this specific example, the reprucussions would not be too bad. However, consider the same example but in a scenario where it wasn't a factorial function but a much more complex library of functions. It is highly unlikely that the newer code (written by the contributor) is better than the previously used code (that was probably written by a whole team of clever people). Even in the case the new code is objectively better (*e.g.* faster), it may produce slightly different results to the original code that is used across the rest of the repository, introducing inconsistencies. 


## Functionality

The best way to check functionality is in testing the code yourself. Consider the below code blocks:

```python title="non-functional code"
# This function will not work as 'result' is not initialised before being
# used in calculations
def factorial(n):
    for i in range(1, n + 1):
        result = result * i
    return result
```

```python title="Good code"
# The function now outputs the correct value but does not handle edge cases
# such as the user inputting a letter instead of a number etc.
def factorial(n):
    for i in range(1, n + 1):
        result = result * i
    return result
```

```python title="Better code"
# The function now outputs the correct value and accounts for the case
# where the user does not input a non-negative integer
def factorial(n):
    if not isinstance(n, int) or n < 0:
        return None
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```


## Complexity

Code that is more complex than it needs to be is generally bad.

```python title="highly complex code"
import math

# Don't worry about what this code does. In a nutshell the gamma function will 
# return the factorial function for x greater than or equal to 1. 
# However it is defined for all real numbers, not just non-negative integers.

# Also note that this approximation of the gamma function isn't even very
# accurate.

def gamma_function(x):
    if x == 1:
        return 1
    elif x < 0.5:
        return math.pi / (math.sin(math.pi * x) * gamma_function(1 - x))
    else:
        x -= 1
        a = [0.99999999999980993, 676.5203681218851, -1259.1392167224028,
             771.32342877765313, -176.61502916214059, 12.507343278686905,
             -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7]
        t = a[0]
        for i in range(1, len(a)):
            t += a[i] / (x + i)
        return math.sqrt(2 * math.pi) * (x + 4.5) ** (x - 0.5) * math.exp(-x - 4.5) * t
```

```python title="Better code"
# For most purposes, the below code does the exact same thing as the 'bad' code
# Considering the below code is much less complex, it is viewed as better
def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```


## Understandability

Understandability is inherently subjective. If the reviewer is unable to understand some code, it is likely not their fault, but the code's.

```python title="Bad code"
# This code is still doing the same computation, but the loop is now
# much less understandable
def factorial(n):
    result = 1
    [result := result * i for i in range(1, n+1)]
    return result
```

```python title="Better code"
# In general, this is easier to understand
def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```

## Readability

Readability is inherently subjective. If the reviewer is unable to read some code (or gets lost), it is likely not their fault, but the code's. 

Typically understandability and readability go hand in hand, so a combinatation of the two should be considered.

```python title="Bad code"
# In general, this is hard to read as one does not have any context for what
# f, n, i and r are
def f(n):
    r = 1
    for i in range(1, n + 1):
        r *= i
    return r
```

```python title="Better code"
# In general, this is easier to understand
def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```

:::info[too far the other way]
It should be noted that you can go too far the other way with 'readability' when it comes to naming code objects. Using too many words or being overly descriptive can paradoxically decrease code readability and understandability.
\
The classic 'joke' is:
> There are only two hard things in Computer Science: cache invalidation and naming things. - Phil Karlton
:::

## Maintainability

Code that is maintainable would be easy for someone else to come along to the script (or source code) and make ammendments to it. We are going to stick with the factorial function here for consistency across examples. In order for this to make sense, we are going to reintroduce the error checking that was seen in [functionality](#functionality). Recall that factorial is not defined for strings, non-integers and negative integers. Below we will make two functions that implement error handling in two ways.

```python title="Bad code"
# This is not only hard to read, but would you be comfortable changing this
# code? In this scenario it is hard to follow the indentation and so making
# changes might break the logic if you are not careful.

# You can imagine that if the function was doing something more complex than
# simply calculating the factorial, this would be even harder to maintain.
def factorial(n):
    if isinstance(n, str):
        if isinstance(n, int):
            if n>=0:
                result = 1
                for i in range(1, n + 1):
                    result = result * i
                return result
            else:
                print("input must be non-negative")
        else:
            print("input must be an integer")
    else:
        print("input must be an integer")
```

```python title="Better code"
# This code uses negation to remove the edge cases first, instead of nesting
# the code deep into the function. With comments, this code might be even better
def factorial(n):
    if not isinstance(n, int) or n < 0:
        return None
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```

Maintainability is not limited to deeply nested code, this is just one example of how code can be difficult to maintain due to it being poorly written.

## Comments

A good rule of thumb is comments explain *why, not what*. To see the difference, the examples below have the same code, but the comments will either explain *what* or *why*.

```python title="Commenting on what"
# Define the factorial function
def factorial(n):
    # Checks if the input is not an integer or less than zero
    if not isinstance(n, int) or n < 0:
        # Returns None if the input is not an integer or less than zero
        return None
    # Initialises 'result' to be 1
    result = 1
    # Creates a loop that starts at index 1 and ends at index n+1
    for i in range(1, n + 1):
        # Multiplies 'result' by the current index of the for loop
        result = result * i
    # returns the result of the function
    return result
```

```python title="Commenting on why"
def factorial(n):
    # The factorial function is not defined for non-negative integers,
    # error handling is required in case the user does not enter such. 
    if not isinstance(n, int) or n < 0:
        return None
    result = 1
    for i in range(1, n + 1):
        # The factorial function is defined as: "The product of all integers
        # up to the input".
        result = result * i
    return result
```

The above is obviously an egregious use of comments. Python is readable enough that none of the above comments are particularly required. However, hopefully you can see that explaining "*what*" really does not help the reader.

Sometimes, explaining *what* the code does feel like a necessity. In such cases we refer back to our stance on [complex code](#complexity). If the code is complex enough that *what*-based comments are required, then the code would likely benefit from being refactored. The exception to this is *regular expressions*. If the code uses regex, then the comments should explain what the code is doing. Computers read regex, humans do not (well, not easily).

## Scalability/Expandability

Can the code be repurposed in a new context? Again, we are going to stick with the factorial function for simplicity, spotting scalability problems is usually rather difficult (without lots of domain knowledge).

```python title="Bad code"
# ... Code ... #
# ... Code ... #

factorial = 1
for i in range(1, 6):
    factorial = factorial * i
print(factorial)


# ... Code ... #
# ... Code ... #

factorial = 1
for i in range(1, 11):
    factorial = factorial * i
print(factorial)
```

The above code is the same factorial function, but it is just being used inside the codebase for specific scenarios (here, it calculates 5! and 10!). This is not desirable as we can't reuse the factorial in other parts of the codebase easily. Here, you would need to copy and paste the code block and know that the `range` needed to be changed to n+1 (to calculate n!).

```python title="Better code"
# This function can be called from anywhere with 'factorial(n)' and is more
# versatile
def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```

:::info[Going too far the other way]
Sometimes a line of code is very niche and is only required once in an entire codebase. You do not need to bundle every line of code into a function or a class method. 
> A good rule of thumb is the *DRY* principle. That is *Don't repeat yourself*.

If code is repeated across a codebase (or *will* be repeated in the **near** future), making such code scalable is beneficial for satisfying the other points on this list.
:::

## Style

As mentioned [here](./Conducting-a-code-review.md#style), style is usually down to personal preference and it can sometimes be difficult to separate views on style from readability. Style should only be taken into account if it helps with consistency across the codebase.

```python title="Style inconsistencies"
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

Above is a simple case where indentation changes throughout the codebase. Python does not actually care how many spaces are used for indentation (so the user could use 1 space, 2 spaces, 4, 5 *etc.*). No number of spaces is technically correct, but consistency here makes the codebase more readable.

```python title="Better code"
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

The above is certainly more readable despite the individual lines of code being exactly the same between examples.

## Documentation

As mentioned [here](./Conducting-a-code-review.md#documentation), documentation is not really a requirement in a pull request. The relevant documentation can be changed at a later point in time if needs be. What does matter is the internal documentation of scripts (or source code). What we mean by this is a description of what a file does in the file preamble, or what a function does in its docstring. These are generally very useful to include in a codebase as it helps dramastically with its future use. In the event that a user or developer wants to find the correct file for their needs, high level descriptions of files, classes and functions will likely be incredibly useful to them.

Suppose that our factorial function is amongst a sea of other useful mathematical functions in a bigger python file. If the top of said file included a preamble, explaining the purpose of the file and who created it, future developers will have a much easier time working with the codebase. If no preamble was given, the factorial function might not be seen by future developers and they may end up creating the function themselves (which can be a waste of time with more complex functions).