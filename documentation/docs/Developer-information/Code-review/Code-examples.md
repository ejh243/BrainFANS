---
sidebar_position: 3
title: "Examples of good and bad code"
description: "Examples of good and bad code and what to comment on"
---

# Examples of good and bad code

When conducting a code review, there are many aspects the reviewer should be looking out for. These aspects have been laid out in [this list](./Conducting-a-code-review.md#what-to-look-for-in-a-code-review). This page aims to give good and bad examples of code for each of the points on the list. Code is written in python as *generally* python is easy to read for those new to programming.

## Design

This is a rather abstract concept and goes beyond the following example, the full description given in the [list](./Conducting-a-code-review.md#what-to-look-for-in-a-code-review) is more general than the example given.

Suppose you have the file `main.py`. At some point in `main.py`, the factorial function is required (for some calculation). Suppose further that `main.py` is currently using an external package (like numpy) to source this function.

Now assume that a new pull request is made where the contributor has made their own version of the factorial function, claiming that their version is better. This is an example of bad design. If the external, original version of the factorial function is used in other files, accepting this pull request would result in inconsistencies. In this specific example, the reprucussions would not be too bad. However, consider the same example but in a scenario where it wasn't a factorial function but a much more complex library of functions. It is highly unlikely that the newer code written by the contributor is better than the previously used code (that was likely written by a whole team of people). Even in the case the new code is objectively better (faster), it may produce slightly different results to the original code that is used across the rest of the repository, introducing inconsistencies. 


## Functionality

The best way to check functionality is in testing the code yourself. Consider the below code blocks:

```python title="non-functional code"
# This function will not work as 'result' is not initialised before being
# used in calculations
def factorial(n):
    for i in range(1, n + 1):
        result *= i
    return result
```

```python title="Good code"
# The function now outputs the correct value 
# (provided n is a non-negative integer)
def factorial(n):
    for i in range(1, n + 1):
        result *= i
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
        result *= i
    return result
```


## Complexity

Code that is more complex than it needs to be is generally bad.

```python title="highly complex code"
import math

# Note that the gamma function will return the factorial function for 
# x greater than or equal to 1. However it covers all real numbers, not just
# non-negative integers.

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
        result *= i
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
        result *= i
    return result
```

## Readability

Readability is inherently subjective. If the reviewer is unable to read some code (or gets lost), it is likely not their fault, but the code's.

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
        result *= i
    return result
```