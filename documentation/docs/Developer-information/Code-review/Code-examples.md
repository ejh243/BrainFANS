---
sidebar_position: 3
title: "Examples of good and bad code"
description: "Examples of good and bad code and what to comment on"
---

# Examples of good and bad code

When conducting a code review, there are many aspects the reviewer should be looking out for. These aspects have been laid out in [this list](./Conducting-a-code-review.md#what-to-look-for-in-a-code-review). This page aims to give good and bad examples of code for each of the points on the list. Code is written in python as, *generally*, python is easy to read for those new to programming. To keep consistency across code examples, all code will centre around the factorial function. If you have not heard of this function (or want a reminder), please follow [this link](https://en.wikipedia.org/wiki/Factorial).

## Design

Suppose the file `main.py` currently utilizes an external package, such as numpy, to source the factorial function. Now, imagine a developer creates a pull request where the developer has written their own factorial function, claiming that it is superior to the original implementation.

This scenario raises concerns about code design and consistency. If the new factorial function is accepted without considering its impact on other parts of the codebase, inconsistencies may arise. To see this, imagine lots of other files in the codebase are still using the external (numpy) version of the factorial function. Accepting this pull request would result in different files using different versions of the factorial function (either the newly written one or the external numpy one). If the two versions of the function differ in terms of output, input (*etc.*) then the codebase becomes more difficult to work with. When the factorial function is needed in a future file, which one should the developer pick? If another file relies on two separate scripts that each use a different version of the factorial function, will the output be as expected (perhaps the differences between the versions of the function lead to incompatibility)?

While this example focuses on the factorial function, the same principle applies to (and is more relevant for) more complex libraries, functions, classes (*etc.*). Even if the content of the pull request is deemed objectively better—such as being faster—it may produce results that differ from the original code, leading to inconsistencies across the repository.

By carefully evaluating proposed code changes in the context of the existing codebase, we can avoid introducing inconsistencies and maintain the overall stability and coherence of the system.

Note that this is not to say that any fundamental changes to the repository should be rejected. Maybe a 'better' version of a function *is* required for optimisation. The important aspect of 'design' is that the new code will integrate properly with the rest of the codebase. In the example scenario, if all instances of the factorial function across the codebase were changed (and it was genuinely 'better'), then the pull request can still be accepted.

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
    result = 1
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


## Readability/Understandability/Maintainability

These topics all roll into one for the majority of purposes. What they come down to is: Code is not just for the developer who wrote it; but for others to read, use and maintain. If the reviewer is unable to read/understand some code, it will likely be harder for other developers to work on the file in the future. 

```python title="Bad code"
# If you didn't have the preconception that this was the factorial funciton
# would you instantly recognise what this is doing? Would you be comfortable in
# making any changes? Can you spot the bug in this code?
def f(n):
    if isinstance(n, int):
        if n < 0:
            r = 1
            [r := r * i for i in range(1, n+1)]
            return r
        else:
            print("invalid n")
    else:
        print("invalid n")
```

```python title="Better code"
# In general, this is (hopefully) easier to understand
def factorial(n):
    if not isinstance(n, int) or n < 0:
        return None
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```

Note that readability/understandability/maintainability can be subjective in lots of scenarios. For now, if you, the reviewer, understands the code to a 'good enough' degree (and you are not completely lost), the code is fine and does not warrant a refactoring. There is [no such thing as perfect code](./Best-practices.md#there-is-no-perfect-code). It is very easy to go back and forth on what makes code cleaner and easier to read, but all this really ends up doing is wasting time. If the proposed changes demonstrate reasonable readability and improve the existing codebase, approve the request.

:::info[Clean code]
Clean code is a set of programming principles that aims for consistent naming conventions, function structure and minimal complexity. To a degree, clean code is great for achieving the points on this list. But it also comes with the drawbacks of performance hits and sometimes (paradoxically) decreased readability and understandability. A developer can go too far the other way with readability/understandability/maintainability. Strive for 'good enough', a compromise between easy to read and easy to write.
:::

## Comments

A good rule of thumb is: "comments explain *why, not what*". To see the difference, the examples below have the same code, but the first comments the *what* and the seconf the *why*.

```python title="Commenting on what"
# Define the factorial function
def factorial(n):
    # Checks if the input is not an integer or less than zero
    if not isinstance(n, int) or n < 0:
        # Returns None if the input is not an integer or less than zero
        return None
    # Initialises 'result' to be 1
    result = 1
    # Creates a loop that starts at index 1 and ends at index n
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

Sometimes, explaining *what* the code does feel like a necessity. In such cases we refer back to our stance on [complex code](#complexity). If the code is complex enough that *what*-based comments are required, then the code would likely benefit from being refactored. The exception to this is circumstances like regular expressions. If the code uses something like regex, then the comments should explain what the code is doing. Computers read regex, humans do not (well, not easily).

## Scalability/Expandability

Can the code be repurposed in a new context? We appreciate that spotting scalability problems is usually rather difficult in more complex codebases without lots of domain knowledge. Generally, you should be looking to see if some code can be modularised so that it can be repurposed elsewhere.

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

:::info[DRY]
A common term you might hear in programming/software engineering spaces is DRY:
> The *DRY* principle: *Don't repeat yourself*.

This principle comes from the idea that: if you repeat yourself a lot in code, making changes in the future becomes much more cumbersome. Taking on this principle can help a lot with making code scalable as it often results in more modular code.

However, one can go too far with DRY and end up with difficult to understand code. You do not *need* code to be heavily abstracted wherever possible. The DRY principle should be used when the repetition becomes a problem, not when it "might become a problem".
:::

Alongside repetition, another aspect of code that can increase scalability is: "Functions/Classes (constructs) should do *one* thing". This is also known as the 'single-responsibility principle'. Functions that do lots of things are harder to reimplement elsewhere. Functions that do a single thing are more modular and can be repurposed. Consider the two following examples:

```python title="Bad code"
# This seems harmless enough, you call the factorial function and it outputs the
# answer. However, this function is doing TWO things. First, it evaluates the 
# factorial then second, it prints the answer. What if we just wanted the
# answer and not the printed sentence? This function is not very modular.
def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result = result * i
    print("The factorial of", n, "is", result)
```

```python title="Better code"
# This function now does one thing, it returns the value of n!.
# If we wanted to print the answer in a sentence, we could achieve this with
# a different line of code (or function) outside of the factorial function.
# This way the function can be used for a variety of purposes.
def factorial(n):
    result = 1
    for i in range(1, n + 1):
        result = result * i
    return result
```