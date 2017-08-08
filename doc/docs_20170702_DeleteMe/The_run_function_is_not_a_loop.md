---
title: The run function is not a loop
permalink: /The_run_function_is_not_a_loop/
---

In Meep, there are functions `run-until` and similar that run the simulation, and take arguments allowing custom actions to be performed on every time step, or on some subset of the time steps. Many users misunderstand this, however, and make the same mistake: they think the `run` function is a "looping" construct of some kind, and that you can just put any code you want into it and it will get executed for every time step. This mistake and how to correct it are described in this article.

### Hello world

Let's consider a "Hello world" example. Suppose we start with a control file that runs for 200 time units and outputs \(E_z\) on each time step:

`(run-until 200 output-efield-z)`

and now we want to modify it to also print "Hello world!" for every time step, as it is running.

#### The wrong way

Many users will naively write:

`(run-until 200`
`   output-efield-z`
`   (print "Hello world!\n"))`

**This is wrong**. It will output "Hello world!" *once*, then give an error. What is going on? The problem is that you are thinking of `run-until` in the wrong way, as if it were a loop:

` for time ≤ 200 do`
`   output-efield-z`
`   (print "Hello world!\n")`

This is **not** what is happening. Instead, `run-until` is just a *function* that runs the simulation, and its arguments should be *functions* that are called for each time step. That is, it is really doing something like:

` evaluate the arguments: 200: a number`
`                         output-efield-z: a function`
`                         (print "Hello world!\n"): prints output and returns #`<unspecified>
` call run-until:`
`    time-step until t=200`
`    at each time step, call the arguments:`
`       call (output-efield-z)`
`       call (#`<unspecified>`)`

Two things went wrong. First, the arguments are evaluated *before* calling the function, which means that the `print` statement is executed before `run-until` even starts! Second, `run-until` then tries to call the *result* of `(print` `...)` as if it were a function, which causes an error because `(print` `...)` does not return a function. (The `print` returns a special Scheme code `#`<unspecified> that means it doesn't really return anything at all.)

#### The right way

What we *should* have passed to `run-until`, instead of the *result* of calling `(print` `...)`, is a *function* that calls `(print` `...)`. There are two ways to do it.

First, we could explicitly define a function, call it `my-hello`, that does what we want:

`(define (my-hello) (print "Hello world!\n"))`
`(run-until 200 output-efield-z my-hello)`

Notice two things. First, `my-hello` is a function of no arguments, which means that it is just called at every time step (another, more complicated, possibility is described in the [Meep Reference](/Meep_Reference#Writing_your_own_step_functions "wikilink")). Second, when we call `run-until`, we just pass the *name* of the function `my-hello`, and not the *result* of calling the function `(my-hello)`.

A second possibility is that we could use Scheme's `lambda` construct to define our function in-line. The `lambda` syntax in Scheme allows you to define "anonymous" functions without assigning them a name via `define`, and to stick the function definition right into another expression. It works like this:

`(run-until 200 output-efield-z (lambda () (print "Hello world!\n")))`

Here, the `(lambda` `()` `...)` defines a function of no arguments `()` that, when called, executes the `...` statements.

[Category:Meep examples](/Category:Meep_examples "wikilink")