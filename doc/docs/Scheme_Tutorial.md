---
Scheme Tutorial
---

In this section, we'll walk through the process of computing the band structure and outputting some fields for a two-dimensional photonic crystal using MPB. This should give you the basic idea of how it works and some of the things that are possible. Here, we tell you the truth, but not the whole truth. The [User Interface](Scheme_User_Interface.md) section gives a more complete, but less colloquial, description of the features supported by Photonic-Bands. See also the [data analysis tutorial](Data_Analysis_Tutorial.md) for more examples, focused on analyzing and visualizing the results of MPB calculations.

[TOC]

The ctl File
------------

The use of the MPB revolves around the control file, abbreviated "ctl" and typically called something like `foo.ctl` (although you can use any filename extension you wish). The ctl file specifies the geometry you wish to study, the number of eigenvectors to compute, what to output, and everything else specific to your calculation. Rather than a flat, inflexible file format, however, the ctl file is actually written in a scripting language. This means that it can be everything from a simple sequence of commands setting the geometry, etcetera, to a full-fledged program with user input, loops, and anything else that you might need.

Don't worry, though—simple things are simple and even there you will appreciate the flexibility that a scripting language gives you. (e.g. you can input things in any order, without regard for whitespace, insert comments where you please, omit things when reasonable defaults are available...)

The ctl file is actually implemented on top of the libctl library, a set of utilities that are in turn built on top of the Scheme language. Thus, there are three sources of possible commands and syntax for a ctl file:

-   [Scheme](https://en.wikipedia.org/wiki/Scheme_(programming_language)), a powerful and beautiful programming language developed at MIT, which has a particularly simple syntax: all statements are of the form `(function arguments...)`. We run Scheme under the GNU Guile interpreter which is designed to be plugged into programs as a scripting and extension language. You don't need to learn much Scheme for a basic ctl file, but it is always there if you need it. You can learn more about [Guile and Scheme](Guile_and_Scheme_Information.md).
-   libctl, a library that we built on top of Guile to simplify communication between Scheme and scientific computation software. libctl sets the basic tone of the user interface and defines a number of useful functions. See the [libctl](https://libctl.readthedocs.io) pages.
-   MPB which defines all the interface features that are specific to photonic band structure calculations. This manual is primarily focused on documenting these features.

It would be an excellent idea at this point for you to go read the [libctl manual](https://libctl.readthedocs.io), particularly the [Basic User Experience](https://libctl.readthedocs.io/en/latest/Libctl_Basic_User_Experience/), which will give you an overview of what the user interface is like, provide a crash course in the Scheme features that are most useful here, and describe some useful general features. We're not going to repeat this material much, so learn it now!

Okay, let's continue with our tutorial. MPB is normally invoked by running something like:

```
unix% mpb foo.ctl >& foo.out
```

which reads the ctl file `foo.ctl` and executes it, saving the output to the file `foo.out`. Some sample ctl files are provided in the `mpb-ctl/examples/` directory. However, as you already know (since you obediently read the libctl manual, right?), if you invoke `mpb` with no arguments, you are dropped into an *interactive* mode in which you can type commands and see their results immediately. Why don't you do that right now, in your terminal window? Then, you can paste in the commands from the tutorial as you follow it and see what they do.

Our First Band Structure
------------------------

As our beginning example, we'll compute the band structure of a two-dimensional square lattice of dielectric rods in air. See [our online textbook](http://ab-initio.mit.edu/book), ch. 5. In our control file, we'll first specify the parameters and geometry of the simulation, and then tell it to run and give us the output.

All of the parameters, each of which corresponds to a Scheme variable, have default setting, so we only need to specify the ones we need to change. For a complete listing of the parameter variables and their current values, along with some other info, type `(help)` at the `guile>` prompt. One of the parameters, `num-bands`, controls how many bands (eigenstates) are computed at each k point. If you type `num-bands` at the prompt, it will return the current value, `1`--this is too small; let's set it to a larger value:

```
(set! num-bands 8)
```

This is how we change the value of variables in Scheme. If you type `num-bands` now, it will return `8`. The next thing we want to set although the order really doesn't matter, is the set of k points (Bloch wavevectors) we want to compute the bands at. This is controlled by the variable `k-points`, a list of 3-vectors which is initially empty. We'll set it to the corners of the irreducible Brillouin zone of a square lattice, Gamma, X, M, and Gamma again:

```
(set! k-points (list (vector3 0 0 0)     ; Gamma
                     (vector3 0.5 0 0)   ; X
                     (vector3 0.5 0.5 0) ; M
                     (vector3 0 0 0)))   ; Gamma
```

Notice how we construct a list, and how we make 3-vectors; notice also how we can break things into multiple lines if we want, and that a semicolon (';') marks the start of a comment. Typically, we'll want to also compute the bands at a lot of intermediate k points, so that we see the continuous band structure. Instead of manually specifying all of these intermediate points, however, we can just call one of the functions provided by libctl to interpolate them for us:

```
(set! k-points (interpolate 4 k-points))
```

This takes the `k-points` and linearly interpolates four new points between each pair of consecutive points. If we type `k-points` now at the prompt, it will show us all 16 points in the new list:

```
(#(0 0 0) #(0.1 0.0 0.0) #(0.2 0.0 0.0) #(0.3 0.0 0.0) #(0.4 0.0 0.0) #(0.5 0 0) #(0.5 0.1 0.0) #(0.5 0.2 0.0) #(0.5 0.3 0.0) #(0.5 0.4 0.0) #(0.5 0.5 0) #(0.4 0.4 0.0) #(0.3 0.3 0.0) #(0.2 0.2 0.0) #(0.1 0.1 0.0) #(0 0 0))
```

Alternatively, you can use `(set! k-points (kinterpolate-uniform 4 k-points))` to interpolate points that are roughly uniformly spaced in <b>k</b> space (i.e. it will use a variable number of points between each pair of vectors to keep the spacing roughly equal).

As is [described below](Scheme_Tutorial.md#a-few-words-on-units), all spatial vectors in the program are specified in the basis of the lattice directions normalized to `basis-size` lengths. The default is unit-normalized. The k points are specified in the basis of the (unnormalized) reciprocal lattice vectors. In this case, we don't have to specify the lattice directions, because we are happy with the defaults--the lattice vectors default to the Cartesian unit axes (i.e. the most common case, a square/cubic lattice). The reciprocal lattice vectors in this case are also the unit axes. We'll see how to change the lattice vectors in later subsections.

Now, we want to set the geometry of the system--we need to specify which objects are in the primitive cell of the lattice, centered on the origin. This is controlled by the variable `geometry`, which is a list of geometric objects. As you know from reading the libctl documentation, objects, which are more complicated, structured data types, are created by statements of the form `(make type (property1 value1) (property2 value2) ...)`. There are various kinds (sub-classes) of geometric object: cylinders, spheres, blocks, ellipsoids, and perhaps others in the future. Right now, we want a square lattice of rods, so we put a single dielectric cylinder at the origin:

```
(set! geometry (list (make cylinder 
                       (center 0 0 0) (radius 0.2) (height infinity)
                       (material (make dielectric (epsilon 12))))))
```

Here, we've set several properties of the cylinder: the `center` is the origin, its `radius` is 0.2, and its `height` (the length along its axis) is `infinity`. Another property, the `material`, it itself an object--we made it a dielectric with the property that its `epsilon` is 12. There is another property of the cylinder that we can set, the direction of its axis, but we're happy with the default value of pointing in the z direction.

All of the geometric objects are ostensibly three-dimensional, but since we're doing a two-dimensional simulation the only thing that matters is their intersection with the xy plane (z=0). Speaking of which, let us set the dimensionality of the system. Normally, we do this when we define the size of the computational cell, controlled by the `geometry-lattice` variable, an object of the `lattice` class: we can set some of the dimensions to have a size `no-size`, which reduces the dimensionality of the system.

```
(set! geometry-lattice (make lattice (size 1 1 no-size)))
```

Here, we define a 1x1 two-dimensional cell (defaulting to square). This cell is *discretized* according to the `resolution` variable, which defaults to `10` (pixels/lattice-unit). That's on the small side, and this is only a 2d calculation, so let's increase the resolution:

```
(set! resolution 32)
```

This results in a 32x32 computational grid. For efficient calculation, it is best to make the grid sizes a power of two, or factorizable into powers of small primes (such as 2, 3, 5 and 7). As a rule of thumb, you should use a resolution of at least `8` in order to obtain reasonable accuracy.

Now, we're done setting the parameters--there are other parameters, but we're happy with their default values for now. At this point, we're ready to go ahead and compute the band structure. The simplest way to do this is to type `(run)`. Since this is a two-dimensional calculation, however, we would like to split the bands up into TE- and TM-polarized modes, and we do this by invoking `(run-te)` and `(run-tm)`.

These produce a lot of output, showing you exactly what the code is doing as it runs. Most of this is self-explanatory, but we should point out one output in particular. Among the output, you should see lines like:

```
tefreqs:, 13, 0.3, 0.3, 0, 0.424264, 0.372604, 0.540287, 0.644083, 0.81406, 0.828135, 0.890673, 1.01328, 1.1124
```

These lines are designed to allow you to easily extract the band-structure information and import it into a spreadsheet for graphing. They comprise the k point index, the k components and magnitude, and the frequencies of the bands, in comma-delimited format. Each line is prefixed by "tefreqs" for TE bands, "tmfreqs" for TM bands, and "freqs" for ordinary bands produced by `(run)`. Using this prefix, you can extract the data you want from the output by passing it through a program like `grep`. For example, if you had redirected the output to a file `foo.out` as described earlier, you could extract the TM bands by running `grep tmfreqs foo.out` at the terminal prompt. Note that the output includes a header line, like:

```
tefreqs:, k index, kx, ky, kz, kmag/2pi, band 1, band 2, band 3, band 4, band 5, band 6, band 7, band 8
```

explaining each column of data. Another output of the `run` is the list of band gaps detected in the computed bands. For example the `(run-tm)` output includes the following gap output:

```
Gap from band 1 (0.282623311147724) to band 2 (0.419334798706834), 38.9514660888911%
Gap from band 4 (0.715673834754345) to band 5 (0.743682920649084), 3.8385522650349%
```

This data is also stored in the variable `gap-list`, which is a list of `(gap-percent gap-min gap-max)` lists. It is important to realize, however, that this band-gap data may include "false positives," from two possible sources:

-   If two bands cross, a false gap may result because the code computes the gap by assuming that bands never cross. Such false gaps are typically quite small (&lt; 1%). To be sure of what's going on, you should either look at the symmetry of the modes involved or compute k points very close to the crossing. Although even if the crossing occurs precisely at one of your k-points, there usually won't be an exact degeneracy for numerical reasons.
-   One typically computes band diagrams by considering k-points around the boundary of the irreducible Brillouin zone. It is possible, though rare, that the band edges may occur at points in the interior of the Brillouin zone. To be absolutely sure you have a band gap and of its size, you should compute the frequencies for points inside the Brillouin zone, too.

You've computed the band structure, and extracted the eigenfrequencies for each k point. But what if you want to see what the fields look like, or check that the dielectric function is what you expect? To do this, you need to output [HDF5 files](https://support.hdfgroup.org/HDF5/) for these functions. HDF5 is a binary format for multi-dimensional scientific data, and can be read by many visualization programs. The output files have filenames with suffixes ".h5".

When you invoke one of the `run` functions, the dielectric function in the unit cell is automatically written to the file `epsilon.h5`. To output the fields or other information, you need to pass one or more arguments to the `run` function. For example:

```
(run-tm output-efield-z)
(run-te (output-at-kpoint (vector3 0.5 0 0) output-hfield-z output-dpwr))
```

This will output the electric (E) field z components for the TM bands at all k-points; and the magnetic (H) field z components and electric field energy densities (D power) for the TE bands at the X point only. The output filenames will be things like `e.k12.b03.z.te.h5`, which stands for the z component (`.z`) of the TE (`.te`) electric field (`e`) for the third band (`.b03`) of the twelfth k point (`.k12`). Each HDF5 file can contain multiple datasets. In this case, it will contain the real and imaginary parts of the field (in datasets "z.r" and "z.i"), and in general it may include other data too (e.g. `output-efield` outputs all the components in one file). See also the [Data Analysis Tutorial](Data_Analysis_Tutorial.md).

There are several other output functions you can pass, described in the [User Interface](Scheme_User_Interface.md), like `output-dfield`, `output-hpwr`, and `output-dpwr-in-objects`. Actually, though, you can pass in arbitrary functions that can do much more than just output the fields--you can perform arbitrary analyses of the bands using functions that we will describe later.

Instead of calling one of the `run` functions, it is also possible to call lower-level functions of the code directly, to have a finer control of the computation. Such functions are described in the reference section.

A Few Words on Units
--------------------

In principle, you can use any units you want with MPB. Maxwell's equations possess an important property--they are *scale-invariant*. See [our online textbook](http://ab-initio.mit.edu/book), ch. 2. If you multiply all of your sizes by 10, the solution scales are simply multiplied by 10 likewise while the frequencies are divided by 10. So, you can solve a problem once and apply that solution to all length-scales. For this reason, we usually pick some fundamental lengthscale *a* of a structure, such as its lattice constant (unit of periodicity), and write all distances in terms of that. That is, we choose units so that *a* is unity. Then, to apply to any physical system, one simply scales all distances by *a*. This is what we have done in the preceding and following examples. This is the default behavior of MPB: the lattice constant is one, and all coordinates are scaled accordingly.

As has been mentioned already, nearly all 3-vectors in the program are specified in the *basis* of the lattice vectors *normalized* to lengths given by `basis-size`, defaulting to the unit-normalized lattice vectors. That is, each component is multiplied by the corresponding basis vector and summed to compute the corresponding Cartesian vector. It is worth noting that a basis is not meaningful for scalar distances such as the cylinder radius. These are just the ordinary cartesian distances in your chosen units of *a*.

Note also that the [k-points](Scheme_User_Interface.md#input-variables), as mentioned above, are an exception: they are in the basis of the reciprocal lattice vectors. See [our online textbook](http://ab-initio.mit.edu/book), appendix B. If a given dimension has size `no-size`, its reciprocal lattice vector is taken to be 2&#960;/*a*.

We provide [conversion functions](Scheme_User_Interface.md#coordinate-conversion-functions) to transform vectors between the various bases.

The frequency eigenvalues returned by the program are in units of *c/a*, where *c* is the speed of light and *a* is the unit of distance. Thus, the corresponding vacuum wavelength is *a* over the frequency eigenvalue.

Bands of a Triangular Lattice
-----------------------------

As a second example, we'll compute the TM band structure of a *triangular* lattice of dielectric rods in air. To do this, we only need to change the lattice, controlled by the variable `geometry-lattice`. We'll set it so that the first two basis vectors (the properties `basis1` and `basis2`) point 30 degrees above and below the x axis, instead of their default value of the x and y axes:

```
(set! geometry-lattice (make lattice (size 1 1 no-size)
                        (basis1 (/ (sqrt 3) 2) 0.5)
                        (basis2 (/ (sqrt 3) 2) -0.5)))
```

We don't specify `basis3`, keeping its default value of the z axis. Notice that Scheme supplies us with all the ordinary arithmetic operators and functions, but they use prefix (Polish) notation, in Scheme fashion. The `basis` properties only specify the directions of the lattice basis vectors, and not their lengths--the lengths default to unity, which is fine here.

The irreducible Brillouin zone of a triangular lattice is different from that of a square lattice, so we'll need to modify the `k-points` list accordingly:

```
(set! k-points (list (vector3 0 0 0)          ; Gamma
                     (vector3 0 0.5 0)        ; M
                     (vector3 (/ -3) (/ 3) 0) ; K
                     (vector3 0 0 0)))        ; Gamma
(set! k-points (interpolate 4 k-points))
```

Note that these vectors are in the basis of the new reciprocal lattice vectors, which are different from before. Notice also the Scheme shorthand `(/` `3)`, which is the same as `(/` `1` `3)` or 1/3.

All of the other parameters (`geometry`, `num-bands`, and `grid-size`) can remain the same as in the previous subsection, so we can now call `(run-tm)` to compute the bands. As it turns out, this structure has an even larger TM gap than the square lattice:

```
Gap from band 1 (0.275065617068082) to band 2 (0.446289918847647), 47.4729292989213%
```

Maximizing the First TM Gap
---------------------------

We will now show you a more sophisticated example utilizing the programming capabilities of Scheme. We will write a script to choose the cylinder radius that maximizes the first TM gap of the triangular lattice of rods from above. All of the Scheme syntax here won't be explained, but this should give you a flavor of what is possible.

First, we will write the function that want to maximize, a function that takes a dielectric constant and returns the size of the first TM gap. This function will change the geometry to reflect the new radius, run the calculation, and return the size of the first gap:

```
(define (first-tm-gap r)
  (set! geometry (list (make cylinder
                        (center 0 0 0) (radius r) (height infinity)
                        (material (make dielectric (epsilon 12))))))
  (run-tm)
  (retrieve-gap 1)) ; return the gap from TM band 1 to TM band 2
```

We'll leave most of the other parameters the same as in the previous example, but we'll also change `num-bands` to 2, since we only need to compute the first two bands:

```
(set! num-bands 2)
```

In order to distinguish small differences in radius during the optimization, it might seem that we have to increase the grid resolution, slowing down the computation. Instead, we can simply increase the *mesh* resolution. This is the size of the mesh over which the dielectric constant is averaged at each grid point, and increasing the mesh size means that the average index better reflects small changes in the structure.

```
(set! mesh-size 7) ; increase from default value of 3
```

Now, we're ready to maximize our function `first-tm-gap`. We could write a loop to do this ourselves, but libctl provides a built-in function `(maximize function tolerance arg-min arg-max)` to do it for us using Brent's algorithm. So, we just tell it to find the maximum, searching in the range of radii from 0.1 to 0.5, with a tolerance of 0.1:

```
(define result (maximize first-tm-gap 0.1 0.1 0.5))
(print "radius at maximum: " (max-arg result) "\n")
(print "gap size at maximum: " (max-val result) "\n")
```

`print` is a function defined by libctl to apply the built-in `display` function to zero or more arguments. After five iterations, the output is:

```
radius at maximum: 0.176393202250021
gap size at maximum: 48.6252611051049
```

The tolerance of 0.1 that we specified means that the true maximum is within 0.1 \* 0.176393202250021, or about 0.02, of the radius found here. It doesn't make much sense here to specify a lower tolerance, since the discretization of the grid means that the code can't accurately distinguish small differences in radius.

Before we continue, let's reset `mesh-size` to its default value:

```
(set! mesh-size 3) ; reset to default value of 3
```

A Complete 2D Gap with an Anisotropic Dielectric
------------------------------------------------

As another example, one which does not require so much Scheme knowledge, let's construct a structure with a complete 2D gap (i.e., in both TE and TM polarizations), in a somewhat unusual way: using an [dielectric](Scheme_User_Interface.md#material-type) structure. An anisotropic dielectric presents a different dielectric constant depending upon the direction of the electric field, and can be used in this case to make the TE and TM polarizations "see" different structures.

We already know that the triangular lattice of rods has a gap for TM light, but not for TE light. The dual structure, a triangular lattice of holes, has a gap for TE light but not for TM light at least for the small radii we will consider. Using an anisotropic dielectric, we can make both of these structures simultaneously, with each polarization seeing the structure that gives it a gap.

As before, our `geometry` will consist of a single cylinder, this time with a radius of 0.3, but now it will have an epsilon of 12 (dielectric rod) for TM light and 1 (air hole) for TE light:

```
(set! geometry (list (make cylinder
                       (center 0 0 0) (radius 0.3) (height infinity)
                       (material (make dielectric-anisotropic
                                   (epsilon-diag 1 1 12))))))
```

Here, `epsilon-diag` specifies the diagonal elements of the dielectric tensor. The off-diagonal elements specified by `epsilon-offdiag` default to zero and are only needed when the principal axes of the dielectric tensor are different from the Cartesian xyz axes.

The background defaults to air, but we need to make it a dielectric (epsilon of 12) for the TE light, so that the cylinder forms a hole. This is controlled via the `default-material` variable:

```
(set! default-material (make dielectric-anisotropic (epsilon-diag 12 12 1)))
```

Finally, we'll increase the number of bands back to eight and run the computation:

```
(set! num-bands 8)
(run) ; just use run, instead of run-te or run-tm, to find the complete gap
```

The result, as expected, is a complete band gap:

```
Gap from band 2 (0.223977612336924) to band 3 (0.274704473679751), 20.3443687933601%
```

If we had computed the TM and TE bands separately, we would have found that the lower edge of the complete gap in this case comes from the TM gap, and the upper edge comes from the TE gap.

Finding a Point-Defect State
----------------------------

Here, we consider the problem of finding a point-defect state in our square lattice of rods. This is a state that is localized in a small region by creating a point defect in the crystal — e.g., by removing a single rod. The resulting mode will have a frequency within, and be confined by, the gap. See [our online textbook](http://ab-initio.mit.edu/book), ch. 5.

To compute this, we need a supercell of bulk crystal, within which to put the defect — we will use a 5x5 cell of rods. To do this, we must first increase the size of the lattice by five, and then add all of the rods. We create the lattice by:

```
(set! geometry-lattice (make lattice (size 5 5 no-size)))
```

Here, we have used the default orthogonal basis, but have changed the size of the cell. To populate the cell, we could specify all 25 rods manually, but that would be tedious. A better approach would be to write a loop but in fact this has already been done for you. MPB provides a function, `geometric-objects-lattice-duplicates`, that duplicates a list of objects over the lattice:

```
(set! geometry (list (make cylinder
                       (center 0 0 0) (radius 0.2) (height infinity)
                       (material (make dielectric (epsilon 12))))))
(set! geometry (geometric-objects-lattice-duplicates geometry))
```

There, now the `geometry` list contains 25 rods — the original `geometry` list, which contained one rod, duplicated over the 5x5 lattice.

To remove a rod, we'll just add another rod in the center of the cell with a dielectric constant of 1. When objects overlap, the later object in the list takes precedence, so we have to put the new rod at the end of `geometry`:

```
(set! geometry (append geometry 
                      (list (make cylinder (center 0 0 0) 
                                  (radius 0.2) (height infinity)
                                  (material air)))))
```

Here, we've used the Scheme `append` function to combine two lists, and have also snuck in the predefined material type `air` which has an epsilon of 1.

We'll be frugal and use only 16 points per lattice unit, resulting in an 80x80 grid, instead of the 32 from before:

```
(set! resolution 16)
```

Only a single k point is needed for a point-defect calculation which, for an infinite supercell, would be independent of k:

```
(set! k-points (list (vector3 0.5 0.5 0)))
```

Unfortunately, for a supercell the original bands are folded many times over, in this case, 25 times, so we need to compute many more bands to reach the same frequencies:

```
(set! num-bands 50)
```

At this point, we can call `(run-tm)` to solve for the TM bands. It will take several seconds to compute. Recall that the gap for this structure was for the frequency range 0.2812 to 0.4174. The bands of the solution include exactly one state in this frequency range: band 25, with a frequency of 0.378166. This is exactly what we should expect--the lowest band was folded 25 times into the supercell Brillouin zone, but one of these states was pushed up into the gap by the defect.

We haven't yet output any of the fields, but we don't have to repeat the run to do so. The fields from the last k-point computation remain in memory and can continue to be accessed and analyzed. For example, to output the electric field z component of band 25, we just do:

```
(output-efield-z 25)
```

That's right, the output functions that we passed to `(run)` in the first example are just functions of the band index that are called on each band. We can do other computations too, like compute the fraction of the electric field energy near the defect cylinder within a radius 1.0 of the origin:

```
(get-dfield 25)  ; compute the D field for band 25
(compute-field-energy)  ; compute the energy density from D
(print
 "energy in cylinder: "
 (compute-energy-in-objects (make cylinder (center 0 0 0)
                                  (radius 1.0) (height infinity)
                                  (material air)))
 "\n")
```

The result is 0.624794702341156, or over 62% of the field energy in this localized region; the field decays exponentially into the bulk crystal. The full range of available functions is described in the [User Interface](Scheme_User_Interface.md), but the typical sequence is to first load a field with a `get-` function and then to call other functions to perform computations and transformations on it.

Note that the `compute-energy-in-objects` returns the energy fraction, but does not itself print this value. This is fine when you are running interactively, in which case Guile always displays the result of the last expression, but when running as part of a script you need to explicitly print the result as we have done above with the `print` function. The "\n" string is newline character (like in C), to put subsequent output on a separate line.

Instead of computing all those bands, we can instead take advantage of a special feature of MPB that allows you to compute the bands closest to a "target" frequency, rather than the bands with the lowest frequencies. One uses this feature by setting the `target-freq` variable to something other than zero (e.g. the mid-gap frequency). In order to get accurate results, it's currently also recommended that you decrease the `tolerance` variable, which controls when convergence is judged to have occurred, from its default value of `1e-7`:

```
(set! num-bands 1)  ; only need to compute a single band, now!
(set! target-freq (/ (+ 0.2812 0.4174) 2))
(set! tolerance 1e-8)
```

Now, we just call `(run-tm)` as before. Convergence requires more iterations this time, both because we've decreased the tolerance and because of the nature of the eigenproblem that is now being solved, but only by about 3-4 times in this case. Since we now have to compute only a single band, however, we arrive at an answer much more quickly than before. The result, of course, is again the defect band, with a frequency of 0.378166.

Tuning the Point-Defect Mode
----------------------------

As another example utilizing the programming capability of Scheme, we will write a script to "tune" the defect mode to a particular frequency. Instead of forming a defect by simply removing a rod, we can decrease the radius or the dielectric constant of the defect rod, thereby changing the corresponding mode frequency. In this case, we'll vary the dielectric constant, and try to find a mode with a frequency of, say, 0.314159 (a random number).

We could write a loop to search for this epsilon, but instead we'll use a root-finding function provided by libctl, `(find-root function tolerance arg-min arg-max)`, that will solve the problem for us using a quadratically-convergent algorithm (Ridder's method). First, we need to define a function that takes an epsilon for the center rod and returns the mode frequency minus 0.314159; this is the function we'll be finding the root of:

```
(define old-geometry geometry) ; save the 5x5 grid with a missing rod
(define (rootfun eps)
  ; add the cylinder of epsilon = eps to the old geometry:
  (set! geometry (append old-geometry
                        (list (make cylinder (center 0 0 0)
                                    (radius 0.2) (height infinity)
                                    (material (make dielectric
                                                (epsilon eps)))))))
  (run-tm)  ; solve for the mode (using the targeted solver)
  (print "epsilon = " eps " gives freq. = " (list-ref freqs 0) "\n")
  (- (list-ref freqs 0) 0.314159))  ; return 1st band freq. - 0.314159
```

Now, we can solve for epsilon, searching in the range 1 to 12, with a fractional tolerance of 0.01, by:

```
(define rooteps (find-root rootfun 0.01 1 12))
(print "root (value of epsilon) is at: " rooteps "\n")
```

The sequence of dielectric constants that it tries, along with the corresponding mode frequencies, is:
<center>

| epsilon          | frequency         |
|------------------|-------------------|
| 1                | 0.378165893321125 |
| 12               | 0.283987088221692 |
| 6.5              | 0.302998920718043 |
| 5.14623274327171 | 0.317371748739314 |
| 5.82311637163586 | 0.309702408341706 |
| 5.41898003340128 | 0.314169110036439 |
| 5.62104820251857 | 0.311893530112625 |

</center>

The final answer that it returns is an epsilon of 5.41986120170136. Interestingly enough, the algorithm doesn't actually evaluate the function at the final point. You have to do so yourself if you want to find out how close it is to the root. Ridder's method successively reduces the interval bracketing the root by alternating bisection and interpolation steps. At the end, it does one last interpolation to give you its best guess for the root location within the current interval. If we go ahead and evaluate the band frequency at this dielectric constant, calling `(rootfun rooteps)`, we find that it is 0.314159008193209, matching our desired frequency to nearly eight decimal places after seven function evaluations! Of course, the computation isn't really this accurate anyway, due to the finite discretization.

A slight improvement can be made to the calculation above. Ordinarily, each time you call the `(run-tm)` function, the fields are initialized to random values. It would speed convergence somewhat to use the fields of the previous calculation as the starting point for the next calculation. We can do this by instead calling a lower-level function, `(run-parity TM false)`. The first parameter is the polarization to solve for, and the second tells it not to reset the fields if possible.

Editors and ctl
-------------

It is useful to have [emacs](https://en.wikipedia.org/wiki/Emacs) use its `scheme-mode` for editing ctl files, so that hitting tab indents nicely, and so on. `emacs` does this automatically for files ending with ".scm"; to do it for files ending with ".ctl" as well, add the following lines to your `~/.emacs` file:

```
 (push '("\\.ctl\\'" . scheme-mode) auto-mode-alist)
```

or if your `emacs` version is 24.3 or earlier and you have other ".ctl" files which are not Scheme:

```
 (if (assoc "\\.ctl" auto-mode-alist)
       nil
         (add-to-list 'auto-mode-alist '("\\.ctl\\'" . scheme-mode))))
```

Incidentally, `emacs` scripts are written in "elisp," a language closely related to Scheme.

If you don't use emacs (or derivatives such as Aquamacs), it would be good to find another editor that supports a Scheme mode.  For example, [jEdit](http://www.jedit.org) is a free/open-source cross-platform editor with Scheme-syntax support.  Another option is [gedit](http://projects.gnome.org/gedit/). One of our users has donated a [hilighting mode for Meep/MPB](http://github.com/hessammehr/meepmpb-highlight) that specially highlights the Meep/MPB keywords.