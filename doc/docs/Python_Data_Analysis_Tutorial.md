---
Data Analysis Tutorial
---

In the [Tutorial](Python_Tutorial.md), we focused on how to perform a calculation. Now, we'll provide a brief description on what you might do with the results of the calculations, and in particular how you might visualize the results. We'll focus on two systems: 2d and 3d.

[TOC]

Triangular Lattice of Rods
--------------------------

First, we'll return to the two-dimensional [triangular lattice of rods](Python_Tutorial.md#bands-of-a-triangular-lattice) in air. See also [our online textbook](http://ab-initio.mit.edu/book), ch. 5. The control file for this calculation, which can also be found in `meep/python/examples/mpb_tri_rods.py`, will consist of:

### The mpb\_tri\_rods.py Control File

```py
import math
import meep as mp
from meep import mpb

num_bands = 8
resolution = 32
geometry_lattice = mp.Lattice(size=mp.Vector3(1, 1),
                              basis1=mp.Vector3(math.sqrt(3) / 2, 0.5),
                              basis2=mp.Vector3(math.sqrt(3) / 2, -0.5))
geometry = [mp.Cylinder(0.2, material=mp.Medium(epsilon=12))]
k_points = [
    mp.Vector3(),               # Gamma
    mp.Vector3(y=0.5),          # M
    mp.Vector3(1 / -3, 1 / 3),  # K
    mp.Vector3(),               # Gamma
]
k_points = mp.interpolate(4, k_points)

ms = mpb.ModeSolver(
    geometry=geometry,
    geometry_lattice=geometry_lattice,
    k_points=k_points,
    resolution=resolution,
    num_bands=num_bands
)
ms.run_tm(mpb.output_at_kpoint(mp.Vector3(1 / -3, 1 / 3), mpb.fix_efield_phase,
          mpb.output_efield_z))
ms.run_te()
```

Notice that we're computing both TM and TE bands where we expect a gap in the TM bands, and are outputting the z component of the electric field for the TM bands at the K point. The `fix_efield_phase` will be explained below. The file `meep/python/examples/mpb_data_analysis.py` imports the `tri-rods` file and does the data analysis explained below.

Now, run the calculation, directing the output to a file, by entering the following command at the Unix prompt:

```
unix% python mpb_tri_rods.py >& tri-rods.out
``` 

It should finish in a few seconds.

### The tri-rods Dielectric Function

In most cases, the first thing we'll want to do is to look at the dielectric function, to make sure that we specified the correct geometry. We can do this by looking at the `epsilon.h5` output file.

The first thing that might come to mind would be to examine `epsilon.h5` directly, say by converting it to a PNG image with `h5topng` from [h5utils](https://github.com/stevengj/h5utils), magnifying it by 3:

```
unix% h5topng -S 3 epsilon.h5
```

<center>![](images/tri-rods-eps-1.gif)</center>

This image initially seems wrong! Why is the rod oval-shaped and not circular? Actually, the dielectric function is correct, but the image is distorted because the primitive cell of our lattice is a rhombus with 60<sup>$\circ$</sup> acute angles. Since the output grid of MPB is defined over the non-orthogonal unit cell, while the image produced by `h5topng` and most other plotting programs is square, the image is skewed.

We can fix the image in a variety of ways, but the best way is probably to use the `MPBData` class included in the `mpb` module. `MPBData` allows us to rearrange the data into a rectangular cell (`rectify=True`) with the same area/volume, expand the data to include multiple periods (`periods=<num_periods>`), and change the resolution per unit distance in each direction to a fixed value (`resolution=<desired_resolution>`). Use the Python `help` utility for more options (`help(mpb.MPBData)`). In this case, we'll create an `MPBData` instance that will rectify the cell, expand it to three periods in each direction, and fix the resolution to 32 pixels per *a*. Next, we call `MPBData.convert`, passing it the epsilon data we wish to convert. The whole process looks like this:

```
md = mpb.MPBData(rectify=True, periods=3, resolution=32)
eps = ms.get_epsilon()
converted_eps = md.convert(eps)
```

It's important to use `resolution` when you use `rectify`, as otherwise the non-square unit cell output by `rectify` will have a different density of grid points in each direction, and appear distorted. The output of `MPBData.convert` is by default a numpy array. We can examine it in `matplotlib`:

```py
plt.imshow(converted_eps.T, interpolation='spline36', cmap='binary')
plt.axis('off')
plt.show()
```

<center>![](images/tri-rods-eps-2.gif)</center>

As you can see, the rods are now circular as desired, and they clearly form a triangular lattice.

### Gaps and Band Diagram for tri-rods

At this point, let's check for band gaps by picking out lines with the word "Gap" in them:

```
unix% grep Gap tri-rods.out
Gap from band 1 (0.275065617068082) to band 2 (0.446289918847647), 47.4729292989213%
Gap from band 3 (0.563582903703468) to band 4 (0.593059066215511), 5.0968516236891%
Gap from band 4 (0.791161222813268) to band 5 (0.792042731370125), 0.111357548663006%
Gap from band 5 (0.838730315053238) to band 6 (0.840305955160638), 0.187683867865441%
Gap from band 6 (0.869285340346465) to band 7 (0.873496724070656), 0.483294361375001%
Gap from band 4 (0.821658212109559) to band 5 (0.864454087942874), 5.07627823271133%
```

The first five gaps are for the TM bands (which we ran first), and the last gap is for the TE bands. Note, however that the &lt; 1% gaps are probably false positives due to band crossings, as described in the [tutorial](Python_Tutorial.md#our-first-band-structure). There are no complete (overlapping TE/TM) gaps, and the largest gap is the 47% TM gap as expected. See [our online textbook](http://ab-initio.mit.edu/book), appendix C. To be absolutely sure of this and other band gaps, we would also check k-points within the interior of the Brillouin zone, but we'll omit that step here.

Next, let's plot out the band structure. To do this, we'll first extract the TM and TE bands as comma-delimited text, which can then be imported and plotted in our favorite spreadsheet/plotting program.

```
unix% grep tmfreqs tri-rods.out > tri-rods.tm.dat
unix% grep tefreqs tri-rods.out > tri-rods.te.dat
```

The TM and TE bands are both plotted below against the "k index" column of the data, with the special k-points labelled. TM bands are shown in blue (filled circles) with the gaps shaded light blue, while TE bands are shown in red (hollow circles) with the gaps shaded light red.

<center>![](images/tri-rods-bands.gif)</center>

Note that we truncated the upper frequencies at a cutoff of 1.0 c/a. Although some of our bands go above that frequency, we didn't compute enough bands to fill in all of the states in that range. Besides, we only really care about the states around the gap(s), in most cases.

### The Source of the TM Gap: Examining the Modes

Now, let's actually examine the electric-field distributions for some of the bands which were saved at the K point, remember. Besides looking neat, the field patterns will tell us about the characters of the modes and provide some hints regarding the origin of the band gap.

As before, we'll run `MPBData.convert` on the field output arrays, and view the results in `matplotlib`. We'll need to define a band function to get the efield at each band, then we call `convert` on the z component of each array:

```py
def get_efields(ms, band):
    efields.append(ms.get_efield(band, bloch_phase=True))

ms.run_tm(mpb.output_at_kpoint(mp.Vector3(1 / -3, 1 / 3), mpb.fix_efield_phase,
          get_efields))

# Create an MPBData instance to transform the efields
md = mpb.MPBData(rectify=True, resolution=32, periods=3)

converted = []
for f in efields:
    # Get just the z component of the efields
    f = f[:, :, 2]
    converted.append(md.convert(f))

for i, f in enumerate(converted):
    plt.subplot(331 + i)
    plt.contour(converted_eps.T, cmap='binary')
    plt.imshow(np.real(f).T, interpolation='spline36', cmap='RdBu', alpha=0.9)
    plt.axis('off')

plt.show()
```

Here, we've used the `bloch_phase=True` keyword argument to `get_efield` to multiply the result by its Bloch phase (see [Stored Fields and Bloch Phases]('Python_user_interface.md#stored-fields-and-bloch-phases')).

Note that we only plot the real part of the z component. Since these are TM fields, the z component is the only non-zero part of the electric field. The real and imaginary parts of the fields correspond to what the fields look like at half-period intervals in time, and in general they are different. However, at K they are redundant, due to the inversion symmetry of that k-point (proof left as an exercise for the reader). Usually, looking at the real parts alone gives you a pretty good picture of the state, especially if you use `fix_efield_phase` (see below), which chooses the phase to maximize the field energy in the real part. Sometimes, though, you have to be careful: if the real part happens to be zero, what you'll see is essentially numerical noise and you should switch to the imaginary part.

The resulting field images are shown below:

| TM band 1                                                    | TM band 2                                                    | TM band 3                                                    | TM band 4                                                    | TM band 5                                                    | TM band 6                                                    | TM band 7                                                    | TM band 8                                                    |
|--------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|--------------------------------------------------------------|
| ![center](images/tri-rods-ez1.gif) | ![center](images/tri-rods-ez2.gif) | ![center](images/tri-rods-ez3.gif) | ![center](images/tri-rods-ez4.gif) | ![center](images/tri-rods-ez5.gif) | ![center](images/tri-rods-ez6.gif) | ![center](images/tri-rods-ez7.gif) | ![center](images/tri-rods-ez8.gif) |

Your images should look the same as the ones above. If we hadn't included `fix_efield_phase` before `output_efield_z` in the .py file, on the other hand, yours would have differed slightly (e.g. by a sign or a lattice shift), because by default the phase is *random*.

When we look at the real parts of the fields, we are really looking at the fields of the modes at a particular instant in time and the imaginary part is half a period later. The point in time relative to the periodic oscillation of the state is determined by the phase of the eigenstate. The `fix_efield_phase` band function picks a canonical phase for the eigenstate, giving us a deterministic picture.

We can see several things from these plots:

First, the origin of the band gap is apparent. The lowest band is concentrated within the dielectric rods in order to minimize its frequency. The next bands, in order to be orthogonal, are forced to have a node within the rods, imposing a large "kinetic energy" (and/or "potential energy") cost and hence a gap. See [our online textbook](http://ab-initio.mit.edu/book), ch. 5. Successive bands have more and more complex nodal structures in order to maintain orthogonality. The contrasting absence of a large TE gap has to do with boundary conditions. The perpendicular component of the displacement field must be continuous across the dielectric boundary, but the parallel component need not be.

We can also see the deep impact of symmetry on the states. The K point has C<sub>3v</sub> symmetry and not quite the full C<sub>6v</sub> symmetry of the dielectric structure. This symmetry group has only one two-dimensional representation--that is what gives rise to the degenerate pairs of states (2/3, 4/5, and 7/8), all of which fall into this "p-like" category where the states transform like two orthogonal dipole field patterns, essentially. The other two bands, 1 and 6, transform under the trivial "s-like" representation with band 6 just a higher-order version of 1.

Diamond Lattice of Spheres
--------------------------

Now, let us turn to a three-dimensional structure, a diamond lattice of dielectric spheres in air. See [our online textbook](http://ab-initio.mit.edu/book), ch. 6. The basic techniques to compute and analyze the modes of this structure are the same as in two dimensions, but of course, everything becomes more complicated in 3d. It's harder to find a structure with a complete gap, the modes are no longer polarized, the computations are bigger, and visualization is more difficult, for starters.

The control file for this calculation, which can also be found in `meep/python/examples/mpb_diamond.py`, consists of:

### Diamond Control File

```py
import math
import meep as mp
from meep import mpb

# Dielectric spheres in a diamond (fcc) lattice.  This file is used in
# the "Data Analysis Tutorial" section of the MPB manual.

sqrt_half = math.sqrt(0.5)
geometry_lattice = mp.Lattice(
    basis_size=mp.Vector3(sqrt_half, sqrt_half, sqrt_half),
    basis1=mp.Vector3(0, 1, 1),
    basis2=mp.Vector3(1, 0, 1),
    basis3=mp.Vector3(1, 1)
)

# Corners of the irreducible Brillouin zone for the fcc lattice,
# in a canonical order:
vlist = [
    mp.Vector3(0, 0.5, 0.5),        # X
    mp.Vector3(0, 0.625, 0.375),    # U
    mp.Vector3(0, 0.5, 0),          # L
    mp.Vector3(0, 0, 0),            # Gamma
    mp.Vector3(0, 0.5, 0.5),        # X
    mp.Vector3(0.25, 0.75, 0.5),    # W
    mp.Vector3(0.375, 0.75, 0.375)  # K
]

k_points = mp.interpolate(4, vlist)

# define a couple of parameters (which we can set from the command_line)
eps = 11.56  # the dielectric constant of the spheres
r = 0.25  # the radius of the spheres

diel = mp.Medium(epsilon=eps)

# A diamond lattice has two "atoms" per unit cell:
geometry = [mp.Sphere(r, center=mp.Vector3(0.125, 0.125, 0.125), material=diel),
            mp.Sphere(r, center=mp.Vector3(-0.125, -0.125, -0.125), material=diel)]

# A simple fcc lattice would have only one sphere/object at the origin.

resolution = 16  # use a 16x16x16 grid
mesh_size = 5
num_bands = 5

ms = mpb.ModeSolver(
    geometry_lattice=geometry_lattice,
    k_points=k_points,
    geometry=geometry,
    resolution=resolution,
    num_bands=num_bands,
    mesh_size=mesh_size
)


# run calculation, outputting electric_field energy density at the U point:
ms.run(mpb.output_at_kpoint(mp.Vector3(0, 0.625, 0.375), mpb.output_dpwr))
```

As before, run the calculation, directing the output to a file. This will take a minute. We'll put it in the background with `nohup` so that it will finish even if we log out:

```
unix% nohup python mpb_diamond.py >& diamond.out &
```

### Important Note on Units for the Diamond/FCC Lattice

[As usual](Python_Tutorial.md#a-few-words-on-units), all distances are in the dimensionless units determined by the length of the lattice vectors. We refer to these units as *a*, and frequencies are given in units of *c/a*. By default, the lattice/basis vectors are unit vectors, but in the case of fcc lattices this conflicts with the convention in the literature. In particular, the canonical *a* for fcc is the edge-length of a cubic supercell containing the lattice.

In order to follow this convention, we set the length of our basis vectors appropriately using the `basis_size` property of the `Lattice` class. The lattice vectors default to the same length as the basis vectors. If the cubic supercell edge has unit length (*a*), then the fcc lattice vectors have length sqrt(0.5).

### Gaps and Band Diagram for the Diamond Lattice

The diamond lattice has a complete band gap:

```
unix% grep Gap diamond.out
Gap from band 2 (0.396348703007373) to band 3 (0.440813418580596), 10.6227251392791%
```

We can also plot its band diagram, much as for the tri-rods case except that now we can't classify the bands by polarization.

```
unix% grep freqs diamond.out > diamond.dat
```

The resulting band diagram, with the complete band gap shaded yellow, is shown below. Note that we only computed 5 bands, so in reality the upper portion of the plot would contain a lot more bands which are of less interest than the bands adjoining the gap.

<center>![](images/diamond-bands.gif)</center>

### Visualizing the Diamond Lattice Structure and Bands

TODO
