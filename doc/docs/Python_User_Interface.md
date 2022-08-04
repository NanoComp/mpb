---
# Python User Interface
---

The Python user interface is documented in this page. It currently depends on and lives in the [Meep repository](https://github.com/NanoComp/meep), but will eventually be migrated to the MPB repo. Installing Meep with Python enabled will automatically build PyMPB if MPB is installed. Detailed instructions can be found [here](http://meep.readthedocs.io/en/latest/Installation/#building-from-source). We do not document the Python language or the functions provided by [meep](https://meep.readthedocs.io). See also the [PyMeep User Interface](https://meep.readthedocs.io/en/latest/Python_User_Interface/) section of the [Meep manual](https://meep.readthedocs.io).

[TOC]

The `ModeSolver` Class
---------------------------

The `ModeSolver` class stores all the variables that you can set to control various parameters of the MPB computation. They are also listed, along with their default values, by the `help` command in a Python interpreter.

```python
from meep.mpb import ModeSolver
help(ModeSolver)
```

Here is the constructor's function signature.

```python
ModeSolver(object):

    def __init__(self,
                 resolution=10,
                 is_negative_epsilon_ok=False,
                 eigensolver_flops=0,
                 is_eigensolver_davidson=False,
                 eigensolver_nwork=3,
                 eigensolver_block_size=-11,
                 eigensolver_flags=68,
                 use_simple_preconditioner=False,
                 force_mu=False,
                 mu_input_file='',
                 epsilon_input_file='',
                 mesh_size=3,
                 target_freq=0.0,
                 tolerance=1.0e-7,
                 num_bands=1,
                 k_points=[],
                 ensure_periodicity=True,
                 geometry=[],
                 geometry_lattice=mp.Lattice(),
                 geometry_center=mp.Vector3(0, 0, 0),
                 default_material=mp.Medium(epsilon=1),
                 dimensions=3,
                 random_fields=False,
                 filename_prefix='',
                 deterministic=False,
                 verbose=False):
```

In brackets after each variable is the type of value that it should hold. The classes, complex datatypes like `GeometricObject`, are described in a later subsection. The basic datatypes, like `bool`, `complex`, and `float`, are Python builtins.

**`geometry` [ list of `GeometricObject` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the geometric objects making up the structure being simulated. When objects overlap, later objects in the list take precedence. Defaults to no objects (empty list).

**`default_material` [ `Medium` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Holds the default material that is used for points not in any object of the geometry list. Defaults to air (epsilon of 1). See also `epsilon_input_file`, below.

**`ensure_periodicity` [`bool`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If true (the default), then geometric objects are treated as if they were shifted by all possible lattice vectors; i.e. they are made periodic in the lattice.

**`geometry_lattice` [`Lattice` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the basis vectors and lattice size of the computational cell which is centered on the origin of the coordinate system. These vectors form the basis for all other 3-vectors in the geometry, and the lattice size determines the size of the primitive cell. If any dimension of the lattice size is 0, then the dimension of the lattice is reduced (i.e. it becomes two- or one-dimensional). That is, the dielectric function becomes two-dimensional; it is still, in principle, a three dimensional system, and the k-point vectors can be three-dimensional. Generally, you should make any 0  dimension(s) perpendicular to the others. Defaults to the orthogonal x-y-z vectors of unit length (i.e. a square/cubic lattice).

**`resolution` [`integer` or `Vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies the computational grid resolution, in pixels per lattice unit (a lattice unit is one basis vector in a given direction). If `resolution` is a `Vector3`, then specifies a different resolution for each direction; otherwise the resolution is uniform. The grid size is then the product of the lattice size and the resolution, rounded up to the next positive integer. Defaults to `10`. You can call `ModeSolver.optimize_grid_size()` *after* setting the `resolution` and `geometry_lattice` to adjust the grid size for maximal performance. This rounds the grid size in each direction to the nearest integer with small factors, to improve FFT speed.

**`dimensions` [`integer`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Explicitly specifies the dimensionality of the simulation; if the value is less than 3, the sizes of the extra dimensions in `grid_size` are ignored (assumed to be one). Defaults to 3. *Deprecated:* the preferred method is to set `geometry_lattice` to have size 0 in any unwanted dimensions.

**`k_points` [ list of `Vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
List of Bloch wavevectors to compute the bands at, expressed in the basis of the reciprocal lattice vectors. The reciprocal lattice vectors are defined as follows (see [our online textbook](http://ab-initio.mit.edu/book), appendix B): Given the lattice vectors R<sub>i</sub> (*not* the basis vectors), the reciprocal lattice vector G<sub>j</sub> satisfies R<sub>i</sub> * G<sub>j</sub> = 2*&#960;*&#948;<sub>i,j</sub>, where &#948;<sub>i,j</sub> is the Kronecker delta (1 for `i=j` and 0 otherwise). R<sub>i</sub> for any `0` dimensions is taken to be the corresponding basis vector. Normally, the wavevectors should be in the first Brillouin zone ([see below](Python_User_Interface.md#coordinate-conversion-functions)). `k_points` defaults to none (empty list).

**`num_bands` [`integer`]**   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Number of bands (eigenvectors) to compute at each k point. Defaults to 1.

**`target_freq` [`float`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If zero, the lowest-frequency `num_bands` states are solved for at each k point (ordinary eigenproblem). If non-zero, solve for the `num_bands` states whose frequencies have the smallest absolute difference with `target_freq` (special, "targeted" eigenproblem). Beware that the targeted solver converges more slowly than the ordinary eigensolver and may require a lower `tolerance` to get reliable results. Defaults to 0.

**`tolerance` [`float`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Specifies when convergence of the eigensolver is judged to have been reached when the eigenvalues have a fractional change less than `tolerance` between iterations. Defaults to 1.0e-7.

**`filename_prefix` [`string`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A string prepended to all output filenames. Defaults to `"FILE-"`, where your script name is FILE.py. You can change this to the empty string `""` to use no prefix.

**`epsilon_input_file` [`string`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
If this string is not `""` (the default), then it should be the name of an HDF5 file whose first/only dataset defines a dielectric function over some discrete grid. This dielectric function is then used in place of `default_material` (*i.e.* where there are no `geometry` objects). The grid of the epsilon file dataset need not match `grid_size`; it is scaled and/or linearly interpolated as needed. The lattice vectors for the epsilon file are assumed to be the same as `geometry_lattice`. Note that, even if the grid sizes match and there are no geometric objects, the dielectric function used by MPB will not be exactly the dielectric function of the epsilon file, unless you also set `mesh_size` to 1 (see above).

**`eigensolver_block_size` [`integer`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The eigensolver uses a "block" algorithm, which means that it solves for several bands simultaneously at each k-point. `eigensolver_block_size` specifies this number of bands to solve for at a time; if it is zero or &gt;= `num_bands`, then all the bands are solved for at once. If `eigensolver_block_size` is a negative number, -*n*, then MPB will try to use nearly-equal block-sizes close to *n*. Making the block size a small number can reduce the memory requirements of MPB, but block sizes &gt; 1 are usually more efficient. There is typically some optimum size for any given problem. Defaults to -11 (i.e. solve for around 11 bands at a time).

**`use_simple_preconditioner` [`bool`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Whether or not to use a simplified preconditioner. Defaults to `False` which is fastest most of the time. Turning this on increases the number of iterations, but decreases the time for each iteration.

**`deterministic` [`bool`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Since the fields are initialized to random values at the start of each run, there are normally slight differences in the number of iterations, etcetera, between runs. Setting `deterministic` to `True` makes things deterministic. The default is `False`.

**`eigensolver_flags` [`integer`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
This variable is undocumented and reserved for use by Jedi Masters only.

**`mesh_size` [`integer`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The mesh size determines the window size over which sub-pixel smoothening happens. Setting the `mesh_size` to `1` disables sub-pixel smoothing.

Predefined Variables
--------------------

Variables predefined for your convenience and amusement, available via the `meep` package.

**`air`, `vacuum` [`Medium` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Two aliases for a predefined material type with a dielectric constant of 1.

**`nothing` [`material-type` class ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A material that, effectively, punches a hole through other objects to the background (`default_material` or `epsilon_input_file`).

**`inf` [`float`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A big number (1e20) to use for "infinite" dimensions of objects.

Output Variables
----------------

Attributes of `ModeSolver` whose values are set upon completion of the eigensolver.

**`freqs` [ list of `float`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A list of the frequencies of each band computed for the last k point. Guaranteed to be sorted in increasing order. The frequency of band `b` can be retrieved via `ms.freqs[b - 1]`.

**`iterations` [`integer`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The number of iterations required for convergence of the last k point.

**`parity` [`string`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A string describing the current required parity/polarization (`te`, `zeven`, etcetera, or "" for none). Useful for prefixing output lines for grepping.

Yet more global variables are set by the `run` function and its variants, for use after `run` completes or by a band function which is called for each band during the execution of `run`.

**`current_k` [`Vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The k point most recently solved from the `k_points` list.

**`gap_list` [ list of (*`percent, freq_min, freq_max`*) tuples ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
This is a list of the gaps found by the eigensolver, and is set by the `run` functions when two or more k-points are solved. It is the empty list if no gaps are found.

**`band_range_data` [ list of ((*`min, kpoint`*), (*`max, kpoint`*)) ]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
For each band, this list contains the minimum and maximum frequencies of the band, and the associated k points where the extrema are achieved. Note that the bands are defined by sorting the frequencies in increasing order, so this can be confused if two bands cross.

Classes
-------

Meep defines several types of classes that are useful to MPB, the most numerous of which are the various geometric object classes. 

### Lattice

The lattice class is normally used only for the `geometry_lattice` variable and specifies the three lattice directions of the crystal and the lengths of the corresponding lattice vectors.

**`basis1, basis2, basis3` [`Vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The three lattice directions of the crystal, specified in the cartesian basis. The lengths of these vectors are ignored--only their directions matter. The lengths are determined by the `basis_size` property, below. These vectors are then used as a basis for all other 3-vectors in the script. They default to the x, y, and z directions, respectively.

**`basis_size` [`Vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The components of `basis_size` are the lengths of the three basis vectors, respectively. They default to unit lengths.

**`size` [`Vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The size of the lattice (i.e. the length of the lattice vectors R<sub>i</sub>, in which the crystal is periodic) in units of the basis vectors. Thus, the actual lengths of the lattice vectors are given by the components of `size` multiplied by the components of `basis_size`. Alternatively, you can think of `size` as the vector between opposite corners of the primitive cell, specified in the lattice basis. Defaults to unit lengths.

If any dimension has a size of `0`, then the dimensionality of the problem is reduced by one. Strictly speaking, the dielectric function is taken to be uniform along that dimension. In this case, the `0` dimension should generally be orthogonal to the other dimensions.

### Medium

This class is used to specify the materials that geometric objects are made of. See the [Meep manual](http://meep.readthedocs.io/en/latest/Python_User_Interface/#medium) for more information.

**`epsilon` [`float`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The dielectric constant (must be positive). Default is 1. You can also use `index=n` as a synonym for `epsilon=n*n`.

**`epsilon-diag` [`Vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The diagonal elements (a b c) of the dielectric tensor. Defaults to (1, 1, 1).

**`epsilon_offdiag` [`Vector3`]**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
The off-diagonal elements (u v w) of the dielectric tensor. Defaults to zero. This is a `Vector3` for which the components may be complex numbers (e.g. `3+0.1j`). If non-zero imaginary parts are specified, then the dielectric tensor is complex-hermitian. This is only supported when MPB is configured with the `--with-hermitian-eps` flag. This is not dissipative (the eigenvalues of epsilon are real), but rather breaks time-reversal symmetry, corresponding to a gyrotropic (magneto-optic) material (see [our online textbook](http://ab-initio.mit.edu/book), ch. 2). Note that [inversion symmetry](Python_User_Interface.md#inversion-symmetry) may not mean what you expect for complex-hermitian epsilon, so be cautious about using `mpbi` in this case.

For example, a material with a dielectric constant of 3.0 for P-polarization and 5.0 for S-polarization would be specified via `m = mp.Medium(epsilon_diag=mp.Vector3(3, 3, 5))`. Please [be aware](Python_User_Interface.md#run-functions) that not all 2d anisotropic dielectric structures will have P- and S-polarized modes, however.

**material functions**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
This material type allows you to specify the material as an arbitrary function of position. For an example of this, see the `mpb_bragg_sine.py` file in the `python/examples/` directory of the [Meep source](http://github.com/NanoComp/meep). See the documentation in the [Meep manual](http://meep.readthedocs.io/en/latest/Python_User_Interface/#medium).

Normally, the dielectric constant is required to be positive or positive-definite, for a tensor. However, MPB does have a somewhat experimental feature allowing negative dielectrics (e.g. in a plasma). To use it, call the function `ModeSolver.allow_negative_epsilon()` before `ModeSolver.run()`. In this case, it will output the (real) frequency *squared* in place of the (possibly imaginary) frequencies. Convergence will be somewhat slower because the eigenoperator is not positive definite.

### GeometricObject

This class, and its descendants, are used to specify the solid geometric objects that form the dielectric structure being simulated. You can read about them in the [Meep manual](http://meep.readthedocs.io/en/latest/Python_User_Interface/#geometricobject).

Here are some examples of geometric objects created using the above classes, assuming that the lattice directions (the basis) are just the ordinary unit axes, and `m` is some material we have defined:

```
import meep as mp

# A cylinder of infinite radius and height 0.25 pointing along the x axis,
# centered at the origin:
c = mp.Cylinder(radius=mp.inf, material=m, height=0.25, axis=mp.Vector3(x=1))
```

```
# An ellipsoid with its long axis pointing along (1,1,1), centered on
# the origin (the other two axes are orthogonal and have equal
# semi-axis lengths):
e = mp.Ellipsoid(material=m, size=mp.Vector3(0.8, 0.2, 0.2), e1=mp.Vector3(1, 1, 1),
                 e2=mp.Vector3(0, 1, -1), e3=mp.Vector3(-2, 1, 1))
```

```
# A unit cube of material m with a spherical air hole of radius 0.2 at
# its center, the whole thing centered at (1,2,3):
from meep import mpb
ms = mpb.ModeSolver()
ms.geometry = [mp.Block(center=mp.Vector3(1, 2, 3), material=m, size=mp.Vector3(1, 1, 1)),
               mp.Sphere(center=mp.Vector3(1, 2, 3), material=mp.air, radius=0.2)]
```

### MPBArray

Most functions in MPB that return numpy arrays actually return `MPBArray`s, which is just a wrapper around a numpy array that stores some additional data:

+ `lattice`: This attribute stores the `ModeSolver`'s geometry_lattice as a convenient way to pass it to `MPBData`. See the [Python Data Analysis Tutorial](Python_Data_Analysis_Tutorial.md) for more information.

+ `kpoint`: The current kpoint at the time this array was created. This is also a convenient way to pass this information to `MPBData`.

+ `bloch_phase`: `True` if the array has been multiplied by the Bloch phase, `False` otherwise.

Functions
---------

Here, we describe the functions that are defined by MPB. There are many types of functions defined ranging from utility functions for duplicating geometric objects to run functions that start the computation.

### Geometry Utilities

Some utility functions are provided to help you manipulate geometric objects. These are available in the `meep` package:

**`geometric_object_duplicates(shift_vector, min_multiple, max_multiple, go)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return a list of duplicates of `go`, shifted by various multiples of `shift_vector` from `min_multiple` to `max_multiple`, inclusive, in steps of 1.

**`geometric_objects_duplicates(shift_vector, min_multiple, max_multiple, go_list)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Same as `geometric_object_duplicates`, except operates on a list of objects, `go_list`. If *A* appears before *B* in the input list, then all the duplicates of *A* appear before all the duplicates of *B* in the output list.

**`geometric_objects_lattice_duplicates(lat, go_list, [ ux, uy, uz ])`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Duplicates the objects in `obj_list` by multiples of the lattice basis vectors in `lat`, making all possible shifts of the "primitive cell" (see below) that fit inside the lattice cell. This is useful for supercell calculations. See the [Tutorial](Python_Tutorial.md). The primitive cell to duplicate is `ux` by `uy` by `uz`, in units of the basis vectors. These three parameters are optional; any that you do not specify are assumed to be `1`.

**`is_point_in_object(point, obj)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns whether or not the given 3-vector `point` is inside the geometric object `obj`.

**`is_point_in_periodic_object(point, obj)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As `point_in_object?`, but also checks translations of the given object by the lattice vectors.

**`display_geometric_object_info(indent_by, obj)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Outputs some information about the given `obj`, indented by `indent_by` spaces.

### Coordinate Conversion Functions

The following functions allow you to easily convert back and forth between the lattice, cartesian, and reciprocal bases. See also the [note on units](Python_Tutorial.md#a-few-words-on-units) in the tutorial.

**`lattice_to_cartesian(x, lat)`, `cartesian_to_lattice(x, lat)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Convert `x` between the lattice basis of `lat`, (the basis of the lattice vectors normalized to `lat.basis_size`) and the ordinary cartesian basis, where `x` is either a `Vector3` or a `Matrix`, returning the transformed vector/matrix. In the case of a matrix argument, the matrix is treated as an operator on vectors in the given basis, and is transformed into the same operator on vectors in the new basis.

**`reciprocal_to_cartesian(x, lat)`, `cartesian_to_reciprocal(x, lat)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Like the above, except that they convert to/from reciprocal space (the basis of the reciprocal lattice vectors). Also, the cartesian vectors output/input are in units of 2&#960;.

**`reciprocal_to_lattice(x, lattice)`, `lattice_to_reciprocal(x, lat)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Convert between the reciprocal and lattice bases, where the conversion again leaves out the factor of 2&#960; (i.e. the lattice-basis vectors are assumed to be in units of 2&#960;).

Also, a couple of rotation functions are defined on the `Vector3` class, for convenience, so that you don't have to explicitly convert to cartesian coordinates in order to use the `Vector3.rotate` function.

**`Vector3.rotate_lattice(axis, theta, lat)`, `Vector3.rotate_reciprocal(axis, theta, lat)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Like `Vector3.rotate`, except that `axis` and the `Vector3` instance are specified in the lattice/reciprocal bases of `lat`.

Usually, k-points are specified in the first Brillouin zone, but sometimes it is convenient to specify an arbitrary k-point. However, the accuracy of MPB degrades as you move farther from the first Brillouin zone due to the choice of a fixed planewave set for a basis. This is easily fixed: simply transform the k-point to a corresponding point in the first Brillouin zone, and a completely equivalent solution (identical frequency, fields, etcetera) is obtained with maximum accuracy. The following method of the `ModeSolver` class accomplishes this:

**`first_brillouin_zone(k)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a k-point *`k`* in the basis of the reciprocal lattice vectors, as usual, return an equivalent point in the first Brillouin zone of the current lattice (`geometry_lattice`).

Note that `first_brillouin_zone` can be applied to the entire `k_points` list with the Python expression: `map(ms.first_brillouin_zone, ms.k_points)`.

### Run Functions

These are functions to help you run and control the simulation. They are all methods of the `ModeSolver` class. The ones you will most commonly use are the `run` function and its variants. The syntax of these functions, and one lower-level function, is:

**`ModeSolver.run(*band_func)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
This runs the simulation described by the input parameters (see above), with no constraints on the polarization of the solution. That is, it reads the input parameters, initializes the simulation, and solves for the requested eigenstates of each k-point. `run` takes as arguments zero or more "band functions" `band_func`. A band function should be a function of one integer argument, the band index, so that `band_func(which_band)` performs some operation on the band `which_band` (e.g. outputting fields). After every k-point, each band function is called for the indices of all the bands that were computed. Alternatively, a band function may be a "thunk" (function of zero arguments), in which case `band_func()` is called exactly once per k-point.

**`ModeSolver.run_zeven(*band_func)`, ``Modesolver.run_zodd(*band_func)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
These are the same as the `run` function except that they constrain their solutions to have even and odd symmetry with respect to the z=0 plane. You should use these functions *only* for structures that are symmetric through the z=0 mirror plane, where the third basis vector is in the z direction (0,0,1) and is orthogonal to the other two basis vectors, and when the k vectors are in the xy plane. Under these conditions, the eigenmodes always have either even or odd symmetry. In two dimensions, even/odd parities are equivalent to P/S polarizations, respectively and are often strongly analogous even in 3d. Such a symmetry classification is useful for structures such as waveguides and photonic-crystal slabs. See the [online book](http://ab-initio.mit.edu/book) (ch. 3).

**`ModeSolver.run_te(*band_func)`, `ModeSolver.run_tm(*band_func)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
These are the same as the `run` function except that they constrain their solutions to be TE- and TM-polarized, respectively, in two dimensions. The TE and TM polarizations are defined has having electric and magnetic fields in the xy plane, respectively. Equivalently, the H/E field of TE/TM light has only a z component (making it easier to visualize).

These functions are actually equivalent to calling `run_zeven` and `run_zodd`, respectively. Note that for the modes to be segregated into TE and TM polarizations, the dielectric function must have mirror symmetry for reflections through the xy plane. If you use [anisotropic dielectrics](Python_User_Interface.md#material_type), you should be aware that they break this symmetry if the z direction is not one of the principle axes. If you use `run_te` or `run_tm` in such a case of broken symmetry, MPB will exit with an error.

**`ModeSolver.run_yeven(*band_func)`, `ModeSolver.run_yodd(*band_func)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
These functions are analogous to `run_zeven` and `run_zodd`, except that they constrain their solutions to have even and odd symmetry with respect to the y=0 plane. You should use these functions *only* for structures that are symmetric through the y=0 mirror plane, where the second basis vector is in the y direction (0,1,0) and is orthogonal to the other two basis vectors, and when the k vectors are in the xz plane.

**`run_yeven_zeven`, `run_yeven_zodd`, `run_yodd_zeven`, `run_yodd_zodd`, `run_te_yeven`, `run_te_yodd`, `run_tm_yeven`, `run_tm_yodd`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
These `run`-like functions combine the `yeven`/`yodd` constraints with `zeven`/`zodd` or `te`/`tm`. See also `run_parity`, below.

**`ModeSolver.run_parity(p, reset_fields, *band-func)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Like the `run` function, except that it takes two extra parameters, a parity `p` and a boolean (`True`/`False`) value `reset_fields`. `p` specifies a parity constraint, and should be one of the predefined variables:

_   `mp.NO_PARITY`: equivalent to `run`
_   `mp.EVEN_Z` (or `TE`): equivalent to `run_zeven` or `run_te`
_   `mp.ODD_Z` (or `TM`): equivalent to `run_zodd` or `run_tm`
_   `mp.EVEN_Y` (like `EVEN_Z` but for y=0 plane)
_   `mp.ODD_Y` (like `ODD_Z` but for y=0 plane)

It is possible to specify more than one symmetry constraint simultaneously by adding them, e.g. `EVEN_Z + ODD_Y` requires the fields to be even through z=0 and odd through y=0. It is an error to specify incompatible constraints (e.g. `EVEN_Z + ODD_Z`. **Important:** if you specify the z/y parity, the dielectric structure *and* the k vector **must** be symmetric about the z/y=0 plane, respectively. If `reset_fields` is `False`, the fields from any previous calculation will be reused as the starting point from this calculation, if possible; otherwise, the fields are reset to random values. The ordinary `run` functions use a default `reset_fields` of`True`. Alternatively, `reset_fields` may be a string, the name of an HDF5 file to load the initial fields from as exported by `save_eigenvectors`, as shown [below](Python_User_Interface.md#manipulating-the-raw-eigenvectors).

**`ModeSolver.display_eigensolver_stats()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Display some statistics on the eigensolver convergence; this function is useful mainly for MPB developers in tuning the eigensolver.

Several band functions for outputting the eigenfields are defined for your convenience, and are described in the [section below](Python_User_Interface.md#bandoutput-functions). You can also define your own band functions, and for this purpose the functions described in the section **Field manipulation functions**, below, are useful. A band function takes the form:

```
def my_band_func(ms, which_band):
...do stuff here with band index which_band and ModeSolver instance ms...
)
```

Note that the output variable `freqs` may be used to retrieve the frequency of the band (see above). Also, a global variable `current_k` is defined holding the current k-point vector from the `k_points` list.

There are also some even lower-level functions that you can call, although you should not need to do most of the time:

**`ModeSolver.init_params(p, reset_fields)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Read the input variables and initialize the simulation in preparation for computing the eigenvalues. The parameters are the same as the first two parameters of `run_parity`. This function *must* be called before any of the other simulation functions below. Note, however, that the `run` functions all call `init_params`.

**`ModeSolver.set_parity(p)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
After calling `init_params`, you can change the parity constraint without resetting the other parameters by calling this function. Beware that this does not randomize the fields (see below); you don't want to try to solve for, say, the TM eigenstates when the fields are initialized to TE states from a previous calculation.

**`ModeSolver.randomize_fields()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Initialize the fields to random values.

**`ModeSolver.solve_kpoint(k)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Solve for the requested eigenstates at the Bloch wavevector `k`.

### The Inverse Problem: k as a Function of Frequency

MPB's `run` function(s) and its underlying algorithms compute the frequency `w` as a function of wavevector `k`. Sometimes, however, it is desirable to solve the inverse problem, for `k` at a given frequency `w`. This is useful, for example, when studying coupling in a waveguide between different bands at the same frequency since frequency is conserved even when wavevector is not. One also uses `k(w)` to construct wavevector diagrams, which aid in understanding diffraction (e.g. negative-diffraction materials and super-prisms). To solve such problems, therefore, we provide the `find_k` function described below, which inverts `w(k)` via a few iterations of Newton's method using the group velocity `dw/dk`. Because it employs a root-finding method, you need to specify bounds on `k` and a *crude* initial guess where order of magnitude is usually good enough.

**`ModeSolver.find_k(p, omega, band_min, band_max, korig_and_kdir, tol, kmag_guess, kmag_min, kmag_max, *band_func)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Find the wavevectors in the current geometry/structure for the bands from `band_min` to `band_max` at the frequency `omega` along the `kdir` direction in k-space. Returns a list of the wavevector magnitudes for each band; the actual wavevectors are `kdir.unit().scale(magnitude)`. The arguments of `find_k` are:

-   `p`: parity (same as first argument to `run_parity`, [above](Python_User_Interface.md#run-functions)).
-   `omega`: the frequency at which to find the bands
_   `band_min`, *`band_max`*: the range of bands to solve for the wavevectors of (inclusive).
-   `korig_and_kdir`: If this is a `list` of `Vector3`, the first element should be the original `k` and the second element is the direction in k-space in which to find the wavevectors. (The magnitude of *`kdir`* is ignored.) If this is a single `Vector3` it represents `kdir`, and `korig` defaults to `Vector3(0, 0, 0)`.
-   `tol`: the fractional tolerance with which to solve for the wavevector; `1e-4` is usually sufficient. (Like the `tolerance` input variable, this is only the tolerance of the numerical iteration...it does not have anything to do with e.g. the error from finite grid `resolution`.)
_   `kmag_guess`: an initial guess for the k magnitude (along *`kdir`*) of the wavevector at *`omega`*. Can either be a list (one guess for each band from *`band_min`* to *`band_max`*) or a single number (same guess for all bands, which is usually sufficient).
_   `kmag_min`, *`kmag_max`*: a range of k magnitudes to search; should be large enough to include the correct k values for all bands.
_   `band_func`: zero or more [band functions](Python_User_Interface.md#bandoutput_functions), just as in `run`, which are evaluated at the computed k points for each band.

The `find_k` routine also prints a line suitable for grepping:

```
kvals: omega, band-min, band-max, korig1, korig2, korig3, kdir1, kdir2, kdir3, k magnitudes...
```

### Band/Output Functions

All of these are functions that, given a `ModeSolver` instance and a band index, output the corresponding field or compute some function thereof in the primitive cell of the lattice. They are designed to be passed as band functions to the `run` routines, although they can also be called directly. See also the section on [field normalizations](Python_User_Interface.md#field-normalization).

The output functions are available in the `mpb` module of the `meep` package and can be called as follows:

```python
from meep import mpb
ms = mpb.ModeSolver()
mpb.output_hfield(ms, band)
```

**`output_hfield(ms, which-band)`**  
**`output_hfield_x(ms, which_band)`**  
**`output_hfield_y(ms, which_band)`**  
**`output_hfield_z(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the magnetic (\(\mathbf{H}\)) field for `which_band`; either all or one of the Cartesian components, respectively.

**`output_dfield(ms, which_band)`**  
**`output_dfield_x(ms, which_band)`**  
**`output_dfield_y(ms, which_band)`**  
**`output_dfield_z(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the electric displacement (\(\mathbf{D}\)) field for `which_band`; either all or one of the Cartesian components, respectively.

**`output_efield(ms, which_band)`**  
**`output_efield_x(ms, which_band)`**  
**`output_efield_y(ms, which_band)`**  
**`output_efield_z(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the electric (\(\mathbf{E}\)) field for `which_band`; either all or one of the Cartesian components, respectively.

**`output_bpwr(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the time-averaged magnetic-field energy density (bpwr = \(\mu|\mathbf{H}|^2\)) for `which_band`.

**`output_dpwr(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the time-averaged electric-field energy density (dpwr = \(\varepsilon|\mathbf{E}|^2\)) for `which_band`.

**`fix_hfield_phase(ms, which_band)`**  
**`fix_dfield_phase(ms, which_band)`**  
**`fix_efield_phase(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Fix the phase of the given eigenstate in a canonical way based on the given spatial field. See also `fix_field_phase`, below. Otherwise, the phase is random. These functions also maximize the real part of the given field so that one can hopefully just visualize the real part. To fix the phase for output, pass one of these functions to `run` before the corresponding output function, e.g. `ms.run_tm(mpb.fix_dfield_phase, mpb.output_dfield_z)`

Although we try to maximize the "real-ness" of the field, this has a couple of limitations. First, the phase of the different field components cannot, of course, be chosen independently, so an individual field component may still be imaginary. Second, if you use `mpbi` to take advantage of [inversion symmetry](Python_User_Interface.md#inversion-symmetry) in your problem, the phase is mostly determined elsewhere in the program; `fix_Xfield_phase` in that case only determines the sign.

See also below for the `output_poynting` and `output_tot_pwr` functions to output the Poynting vector and the total electromagnetic energy density, respectively, and the `output_charge_density` function to output the bound charge density.

Sometimes, you only want to output certain bands. For example, here is a function that, given an band/output function like the ones above, returns a new output function that only calls the first function for bands with a large fraction of their energy in an object(s). This is useful for picking out defect states in supercell calculations.

**`output_dpwr_in_objects(band_func, min_energy, *objects)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a band function `band_func`, returns a new band function that only calls `band_func` for bands having a fraction of their electric-field energy greater than `min_energy` inside the given objects (zero or more geometric objects). Also, for each band, prints the fraction of their energy in the objects in the following form which is suitable for grepping:

```
dpwr:, band-index, frequency, energy-in-objects
```

`output_dpwr_in_objects` only takes a single band function as a parameter, but if you want it to call several band functions, you can easily combine them into one with the following routine:

**`combine_band_functions(*band_funcs)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more band functions, returns a new band function that calls all of them in sequence. When passed zero parameters, returns a band function that does nothing.

It is also often useful to output the fields only at a certain k-point, to let you look at typical field patterns for a given band while avoiding gratuitous numbers of output files. This can be accomplished via:

**`output_at_kpoint(k_point, *band_funcs)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given zero or more band functions, returns a new band function that calls all of them in sequence, but only at the specified `k_point`. For other k-points, does nothing.

### Miscellaneous Functions

**`retrieve_gap(lower_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return the frequency gap from the band *`lower_band`* to the band *`lower_band+1`*, as a percentage of mid-gap frequency. The "gap" may be negative if the maximum of the lower band is higher than the minimum of the upper band. The gap is computed from the `band_range_data` of the previous run.

#### Parity

Given a set of eigenstates at a k-point, MPB can compute their *parities* with respect to the z=0 or y=0 plane. The z/y parity of a state is defined as the expectation value under the usual inner product of the mirror-flip operation through z/y=0, respectively. For true even and odd eigenstates (see e.g. `run_zeven` and `run_zodd`), this will be +1 and -1, respectively; for other states it will be something in between.

This is useful e.g. when you have a nearly symmetric structure, such as a waveguide with a substrate underneath, and you want to tell which bands are even-like (parity &gt; 0) and odd-like (parity &lt; 0). Indeed, any state can be decomposed into purely even and odd functions, with absolute-value-squared amplitudes of (1+parity)/2 and (1-parity)/2, respectively.

**`display_zparities()`, `display_yparities()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
These are band functions, designed to be passed to `run`, which output all of the z/y parities, respectively, at each k-point in comma-delimited format suitable for grepping.

**`ModeSolver.compute_zparities()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns a list of the parities about the z=0 plane, one number for each band computed at the last k-point.

**`ModeSolver.compute_yparities()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns a list of the parities about the y=0 plane, one number for each band computed at the last k-point.

Note that the magnetic field is only a pseudo-vector, and is therefore multiplied by -1 under mirror-flip operations. For this reason, the magnetic field *appears* to have opposite symmetry from the electric field, but is really the same.

#### Group Velocities

Given a set of eigenstates at a given k-point, MPB can compute their group velocities (the derivative \(d\omega/d\mathbf{k}\) of frequency with respect to wavevector) using the Hellman-Feynmann theorem. Three functions are provided for this purpose, and we document them here from highest-level to lowest-level.

**`display_group_velocities()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
This is a band function, designed to be passed to `run`, which outputs all of the group velocity vectors (in the Cartesian basis, in units of *c*) at each k-point.

**`ModeSolver.compute_group_velocities()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns a list of group-velocity vectors (in the Cartesian basis, units of *c*) for the bands at the last-computed k-point.

**`ModeSolver.compute_group_velocity_component(direction)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns a list of the group-velocity components (units of *c*) in the given *direction*, one for each band at the last-computed k-point. *direction* is a vector in the reciprocal-lattice basis like the k-points and its length is ignored. This has the advantage of being three times faster than `compute_group_velocities`.

**`ModeSolver.compute_one_group_velocity(which_band)`, `ModeSolver.compute_one_group_velocity_component(direction, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
As above, but returns the group velocity or component thereof only for band `which_band`.

Field Manipulation
------------------

MPB provides a number of ways to take the field of a band and manipulate, process, or output it. These methods usually work in two stages. First, one loads a field into memory, computing it in position space, by calling one of the `get` functions below. Then, other functions can be called to transform or manipulate the field. These fields can also be obtained as numpy arrays which allows one to take advantage of the capabilities of packages like numpy, scipy, and matplotlib.

The simplest class of operations involve only the currently-loaded field, which we describe in the [second subsection](Python_User_Interface.md#loading-and-manipulating-the-current-field) below. To perform more sophisticated operations, involving more than one field, simply get each field as a numpy array.

### Field Normalization

In order to perform useful operations on the fields, it is important to understand how they are normalized. We normalize the fields in the way that is most convenient for perturbation and coupled-mode theory (see [S.G. Johnson et al., (2002)](https://doi.org/10.1103/PhysRevE.65.066611) and also [our online textbook](http://ab-initio.mit.edu/book), ch. 2), so that their energy densities have unit integral. In particular, we normalize the electric (\(\mathbf{E}\)), displacement (\(\mathbf{D} = \varepsilon \mathbf{E}\)) and magnetic (\(\mathbf{H} = -\frac{i}{\omega} \nabla \times \mathbf{E}\)) fields, so that:

-   \(\int \varepsilon |\mathbf{E}|^2  d^3\mathbf{x} = 1\)
-   \(\int |\mathbf{H}|^2 d^3\mathbf{x}= 1\)

where the integrals are over the computational cell. Note the volume element \(d^3\mathbf{x}\) which is the volume of a grid pixel/voxel. If you simply sum \(|\mathbf{H}|^2\) over all the grid points, therefore, you will get (\# grid points) / (volume of cell).

Note that we have dropped the pesky factors of 1/2, π, etcetera from the energy densities, since these do not appear in e.g. perturbation theory, and the fields have arbitrary units anyway. The functions to compute/output energy densities below similarly use \(\varepsilon |\mathbf{E}|^2\) and \(|\mathbf{H}|^2\) without any prefactors.

### Loading and Manipulating the Current Field

In order to load a field into memory, call one of the `get` methods of the `ModeSolver` class. They should only be called after the eigensolver has run or after `init_params`, in the case of `get_epsilon`. One normally calls them after `run`, or in one of the band functions passed to `run`.

**`ModeSolver.get_hfield(which_band, bloch_phase=True)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Loads the magnetic (\(\mathbf{H}\)) field for the band `which_band` and returns an `MPBArray`. To prevent MPB from multiplying by the Bloch phase, pass `bloch_phase=False`.

**`Modesolver.get_dfield(which_band, bloch_phase=True)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Loads the electric displacement (\(\mathbf{D}\)) field for the band `which_band` and returns an `MPBArray`. To prevent MPB from multiplying by the Bloch phase, pass `bloch_phase=False`.

**`ModeSolver.get_efield(which_band, bloch_phase=True)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Loads the electric (\(\mathbf{E}\)) field for the band `which_band` and returns an `MPBArray`. To prevent MPB from multiplying by the Bloch phase, pass `bloch_phase=False`.

**`ModeSolver.get_charge_density(which_band, bloch_phase=True)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Loads the bound charge density \(\nabla \cdot \mathbf{E}\) for the band `which_band`. To prevent MPB from multiplying by the Bloch phase, pass `bloch_phase=False`.

**`ModeSolver.get_epsilon()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Loads the dielectric function and returns a numpy array. Must call `ModeSolver.init_params` or any `run` function before calling `get_epsilon`.

Once loaded, the field can be transformed into another field or a scalar field:

**`ModeSolver.fix_field_phase()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Fix the currently-loaded eigenstate's phase which is normally random in a canonical way, based on the spatial field (\(\mathbf{H}\), \(\mathbf{D}\), or \(\mathbf{E}\)) that has currently been loaded. The phase is fixed to make the real part of the spatial field as big as possible so that you can hopefully visualize just the real part of the field, and a canonical sign is chosen. See also the `fix_Xfield_phase` band functions, above, which are convenient wrappers around `fix_field_phase`.

**`ModeSolver.compute_field_energy()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given the \(\mathbf{H}\) or \(\mathbf{D}\) fields, computes the corresponding energy density function normalized by the total energy in \(\mathbf{H}\) or \(\mathbf{D}\), respectively. Also prints the fraction of the field in each of its Cartesian components in the following form which is suitable for grepping:

```
f-energy-components:, k-index, band-index, x-fraction, y-fraction, z-fraction
```

where `f` is either `h` or `d`. The return value of `compute_field_energy` is a list of 7 numbers: `(U xr xi yr yi zr zi)`. `U` is the total, unnormalized energy, which is in arbitrary units deriving from the normalization of the eigenstate (e.g. the total energy for \(\mathbf{H}\) is always 1.0). `xr` is the fraction of the energy in the real part of the field's x component, `xi` is the fraction in the imaginary part of the x component, etcetera (`yr + yi = y-fraction`, and so on).

**`ModeSolver.compute_field_divergence()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a vector field, compute its divergence.

Various integrals and other information about the eigenstate can be accessed by the following functions, useful e.g. for perturbation theory. Functions dealing with the field vectors require a field to be loaded, and functions dealing with the energy density require an energy density to be loaded via `compute_field_energy`.

**`ModeSolver.compute_energy_in_dielectric(min_eps, max_eps)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns the fraction of the energy that resides in dielectrics with epsilon in the range `min_eps` to `max_eps`.

**`ModeSolver.compute_energy_in_objects(objects)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Returns the fraction of the energy inside zero or more geometric objects in a list.

**`ModeSolver.compute_energy_integral(f)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
`f` is a function `f(u, eps, r)` that returns a number given three parameters: `u`, the energy density at a point; `eps`, the dielectric constant at the same point; and `r`, the position vector in lattice coordinates of the point. `compute_energy_integral` returns the integral of `f` over the unit cell. The integral is computed simply as the sum over the grid points times the volume of a grid pixel/voxel. This can be useful e.g. for perturbation-theory calculations.

**`ModeSolver.compute_field_integral(f)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Like `compute_energy_integral`, but `f` is a function `f(F, eps, r)` that returns a number, possibly complex, where `F` is the complex field vector at the given point.

**`ModeSolver.get_epsilon_point(r)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a position vector *`r`* (in lattice coordinates), return the interpolated dielectric constant at that point. (Since MPB uses a an effective dielectric tensor internally, this actually returns the mean dielectric constant.)

**`ModeSolver.get_epsilon_inverse_tensor_point(r)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a position vector `r` in lattice coordinates, return the interpolated inverse dielectric tensor (a 3x3 matrix) at that point. Near a dielectric interface, the effective dielectric constant is a tensor even if you input only scalar dielectrics; see the [epsilon overview](Developer_Information.md#dielectric-function-computation) for more information. The returned matrix may be complex-Hermetian if you are employing magnetic materials.

**`ModeSolver.get_energy_point(r)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a position vector `r` in lattice coordinates, return the interpolated energy density at that point.

**`ModeSolver.get_field_point(r)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a position vector `r` in lattice coordinates, return the interpolated (complex) field vector at that point.

**`ModeSolver.get_bloch_field_point(r)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a position vector `r` in lattice coordinates, return the interpolated complex Bloch field vector at that point. This is the field without the exp(ikx) envelope.

Finally, we have the following functions to output fields (either the vector fields, the scalar energy density, or epsilon), with the option of outputting several periods of the lattice.

**`ModeSolver.output_field()`**  
**`ModeSolver.output_field_x()`**  
**`ModeSolver.output_field_y()`**  
**`ModeSolver.output_field_z()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the currently-loaded field. For vector fields, `output_field` outputs all of the Cartesian components, while the other variants output only one component.

**`ModeSolver.output_epsilon()`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
A shortcut for calling `get_epsilon` followed by `output_field`. Note that, because epsilon is a tensor, a number of datasets are outputted in `"epsilon.h5"`:

-   `"data"`: 3/trace(1/epsilon)
-   `"epsilon.{xx,xy,xz,yy,yz,zz}"`: the (Cartesian) components of the (symmetric) dielectric tensor.
-   `"epsilon_inverse.{xx,xy,xz,yy,yz,zz}"`: the (Cartesian) components of the (symmetric) inverse dielectric tensor.

### Storing and Combining Multiple Fields

In order to perform operations involving multiple fields, e.g. computing the Poynting vector \(\mathbf{E}^* \times \mathbf{H}\), they must be stored in numpy arrays. Field variables come in three flavors, real-scalar (rscalar) fields, complex-scalar (cscalar) fields, and complex-vector (cvector) fields.

Once you have stored the fields in variables, you probably want to compute something with them. This can be done in three ways: combining fields into new fields with numpy functions (e.g. combine \(\mathbf{E}\) and \(\mathbf{H}\) to \(\mathbf{E}^* \times \mathbf{H}\)), integrating some function of the fields with `integrate_fields` (e.g. to compute coupling integrals for perturbation theory), and getting the field values at arbitrary points with `*_field_get_point` (e.g. to do a line or surface integral). These three functions are described below:

**`ModeSolver.cvector_field_get_point(f, r)`**  
**`ModeSolver.cvector_field_get_point_bloch(f, r)`**  
**`ModeSolver.rscalar_field_get_point(f, r)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Given a position vector `r` in lattice coordinates, return the interpolated field cvector/rscalar from `f` at that point. `cvector_field_get_point_bloch` returns the field *without* the exp(ikx) Bloch wavevector, in analogue to `get_bloch_field_point`.

We also provide functions, in analogue to e.g `get_efield` and `output_efield` above, to "get" various useful functions as the [current field](Python_User_Interface.md#loading-and-manipulating-the-current-field) and to output them to a file:

**`ModeSolver.get_poynting(which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Loads the Poynting vector \(\mathbf{E}^* \times \mathbf{H}\) for the band `which_band`, the flux density of electromagnetic energy flow, as the current field. 1/2 of the real part of this vector is the time-average flux density which can be combined with the imaginary part to determine the amplitude and phase of the time-dependent flux.

**`output_poynting(ms, which_band)`**  
**`output_poynting_x(ms, which_band)`**  
**`output_poynting_y(ms, which_band)`**  
**`output_poynting_z(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the Poynting vector field for `which_band` and `ModeSolver` `ms`; either all or one of the Cartesian components, respectively.

**`ModeSolver.get_tot_pwr(which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Load the time-averaged electromagnetic-field energy density (\(|\mathbf{H}|^2 + \varepsilon |\mathbf{E}|^2\)) for `which_band`. If you multiply the real part of the Poynting vector by a factor of 1/2, above, you should multiply by a factor of 1/4 here for consistency.

**`output_tot_pwr(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the time-averaged electromagnetic-field energy density (above) for `which_band` and `ModeSolver` `ms`.

**`output_charge_density(ms, which_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Output the bound charge density (above) for `which_band`.

#### Stored Fields and Bloch Phases

Complex vector fields like **E** and **H** as computed by MPB are physically of the Bloch form: exp(ikx) times a periodic function. What MPB actually stores, however, is just the periodic function, the Bloch envelope, and only multiplies by exp(ikx) when the fields are output or passed to the user (e.g. in integration functions). This is mostly transparent, with a few exceptions noted above for functions that do not include the exp(ikx) Bloch phase. It is somewhat faster to operate without including the phase.

On some occasions, however, when you get a field and perform some calculation, the resulting field should *not* have any Bloch phase. For example, for the Poynting vector **E**<sup><small>\*</small></sup>x**H**, the exp(ikx) cancels because of the complex conjugation. When creating this sort of field, we must inform MPB not to include the Bloch phase by passing `bloch_phase=False` to the `get_*field` functions. This tells MPB that the field is purely periodic:

Currently, all fields must be either Bloch or non-Bloch (i.e. periodic), which covers most physically meaningful possibilities.

There is another wrinkle: even for fields in Bloch form, the exp(ikx) phase currently always uses the *current* k-point, even if the field was computed from another k-point. So, if you are performing computations combining fields from different k-points, you should take care to always use the periodic envelope of the field, putting the Bloch phase in manually if necessary.

### Manipulating the Raw Eigenvectors

MPB also includes a few low-level routines to manipulate the raw eigenvectors that it computes in a transverse planewave basis.

The most basic operations involve copying, saving, and restoring the current set of eigenvectors or some subset thereof:

**`ModeSolver.get_eigenvectors(first_band, num_bands)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Return a numpy arrray that is a copy of `num_bands` current eigenvectors starting at `first_band`. e.g. to get a copy of all of the eigenvectors, use `ms.get_eigenvectors(1, num_bands)`.

**`ModeSolver.set_eigenvectors(ev, first_band)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Set the current eigenvectors, starting at `first_band`, to those in the `ev` eigenvector object (as returned by `get_eigenvectors`). Does not work if the grid sizes don't match.

**`ModeSolver.load_eigenvectors(filename)`**  
**`ModeSolver.save_eigenvectors(filename)`**  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
Read/write the current eigenvectors (raw planewave amplitudes) to/from an HDF5 file named `filename`. Instead of using `load_eigenvectors` directly, you can pass the `filename` as the *`reset_fields`* parameter of `run_parity`, as [shown above](Python_User_Interface.md#run-functions). Loaded eigenvectors must be of the same size (same grid size and \#bands) as the current settings.

Parallel MPB
------------

TODO

### MPB with MPI Parallelization

MPI Parallelization is currently not supported for the Python interface.

### Alternative Parallelization: mpb-split

TODO
