Han "Abby" Men (abby.men@gmail.com)
6/6/2016

0. Most of what I have done is limited to two files ./mpb/material_grid_opt_sdp.c and ./mpb/material_grid_opt_lp.c. In addition, configure.ac, ./mpb/Makefile.am, and ./mpb/mpb.scm.in have been modified.

1. To install mosek and obtain an academic license: (a license would expire after about 3 months, but it doesn't seem to have a limit on how many one can keep on requesting)
https://www.mosek.com/resources/downloads

2. To configure and make:

```sh
$ ./configure --with-libctl --with-mosek --with-mpi
$ make
$ make install
```

A few notes on the configuration options:

"--with-mosek" is required for the SDP and LP optimization routines implemented in material_grid_opt_sdp.c material_grid_opt_lp.c;

"--with-mpi" is also (unfortunately) required for the optimization routines because I cannot find the older versions of the code before mpi was enabled. (I also didn't get to implement conditional check because I didn't have the time. In any case, it's probably a much better idea to run the parallel version anyway)

3. To run: `$ mpirun -np 4 mpb-mpi ./examples/SG212.ctl > SG212.out 2>& 1 &`

4: To monitor the optimization status: `$ tail -f SG212.out | grep maximization &`
