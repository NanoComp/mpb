/* Copyright (C) 1999, 2000 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <math.h>

#include <ctl.h>

#include "../src/config.h"

#include <check.h>
#include <matrixio.h>

int verbose = 0;

/* a macro to set x = fractional part of x input, xi = integer part,
   with 0 <= x < 1.0. */
#define MODF_POSITIVE(x, xi) { x=modf(x, &xi); if (x<0){ x +=1.0; xi -=1.0; } }

#define ADJ_POINT(i1, i2, nx, dx, xi, xi2) { \
     if (dx >= 0.0) { \
	  i2 = i1 + 1; \
	  if (i2 >= nx) { \
	       i2 -= nx; \
	       xi2 = xi + 1.0; \
	  } \
	  else \
	       xi2 = xi; \
     } \
     else { \
	  i2 = i1 - 1; \
	  if (i2 < 0) { \
	       i2 += nx; \
	       xi2 = xi - 1.0; \
	  } \
	  else \
	       xi2 = xi; \
          dx = -dx; \
     } \
}

void add_cmplx_times_phase(real *sum_re, real *sum_im,
			   real d_re, real d_im,
			   double ix, double iy, double iz, real *s,
			   real scale_by)
{
     static real phase = 0.0, p_re = 1.0, p_im = 0.0;
     real new_phase;

     new_phase = ix * s[0] + iy * s[1] + iz * s[2];
     if (new_phase != phase) {
	  phase = new_phase;
	  p_re = cos(phase);
	  p_im = sin(phase);
     }
     *sum_re += (d_re * p_re - d_im * p_im) * scale_by;
     *sum_im += (d_re * p_im + d_im * p_re) * scale_by;
}

#define TWOPI 6.2831853071795864769252867665590057683943388

#define MAX2(a,b) ((a) >= (b) ? (a) : (b))
#define MIN2(a,b) ((a) < (b) ? (a) : (b))

void map_data(real *d_in_re, real *d_in_im, int n_in[3], 
	      real *d_out_re, real *d_out_im, int n_out[3], 
	      matrix3x3 coord_map,
	      real *kvector,
	      short pick_nearest)
{
     int i, j, k;
     real s[3]; /* phase difference per cell in each lattice direction */
     real min_out_re = 1e20, max_out_re = -1e20, 
	  min_out_im = 1e20, max_out_im = -1e20;

     CHECK(d_in_re && d_out_re, "invalid arguments");
     CHECK((d_out_im && d_in_im) || (!d_out_im && !d_in_im),
	   "both input and output must be real or complex");

     coord_map.c0 = vector3_scale(1.0 / n_out[0], coord_map.c0);
     coord_map.c1 = vector3_scale(1.0 / n_out[1], coord_map.c1);
     coord_map.c2 = vector3_scale(1.0 / n_out[2], coord_map.c2);

     for (i = 0; i < 3; ++i) {
	  if (kvector)
	       s[i] = kvector[i] * TWOPI;
	  else
	       s[i] = 0;
     }

     for (i = 0; i < n_out[0]; ++i)
	  for (j = 0; j < n_out[1]; ++j)
	       for (k = 0; k < n_out[2]; ++k) {
		    real x, y, z;
		    double xi, yi, zi, xi2, yi2, zi2;
		    double dx, dy, dz, mdx, mdy, mdz;
		    int i1, j1, k1, i2, j2, k2;
		    int ijk = (i * n_out[1] + j) * n_out[2] + k;
		    
		    /* find the point corresponding to d_out[i,j,k] in
		       the input array, and also find the next-nearest
		       points. */
		    x = coord_map.c0.x*i + coord_map.c1.x*j + coord_map.c2.x*k;
		    y = coord_map.c0.y*i + coord_map.c1.y*j + coord_map.c2.y*k;
		    z = coord_map.c0.z*i + coord_map.c1.z*j + coord_map.c2.z*k;
		    MODF_POSITIVE(x, xi);
		    MODF_POSITIVE(y, yi);
		    MODF_POSITIVE(z, zi);
		    i1 = x * n_in[0]; j1 = y * n_in[1]; k1 = z * n_in[2];
		    dx = x * n_in[0] - i1;
		    dy = y * n_in[1] - j1;
		    dz = z * n_in[2] - k1;
		    ADJ_POINT(i1, i2, n_in[0], dx, xi, xi2);
		    ADJ_POINT(j1, j2, n_in[1], dy, yi, yi2);
		    ADJ_POINT(k1, k2, n_in[2], dz, zi, zi2);

		    /* dx, mdx, etcetera, are the weights for the various
		       points in the input data, which we use for linearly
		       interpolating to get the output point. */
		    if (pick_nearest) {
			 /* don't interpolate */
			 dx = dx <= 0.5 ? 0.0 : 1.0;
			 dy = dy <= 0.5 ? 0.0 : 1.0;
			 dz = dz <= 0.5 ? 0.0 : 1.0;
		    }
		    mdx = 1.0 - dx;
		    mdy = 1.0 - dy;
		    mdz = 1.0 - dz;
		    
		    /* Now, linearly interpolate the input to get the
		       output.  If the input/output are complex, we
		       also need to multiply by the appropriate phase
		       factor, depending upon which unit cell we are in. */

#define IN_INDEX(i,j,k) ((i * n_in[1] + j) * n_in[2] + k)
		    if (d_out_im) {
			 d_out_re[ijk] = 0.0;
			 d_out_im[ijk] = 0.0;
			 add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk,
					       d_in_re[IN_INDEX(i1,j1,k1)],
					       d_in_im[IN_INDEX(i1,j1,k1)],
					       xi, yi, zi, s,
					       mdx * mdy * mdz);
			 add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk,
					       d_in_re[IN_INDEX(i1,j1,k2)],
					       d_in_im[IN_INDEX(i1,j1,k2)],
					       xi, yi, zi2, s,
					       mdx * mdy * dz);
			 add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk,
					       d_in_re[IN_INDEX(i1,j2,k1)],
					       d_in_im[IN_INDEX(i1,j2,k1)],
					       xi, yi2, zi, s,
					       mdx * dy * mdz);
			 add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk,
					       d_in_re[IN_INDEX(i1,j2,k2)],
					       d_in_im[IN_INDEX(i1,j2,k2)],
					       xi, yi2, zi2, s,
					       mdx * dy * dz);
			 add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk,
					       d_in_re[IN_INDEX(i2,j1,k1)],
					       d_in_im[IN_INDEX(i2,j1,k1)],
					       xi2, yi, zi, s,
					       dx * mdy * mdz);
			 add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk,
					       d_in_re[IN_INDEX(i2,j1,k2)],
					       d_in_im[IN_INDEX(i2,j1,k2)],
					       xi2, yi, zi2, s,
					       dx * mdy * dz);
			 add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk,
					       d_in_re[IN_INDEX(i2,j2,k1)],
					       d_in_im[IN_INDEX(i2,j2,k1)],
					       xi2, yi2, zi, s,
					       dx * dy * mdz);
			 add_cmplx_times_phase(d_out_re + ijk, d_out_im + ijk,
					       d_in_re[IN_INDEX(i2,j2,k2)],
					       d_in_im[IN_INDEX(i2,j2,k2)],
					       xi2, yi2, zi2, s,
					       dx * dy * dz);
			 min_out_im = MIN2(min_out_im, d_out_im[ijk]);
			 max_out_im = MAX2(min_out_im, d_out_im[ijk]);
		    }
		    else {
			 d_out_re[ijk] =
			      d_in_re[IN_INDEX(i1,j1,k1)] * mdx * mdy * mdz +
			      d_in_re[IN_INDEX(i1,j1,k2)] * mdx * mdy * dz +
			      d_in_re[IN_INDEX(i1,j2,k1)] * mdx * dy * mdz +
			      d_in_re[IN_INDEX(i1,j2,k2)] * mdx * dy * dz +
			      d_in_re[IN_INDEX(i2,j1,k1)] * dx * mdy * mdz +
			      d_in_re[IN_INDEX(i2,j1,k2)] * dx * mdy * dz +
			      d_in_re[IN_INDEX(i2,j2,k1)] * dx * dy * mdz +
			      d_in_re[IN_INDEX(i2,j2,k2)] * dx * dy * dz;
		    }
		    min_out_re = MIN2(min_out_re, d_out_re[ijk]);
		    max_out_re = MAX2(min_out_re, d_out_re[ijk]);
#undef IN_INDEX
	       }

     if (verbose) {
	  printf("real part range: %g .. %g\n", min_out_re, max_out_re);
	  if (d_out_im)
	       printf("imag part range: %g .. %g\n", min_out_im, max_out_im);
     }
}

void handle_dataset(matrixio_id in_file, matrixio_id out_file, 
		    const char *name_re, const char *name_im,
		    matrix3x3 Rout, matrix3x3 coord_map,
		    real *kvector, int resolution, real multiply_size[3],
		    int pick_nearest)
{
     real *d_in_re = NULL, *d_in_im = NULL, *d_out_re = NULL, *d_out_im = NULL;
     int in_dims[3] = {1,1,1}, out_dims[3] = {1,1,1}, rank = 3, i, N;
     int start[3] = {0,0,0};
     matrixio_id data_id;
     char out_name[1000];

     d_in_re = matrixio_read_real_data(in_file, name_re, &rank, in_dims,
				       0, 0, 0, NULL);
     if (!d_in_re)
	  goto done;

     if (verbose)
	  printf("Found dataset %s...\n", name_re);

     if (name_im) {
	  d_in_im = matrixio_read_real_data(in_file, name_im, &rank, out_dims,
					    0, 0, 0, NULL);
	  if (!d_in_im) {
	       fprintf(stderr, "mpb-data: found %s dataset but not %s\n",
		       name_re, name_im);
	       goto done;
	  }
	  
	  for (i = 0; i < 3; ++i) {
	       CHECK(out_dims[i] == in_dims[i],
		     "re/im datasets must have same size!");
	  }

	  if (verbose)
	       printf("   and imaginary part dataset %s...\n", name_im);

     }

     if (verbose)
	  printf("Input data is rank %d, size %dx%dx%d.\n",
		 rank, in_dims[0], in_dims[1], in_dims[2]);

     if (resolution > 0) {
	  out_dims[0] = vector3_norm(Rout.c0) * resolution + 0.5;
	  out_dims[1] = vector3_norm(Rout.c1) * resolution + 0.5;
	  out_dims[2] = vector3_norm(Rout.c2) * resolution + 0.5;
     }
     else {
	  for (i = 0; i < 3; ++i)
	       out_dims[i] = in_dims[i] * multiply_size[i];
     }
     for (i = rank; i < 3; ++i)
	  out_dims[i] = 1;
     for (N = 1, i = 0; i < 3; ++i)
	  N *= (out_dims[i] = MAX2(out_dims[i], 1));

     if (verbose)
	  printf("Output data %dx%dx%d.\n",
		 out_dims[0], out_dims[1], out_dims[2]);

     CHK_MALLOC(d_out_re, real, N);
     if (d_in_im) {
	  CHK_MALLOC(d_out_im, real, N);
     }

     map_data(d_in_re, d_in_im, in_dims, d_out_re, d_out_im, out_dims,
	      coord_map, kvector, pick_nearest);

     strcpy(out_name, name_re);
     if (out_file == in_file)
	  strcat(out_name, "-new");
     if (verbose)
	  printf("Writing dataset to %s...\n", out_name);
     data_id = matrixio_create_dataset(out_file, out_name,"", rank, out_dims);
     matrixio_write_real_data(data_id, out_dims, start, 1, d_out_re);
     matrixio_close_dataset(data_id);

     if (d_out_im) {
	  strcpy(out_name, name_im);
	  if (out_file == in_file)
	       strcat(out_name, "-new");
	  if (verbose)
	       printf("Writing dataset to %s...\n", out_name);
	  data_id = matrixio_create_dataset(out_file, out_name, "",
					    rank, out_dims);
	  matrixio_write_real_data(data_id, out_dims, start, 1, d_out_im);
	  matrixio_close_dataset(data_id);
     }

     if (verbose)
	  printf("Successfully wrote out data.\n");

 done:
     free(d_in_re);
     free(d_in_im);
     free(d_out_re);
     free(d_out_im);
}

void handle_file(const char *fname, const char *out_fname,
		 const char *data_name,
		 int rectify, int resolution, real multiply_size[3],
		 int pick_nearest)
{
     matrixio_id in_file, out_file;
     real *R, *kvector, *copies;
     int dims[2], rank;
     matrix3x3 Rin = {{1,0,0},{0,1,0},{0,0,1}}, Rout, coord_map;
#define NUM_DATANAMES 4
     char datanames[NUM_DATANAMES][30] = {
	  "data", "x", "y", "z",
     };
     int i;

     if (verbose)
	  printf("Reading file %s...\n", fname);
     in_file = matrixio_open(fname, out_fname != NULL);

     if (data_name && !data_name[0])
	  data_name = NULL;

     R = matrixio_read_data_attr(in_file, "lattice vectors",
				 &rank, 2, dims);
     if (R && rank == 2 && dims[0] == 3 && dims[1] == 3) {
	  Rin.c0.x = R[0*3+0]; Rin.c0.y = R[0*3+1]; Rin.c0.z = R[0*3+2];
	  Rin.c1.x = R[1*3+0]; Rin.c1.y = R[1*3+1]; Rin.c1.z = R[1*3+2];
	  Rin.c2.x = R[2*3+0]; Rin.c2.y = R[2*3+1]; Rin.c2.z = R[2*3+2];
	  if (verbose)
	       printf("Read lattice vectors.\n");
     }
     free(R);

     kvector = matrixio_read_data_attr(in_file, "Bloch wavevector",
				       &rank, 1, dims);
     if (rank != 1 || dims[0] != 3) {
	  free(kvector);
	  kvector = NULL;
     }
     else if (verbose)
	  printf("Read Bloch wavevector (%g, %g, %g)\n",
		 kvector[0], kvector[1], kvector[2]);
     
     copies = matrixio_read_data_attr(in_file, "lattice copies",
				      &rank, 1, dims);
     if (copies && rank == 1 && dims[0] == 3) {
	  Rin.c0 = vector3_scale(copies[0], Rin.c0);
	  Rin.c1 = vector3_scale(copies[1], Rin.c1);
	  Rin.c2 = vector3_scale(copies[2], Rin.c2);
	  if (kvector) {
	       kvector[0] *= copies[0];
	       kvector[1] *= copies[1];
	       kvector[2] *= copies[2];
	  }
	  if (verbose)
	       printf("Read lattice copies (%g, %g, %g)\n",
		      copies[0], copies[1], copies[2]);
     }
     free(copies);

     if (verbose)
	  printf("Input lattice = (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n",
		 Rin.c0.x, Rin.c0.y, Rin.c0.z,
		 Rin.c1.x, Rin.c1.y, Rin.c1.z,
		 Rin.c2.x, Rin.c2.y, Rin.c2.z);

     Rout.c0 = vector3_scale(multiply_size[0], Rin.c0);
     Rout.c1 = vector3_scale(multiply_size[1], Rin.c1);
     Rout.c2 = vector3_scale(multiply_size[2], Rin.c2);
     if (rectify) {
	  /* orthogonalize the output lattice vectors; note that this
	     preserves the volume of the output cell. */
	  Rout.c1 = vector3_minus(Rout.c1, 
				vector3_scale(vector3_dot(Rout.c0, Rout.c1) /
					      vector3_dot(Rout.c0, Rout.c0),
					      Rout.c0));
	  Rout.c2 = vector3_minus(Rout.c2, 
				vector3_scale(vector3_dot(Rout.c0, Rout.c2) /
					      vector3_dot(Rout.c0, Rout.c0),
					      Rout.c0));
	  Rout.c2 = vector3_minus(Rout.c2, 
				vector3_scale(vector3_dot(Rout.c1, Rout.c2) /
					      vector3_dot(Rout.c1, Rout.c1),
					      Rout.c1));
     }

     if (verbose)
	  printf("Output lattice = (%g,%g,%g), (%g,%g,%g), (%g,%g,%g)\n",
		 Rout.c0.x, Rout.c0.y, Rout.c0.z,
		 Rout.c1.x, Rout.c1.y, Rout.c1.z,
		 Rout.c2.x, Rout.c2.y, Rout.c2.z);

     coord_map = matrix3x3_mult(matrix3x3_inverse(Rin), Rout);

     if (out_fname) {
	  if (verbose)
	       printf("Creating output file %s...\n", out_fname);
	  out_file = matrixio_create(out_fname);
     }
     else {
	  if (verbose)
	       printf("Writing output datasets to input file %s...\n", fname);
	  out_file = in_file;
     }

     for (i = 0; i < NUM_DATANAMES; ++i) {
	  const char *dname = datanames[i];
	  char name_re[300], name_im[300];

	  if (data_name)
	       dname = data_name;

	  strcpy(name_re, dname);
	  handle_dataset(in_file, out_file, name_re, NULL,
			 Rout, coord_map, kvector, resolution,
			 multiply_size, pick_nearest);

	  sprintf(name_re, "%s.r", dname);
	  sprintf(name_im, "%s.i", dname);
	  handle_dataset(in_file, out_file, name_re, name_im,
			 Rout, coord_map, kvector, resolution,
			 multiply_size, pick_nearest);

	  if (data_name)
	       break;
     }

     free(kvector);

     matrixio_close(in_file);
     if (out_file != in_file)
	  matrixio_close(out_file);
}

void usage(FILE *f)
{
     fprintf(f, "Usage: mpb-data [options] [<filenames>]\n"
	     "Options:\n"
             "         -h : this help message\n"
             "         -V : print version number and copyright\n"
             "         -v : verbose output\n"
	     "  -o <file> : output to <file> (first input file only)\n"
	     "         -r : output rectangular cell\n"
	     "     -n <n> : output resolution of n grid points per a\n"
	     "    -x <mx>\n"
	     "    -y <my>\n"
	     "    -z <mx> : output mx/my/mz periods in the x/y/z directions\n"
	     "     -m <s> : same as -x <s> -y <s> -z <s>\n"
	     "         -p : pixellized output (no grid interpolation)\n"
	     "  -d <name> : use dataset <name> in the input files (default: all mpb datasets)\n"
	     "           -- you can also specify a dataset via <filename>:<name>\n"
	  );
}

/* given an fname of the form <filename>:<data_name>, return a pointer
   to a newly-allocated string containing <filename>, and point data_name
   to the position of <data_name> in fname.  The user must free() the
<filename> string. */
static char *split_fname(char *fname, char **data_name)
{
     int fname_len;
     char *colon, *filename;

     fname_len = strlen(fname);
     colon = strchr(fname, ':');
     if (colon) {
          int colon_len = strlen(colon);
          filename = (char*) malloc(sizeof(char) * (fname_len-colon_len+1));
          CHECK(filename, "out of memory");
          strncpy(filename, fname, fname_len-colon_len+1);
          filename[fname_len-colon_len] = 0;
          *data_name = colon + 1;
     }
else { /* treat as if ":" were at the end of fname */
          filename = (char*) malloc(sizeof(char) * (fname_len + 1));
          CHECK(filename, "out of memory");
          strcpy(filename, fname);
          *data_name = fname + fname_len;
     }
     return filename;
}

int main(int argc, char **argv)
{
     char *out_fname = NULL, *data_name = NULL;
     int rectify = 0, resolution = 0;
     real multiply_size[3] = {1,1,1};
     int pick_nearest = 0;
     int ifile, c;
     extern char *optarg;
     extern int optind;

     while ((c = getopt(argc, argv, "hVvo:x:y:z:m:d:n:pr")) != -1)
          switch (c) {
              case 'h':
                   usage(stdout);
                   return EXIT_SUCCESS;
              case 'V':
                   printf("mpb-data " MPB_VERSION " by Steven G. Johnson.\n"
"Copyright (C) 1999, 2000 Massachusetts Institute of Technology.\n"
"This is free software, and you are welcome to redistribute it under the\n"
"terms of the GNU General Public License (GPL).  mpb-data comes with\n"
"ABSOLUTELY NO WARRANTY; see the GPL for more details.\n");
                   return EXIT_SUCCESS;
              case 'v':
                   verbose = 1;
                   break;
              case 'o':
		   free(out_fname);
                   out_fname = (char*) malloc(sizeof(char) *
                                              (strlen(optarg) + 1));
                   CHECK(out_fname, "out of memory");
                   strcpy(out_fname, optarg);
                   break;
              case 'd':
		   free(data_name);
                   data_name = (char*) malloc(sizeof(char) *
                                              (strlen(optarg) + 1));
                   CHECK(data_name, "out of memory");
                   strcpy(data_name, optarg);
                   break;
              case 'x':
                   multiply_size[0] = atof(optarg);
                   break;
              case 'y':
                   multiply_size[1] = atof(optarg);
                   break;
              case 'z':
                   multiply_size[2] = atof(optarg);
                   break;
              case 'm':
                   multiply_size[0] = atof(optarg);
                   multiply_size[1] = atof(optarg);
                   multiply_size[2] = atof(optarg);
                   break;
              case 'n':
                   resolution = atoi(optarg);
                   break;
              case 'p':
                   pick_nearest = 1;
                   break;
              case 'r':
                   rectify = 1;
                   break;
              default:
                   fprintf(stderr, "Invalid argument -%c\n", c);
                   usage(stderr);
                   return EXIT_FAILURE;
          }
     if (optind == argc) {  /* no parameters left */
          usage(stderr);
          return EXIT_FAILURE;
     }
     
     for (ifile = optind; ifile < argc; ++ifile) {
	  char *dname, *h5_fname;
          h5_fname = split_fname(argv[ifile], &dname);
	  if (!dname[0])
               dname = data_name;

	  handle_file(h5_fname, out_fname, dname, rectify, resolution, 
		      multiply_size, pick_nearest);
	  
	  if (out_fname)
               free(out_fname);
          out_fname = NULL;
          free(h5_fname);
     }
     free(data_name);

     return EXIT_SUCCESS;
}
