#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "config.h"
#include <check.h>
#include <sphere-quad.h>

#include <maxwell.h>
#include <ctl.h>
#include <ctlgeom.h>

/* return a random number in [0,1]: */
static double mydrand(void)
{
     double d = rand();
     return (d / (double) RAND_MAX);
}

/* return a uniform random number in [a,b] */
static double myurand(double a, double b)
{
     return ((b - a) * mydrand() + a);
}

#define MAX_NSQ_PTS 72
#define NUM_PLANES 100000
#define NUM_OBJECTS 50

#define K_PI 3.141592653589793238462643383279502884197

static double pin(double x, double a, double b)
{
     return (x < a ? a : (x > b ? b : x));
}

/* return the angle, in degrees, between two unit-normalized vectors */
static double angle(vector3 v1, vector3 v2)
{
     double theta = 180/K_PI * acos(pin(vector3_dot(v1,v2), -1,1));
     return (theta > 90 ? 180 - theta : theta);
}

/* return a random unit vector, uniformly distributed over a sphere */
vector3 random_unit_vector3(void)
{
     double z, t, r;
     vector3 v;

     z = 2*mydrand() - 1;
     t = 2*K_PI*mydrand();
     r = sqrt(1 - z*z);
     v.x = r * cos(t);
     v.y = r * sin(t);
     v.z = z;
     return v;
}

double find_edge(geometric_object o, vector3 dir, double max, double tol)
{
     double min = 0;
     CHECK(point_in_fixed_objectp(vector3_scale(min, dir), o) &&
	   !point_in_fixed_objectp(vector3_scale(max, dir), o),
	   "object out of bounds in find_edge");
     do {
	  double d = (min + max) / 2;
          if (point_in_fixed_objectp(vector3_scale(d, dir), o))
	       min = d;
	  else
	       max = d;
     } while (max - min > tol);
     return (min + max) / 2;
}

double normal_err_to_object(geometric_object o, double r, vector3 dir)
{
     int i;
     double d;
     vector3 c, nsum = {0,0,0};
     dir = unit_vector3(dir);
     d = find_edge(o, dir, 2, r * 0.01);
     c = vector3_scale(d, dir);
     for (i = 0; i < num_sphere_quad[2]; ++i) {
	  vector3 v;
	  v.x = sphere_quad[2][i][0] * r;
	  v.y = sphere_quad[2][i][1] * r;
	  v.z = sphere_quad[2][i][2] * r;
	  if (point_in_fixed_objectp(vector3_plus(v, c), o))
	       nsum = vector3_plus(nsum,
				   vector3_scale(sphere_quad[2][i][3], v));
     }
     CHECK(vector3_norm(nsum) > 1e-6, "couldn't get normal vector");
     return angle(unit_vector3(normal_to_object(c, o)), unit_vector3(nsum));
}

void normals_to_object(geometric_object o, double r, int ntrials)
{
     int i;
     double err_mean = 0, err_std = 0, err_max = 0;
     display_geometric_object_info(0, o);
     for (i = 0; i < ntrials; ++i) {
	  double dev;
	  double e = normal_err_to_object(o, r, random_unit_vector3());
	  if (e > err_max) err_max = e;
	  /* stable one-pass formula for mean and std. deviation: */
	  dev = (e - err_mean) / (i + 1);
	  err_mean += dev;
	  err_std += i*(i+1) * dev*dev;
     }
     err_std = sqrt(err_std / (ntrials - 1));
     printf("mean error angle for %d-pt formula = "
	    "%g +/- %g degrees, max error = %g degrees\n\n",
	    num_sphere_quad[2], err_mean, err_std, err_max);
}

static vector3 make_vector3(double x, double y, double z)
{
     vector3 v;
     v.x = x; v.y = y; v.z = z;
     return v;
}

/* return a random geometric object, centered at the origin, with
   diameter roughly 1 */
geometric_object random_object(void)
{
     void* m = NULL;
     vector3 c = { 0, 0, 0 };
     geometric_object o;
     switch (rand() % 5) {
	 case 0:
	      o = make_sphere(m, c, myurand(0.5,1.5));
	      break;
	 case 1:
	      o = make_cylinder(m, c, myurand(0.5,1.5), myurand(0.5,1.5),
				random_unit_vector3());
	      break;
	 case 2:
	      o = make_cone(m, c, myurand(0.5,1.5), myurand(0.5,1.5),
			    random_unit_vector3(), myurand(0.5,1.5));
	      break;
	 case 3:
	      o = make_block(m, c, 
			     random_unit_vector3(),
			     random_unit_vector3(),
			     random_unit_vector3(),
			     make_vector3(myurand(0.5,1.5),
					  myurand(0.5,1.5),
					  myurand(0.5,1.5)));
	      break;
	 case 4:
	      o = make_ellipsoid(m, c, 
				 random_unit_vector3(),
				 random_unit_vector3(),
				 random_unit_vector3(),
				 make_vector3(myurand(0.5,1.5),
					      myurand(0.5,1.5),
					      myurand(0.5,1.5)));
	      break;
     }
     return o;
}

int main(void)
{
     int i, j;
     double err_mean, err_std, err_max;
     double min_angle = 360;
     int missed;

     srand(time(NULL));

     printf("Testing spherical quadratures for normals to %d surfaces.\n\n",
	    NUM_PLANES);
     
     /* compute the minimum angle between pairs of points: */
     for (i = 0; i < num_sphere_quad[2]; ++i)
	  for (j = i + 1; j < num_sphere_quad[2]; ++j) {
	       vector3 v1, v2;
	       double a;
	       v1.x = sphere_quad[2][i][0];
	       v1.y = sphere_quad[2][i][1];
	       v1.z = sphere_quad[2][i][2];
	       v2.x = sphere_quad[2][j][0];
	       v2.y = sphere_quad[2][j][1];
	       v2.z = sphere_quad[2][j][2];
	       a = angle(v1,v2);
	       if (a < min_angle && a > 1e-6)
		    min_angle = a;
	  }
     printf("%d-point formula: minimum angle is %g degrees.\n",
	    num_sphere_quad[2], min_angle);
     
	  /* Normals to planes: */
     err_mean = err_std = err_max = 0.0;
     missed = 0;
     for (i = 0; i < NUM_PLANES; ++i) {
	  vector3 n, nsum = {0,0,0};
	  double d;
	  
	  n = random_unit_vector3();
	  d = mydrand();
	  for (j = 0; j < num_sphere_quad[2]; ++j) {
	       vector3 v;
	       real val;
	       v.x = sphere_quad[2][j][0];
	       v.y = sphere_quad[2][j][1];
	       v.z = sphere_quad[2][j][2];
	       val = vector3_dot(v,n) >= d ? 12.0 : 1.0;
	       val *= sphere_quad[2][j][3];
	       nsum = vector3_plus(nsum, vector3_scale(val, v));
	  }
	  if (vector3_norm(nsum) < 1e-6) {
	       ++missed; --i;
	       continue;
	  }
	  nsum = unit_vector3(nsum);
	  {  /* stable one-pass formula for mean and std. deviation: */
	       double e, dev;
	       e = angle(n, nsum);
	       if (e > err_max) err_max = e;
	       dev = (e - err_mean) / (i + 1);
	       err_mean += dev;
	       err_std += i*(i+1) * dev*dev;
	  }
     }
     err_std = sqrt(err_std / (NUM_PLANES - 1));	  
     printf("planes: mean error angle for %d-pt formula = "
	    "%g +/- %g degrees, max error = %g degrees\n", 
	    num_sphere_quad[2], err_mean, err_std, err_max);
     printf("(Fraction missed = %g)\n", 
	    missed * 1.0 / (NUM_PLANES + missed));
     
     /* Normals to spheres: */
     err_mean = err_std = 0.0;
     missed=0;
     for (i = 0; i < NUM_PLANES; ++i) {
	  vector3 n, nsum = {0,0,0}, c;
	  double r, d;
	  int j;
	  
	  n = random_unit_vector3();
	  d = mydrand() * 0.8 + 0.1;
	  r = 1.0 + mydrand() * 10; /* radius of the sphere */
	  c = vector3_scale(r + d, n);  /* center of the sphere */
	  for (j = 0; j < num_sphere_quad[2]; ++j) {
	       vector3 v;
	       real val;
	       v.x = sphere_quad[2][j][0];
	       v.y = sphere_quad[2][j][1];
	       v.z = sphere_quad[2][j][2];
	       val = vector3_norm(vector3_minus(c,v)) <= r ? 12.0 : 1.0;
	       val *= sphere_quad[2][j][3];
	       nsum = vector3_plus(nsum, vector3_scale(val, v));
	  }
	  nsum =  unit_vector3(nsum);
	  if (vector3_norm(nsum) < 1e-6) {
	       --i;
	       continue;
	  }
	  {  /* stable one-pass formula for mean and std. deviation: */
	       double e, dev;
	       e = angle(n, nsum);
	       if (e > err_max) err_max = e;
	       dev = (e - err_mean) / (i + 1);
	       err_mean += dev;
	       err_std += i*(i+1) * dev*dev;
	  }
     }
     err_std = sqrt(err_std / (NUM_PLANES - 1));	  
     printf("spheres: mean error angle for %d-pt formula = "
	    "%g +/- %g degrees, max error = %g degrees\n", 
	    num_sphere_quad[2], err_mean, err_std, err_max);

     printf("\n");
     
     for (i = 0; i < NUM_OBJECTS; ++i) {
	  geometric_object o = random_object();
	  normals_to_object(o, 0.01, NUM_PLANES/100);
	  geometric_object_destroy(o);
     }
     
     return EXIT_SUCCESS;
}
