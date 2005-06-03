#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "config.h"
#include <check.h>

#include <maxwell.h>
#include <ctl.h>

/* return a random number in [0,1]: */
double mydrand(void)
{
     double d = rand();
     return (d / (double) RAND_MAX);
}

#define MAX_NSQ_PTS 72
#define NUM_PLANES 100000

#define K_PI 3.141592653589793238462643383279502884197

/* return the angle, in degrees, between two unit-normalized vectors */
double angle(vector3 v1, vector3 v2)
{
     double theta = 180/K_PI * acos(vector3_dot(v1,v2));
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


int main(void)
{
     int i, j, nsq, num_sq_pts[] = { 12, 50, 72 };
     real x[MAX_NSQ_PTS], y[MAX_NSQ_PTS], z[MAX_NSQ_PTS], w[MAX_NSQ_PTS];

     srand(time(NULL));

     printf("Testing spherical quadratures for normals to %d surfaces.\n\n",
	    NUM_PLANES);

     for (nsq = 0; nsq < sizeof(num_sq_pts) / sizeof(int); ++nsq) {
	  double err_mean, err_std, err_max;
	  double min_angle = 360;
	  int missed;

	  spherical_quadrature_points(x,y,z,w, num_sq_pts[nsq]);

	  /* compute the minimum angle between pairs of points: */
	  for (i = 0; i < num_sq_pts[nsq]; ++i)
	       for (j = i + 1; j < num_sq_pts[nsq]; ++j) {
		    vector3 v1, v2;
		    double a;
		    v1.x = x[i]; v1.y = y[i]; v1.z = z[i];
		    v2.x = x[j]; v2.y = y[j]; v2.z = z[j];
		    a = angle(v1,v2);
		    if (a < min_angle && a > 1e-6)
			 min_angle = a;
	       }
	  printf("%d-point formula: minimum angle is %g degrees.\n",
		 num_sq_pts[nsq], min_angle);

	  /* Normals to planes: */
	  err_mean = err_std = err_max = 0.0;
	  missed = 0;
	  for (i = 0; i < NUM_PLANES; ++i) {
	       vector3 n, nsum = {0,0,0};
	       double d;
	       
	       n = random_unit_vector3();
	       d = mydrand();
	       for (j = 0; j < num_sq_pts[nsq]; ++j) {
		    vector3 v;
		    real val;
		    v.x = x[j]; v.y = y[j]; v.z = z[j];
		    val = vector3_dot(v,n) >= d ? 12.0 : 1.0;
		    val *= w[j];
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
		 num_sq_pts[nsq], err_mean, err_std, err_max);
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
	       for (j = 0; j < num_sq_pts[nsq]; ++j) {
		    vector3 v;
		    real val;
		    v.x = x[j]; v.y = y[j]; v.z = z[j];
		    val = vector3_norm(vector3_minus(c,v)) <= r ? 12.0 : 1.0;
		    val *= w[j];
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
		 num_sq_pts[nsq], err_mean, err_std, err_max);

	  printf("\n");
     }

     return EXIT_SUCCESS;
}
