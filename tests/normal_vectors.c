#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "../src/config.h"
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

int main(void)
{
     int i, nsq, num_sq_pts[] = { 12, 50, 72 };
     real x[MAX_NSQ_PTS], y[MAX_NSQ_PTS], z[MAX_NSQ_PTS], w[MAX_NSQ_PTS];

     srand(time(NULL));

     printf("Testing spherical quadratures to find normals to %d surfaces.\n",
	    NUM_PLANES);

     for (nsq = 0; nsq < sizeof(num_sq_pts) / sizeof(int); ++nsq) {
	  double err_mean, err_std;

	  spherical_quadrature_points(x,y,z,w, num_sq_pts[nsq]);

	  /* Normals to planes: */
	  err_mean = err_std = 0.0;
	  for (i = 0; i < NUM_PLANES; ++i) {
	       vector3 n, nsum = {0,0,0};
	       double d;
	       int j;
	       
	       n.x = mydrand() - 0.5;
	       n.y = mydrand() - 0.5;
	       n.z = mydrand() - 0.5;
	       n = unit_vector3(n);
	       d = mydrand();
	       for (j = 0; j < num_sq_pts[nsq]; ++j) {
		    vector3 v;
		    real val;
		    v.x = x[j]; v.y = y[j]; v.z = z[j];
		    val = vector3_dot(v,n) >= d ? 12.0 : 1.0;
		    val *= w[j];
		    nsum = vector3_plus(nsum, vector3_scale(val, v));
	       }
	       nsum =  unit_vector3(nsum);
	       nsum = vector3_minus(n, nsum);
	       {  /* stable one-pass formula for mean and std. deviation: */
		    double e, dev;
		    e = sqrt(vector3_dot(nsum, nsum) / 3); /* rms error */
		    dev = (e - err_mean) / (i + 1);
		    err_mean += dev;
		    err_std += i*(i+1) * dev*dev;
	       }
	  }
	  err_std = sqrt(err_std / (NUM_PLANES - 1));	  
	  printf("planes: mean error for %d-point formula = %g +/- %g\n", 
		 num_sq_pts[nsq], err_mean, err_std);

	  /* Normals to spheres: */
	  err_mean = err_std = 0.0;
	  for (i = 0; i < NUM_PLANES; ++i) {
	       vector3 n, nsum = {0,0,0}, c;
	       double r, d;
	       int j;
	       
	       n.x = mydrand() - 0.5;
	       n.y = mydrand() - 0.5;
	       n.z = mydrand() - 0.5;
	       n = unit_vector3(n);
	       d = mydrand();
	       do {
		    r = mydrand() * 10; /* radius of the sphere */
	       } while (r + d < 1.0); /* require sphere to intersect surface */
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
	       nsum = vector3_minus(n, nsum);
	       {  /* stable one-pass formula for mean and std. deviation: */
		    double e, dev;
		    e = sqrt(vector3_dot(nsum, nsum) / 3); /* rms error */
		    dev = (e - err_mean) / (i + 1);
		    err_mean += dev;
		    err_std += i*(i+1) * dev*dev;
	       }
	  }
	  err_std = sqrt(err_std / (NUM_PLANES - 1));	  
	  printf("spheres: mean error for %d-point formula = %g +/- %g\n", 
		 num_sq_pts[nsq], err_mean, err_std);
     }

     return EXIT_SUCCESS;
}
