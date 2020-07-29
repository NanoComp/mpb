/* This is a semi-horrifying macro to perform we have different loops over
   the coordinates, depending upon whether we are using complex or real and serial or
   parallel transforms.  Each loop will define, in its body,
   variables (i1,i2,i3) describing the global coordinate of the current
   point, and xyz_index describing the corresponding index in
      the array md->eps_inv[] or similar.

We use it like this:

LOOP_XYZ(md) { // {{ implied open braces (yuck!)
     body...
}}}

where md is the maxwell_data pointer.
 */

#ifdef SCALAR_COMPLEX
#ifndef HAVE_MPI
#define LOOP_XYZ(md)                                                                               \
  {                                                                                                \
    int n1 = md->nx, n2 = md->ny, n3 = md->nz, i1, i2, i3;                                         \
    for (i1 = 0; i1 < n1; ++i1)                                                                    \
      for (i2 = 0; i2 < n2; ++i2)                                                                  \
        for (i3 = 0; i3 < n3; ++i3) {                                                              \
          int xyz_index = ((i1 * n2 + i2) * n3 + i3);
#else /* HAVE_MPI */
/* first two dimensions are transposed in MPI output: */
#define LOOP_XYZ(md)                                                                               \
  {                                                                                                \
    int n1 = md->nx, n2 = md->ny, n3 = md->nz, i1, i2_, i3;                                        \
    int local_n2 = md->local_ny, local_y_start = md->local_y_start;                                \
    for (i2_ = 0; i2_ < local_n2; ++i2_)                                                           \
      for (i1 = 0; i1 < n1; ++i1)                                                                  \
        for (i3 = 0; i3 < n3; ++i3) {                                                              \
          int i2 = i2_ + local_y_start;                                                            \
          int xyz_index = ((i2_ * n1 + i1) * n3 + i3);
#endif /* HAVE_MPI */
#else  /* not SCALAR_COMPLEX */
#ifndef HAVE_MPI
#define LOOP_XYZ(md)                                                                               \
  {                                                                                                \
    int n1 = md->nx, n2 = md->ny, n3 = md->nz, i1_, i2_, i1, i2, i3;                               \
    int n_other = md->other_dims;                                                                  \
    int n_last = md->last_dim_size / 2;                                                            \
    int rank = (n3 == 1) ? (n2 == 1 ? 1 : 2) : 3;                                                  \
    for (i1_ = 0; i1_ < n_other; ++i1_)                                                            \
      for (i2_ = 0; i2_ < n_last; ++i2_) {                                                         \
        int xyz_index = i1_ * n_last + i2_;                                                        \
        switch (rank) {                                                                            \
          case 2:                                                                                  \
            i1 = i1_;                                                                              \
            i2 = i2_;                                                                              \
            i3 = 0;                                                                                \
            break;                                                                                 \
          case 3:                                                                                  \
            i1 = i1_ / n2;                                                                         \
            i2 = i1_ % n2;                                                                         \
            i3 = i2_;                                                                              \
            break;                                                                                 \
          default:                                                                                 \
            i1 = i2_;                                                                              \
            i2 = i3 = 0;                                                                           \
            break;                                                                                 \
        }

#else /* HAVE_MPI */
#define LOOP_XYZ(md)                                                                               \
  {                                                                                                \
    int n1 = md->nx, n2 = md->ny, n3 = md->nz, i1, i2_, i3;                                        \
    int local_n2 = md->local_ny, local_y_start = md->local_y_start;                                \
    int local_n3 = n3 > 1 ? md->last_dim_size / 2 : 1;                                             \
    for (i2_ = 0; i2_ < local_n2; ++i2_)                                                           \
      for (i1 = 0; i1 < n1; ++i1)                                                                  \
        for (i3 = 0; i3 < local_n3; ++i3) {                                                        \
          int i2 = i2_ + local_y_start;                                                            \
          int xyz_index = ((i2_ * n1 + i1) * local_n3 + i3);
#endif /* HAVE_MPI */

#endif /* not SCALAR_COMPLEX */
