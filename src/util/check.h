#ifndef CHECK_H
#define CHECK_H

#ifdef CHECK_DISABLE

#define CHECK(condition, message) /* do nothing */

#else

#include <mpi_utils.h>

#define CHECK(condition, message) \
     if (!(condition)) \
          mpi_die("CHECK failure on line %d of " __FILE__ ": " \
		  message "\n", __LINE__)

#endif

#ifdef DEBUG
#define malloc debug_malloc
#define free debug_free
#endif

#endif /* CHECK_H */
