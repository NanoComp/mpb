/* Copyright (C) 1999 Massachusetts Institute of Technology.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "../src/config.h"
#include <check.h>

#define NUM_POINTERS 1000
#define MAX_SIZE 1024
#define NUM_MALLOCS 1000000

int main(void)
{
     void *pointers[NUM_POINTERS];
     int i, iter;

#ifdef DEBUG_MALLOC
     printf("Using debug_malloc and debug_free routines.\n");
#else
     fprintf(stderr, 
	     "***** NOTE: malloctest is designed to be run when the package\n"
	     "            is configured with --enable-debug, to test the\n"
	     "            debugging malloc/free routines.\n");
#endif

     srand(time(NULL));

     for (i = 0; i < NUM_POINTERS; ++i)
	  pointers[i] = NULL;

     printf("Doing %d malloc/free calls...\n", NUM_MALLOCS);
     for (iter = 0; iter < NUM_MALLOCS; ++iter) {
	  i = rand() % NUM_POINTERS;
	  if (pointers[i])
	       free(pointers[i]);
	  pointers[i] = malloc(rand() % MAX_SIZE + 1);
	  CHECK(pointers[i], "out of memory!");
	  if ((iter + 1) % (NUM_MALLOCS / 20) == 0)
	       printf("...completed %d...\n", iter + 1);
     }
     printf("Done.\n");

     for (i = 0; i < NUM_POINTERS; ++i)
	  if (pointers[i])
	       free(pointers[i]);

     debug_check_memory_leaks();
     
     printf("Okay.\n");

     return EXIT_SUCCESS;
}
