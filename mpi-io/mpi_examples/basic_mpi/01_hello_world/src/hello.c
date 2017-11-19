/******************************************************************************
 *                                                                            *
 *  Basic MPI Example - Hello World                                           *
 *                                                                            *
 *  Find the size of the communicator, your rank within it, and the           *
 *  name of the processor you run on.                                         *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *  The original code was written by Gustav at University of Indiana in 2003. *
 *                                                                            *
 *  The current version has been tested/updated by the HPC department at      *
 *  the Norwegian University of Science and Technology in 2011.               *
 *                                                                            *
 ******************************************************************************/
#include <stdio.h>  /* printf and BUFSIZ defined there */
#include <stdlib.h> /* exit defined there */
#include <mpi.h>    /* all MPI-2 functions defined there */

int main(argc, argv)
int argc;
char *argv[];
{
   int rank, size, length;
   char name[BUFSIZ];

   MPI_Init(&argc, &argv);
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Get_processor_name(name, &length);

   printf("%s: hello world from process %d of %d\n", name, rank, size);

   MPI_Finalize();

   exit(0);
}
