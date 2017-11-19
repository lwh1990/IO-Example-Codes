/******************************************************************************
 *                                                                            *
 *  Basic MPI Example - Greetings, Master                                     *
 *                                                                            *
 *  The master process broadcasts the name of the CPU it runs on to           *
 *  the pool. All other processes respond by sending greetings                *
 *  to the master process, which collects the messages and displays           *
 *  them on standard output.                                                  *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *  The original code was written by Gustav at University of Indiana in 2003. *
 *                                                                            *
 *  The current version has been tested/updated by the HPC department at      *
 *  the Norwegian University of Science and Technology in 2011.               *
 *                                                                            *
 ******************************************************************************/
#include <stdio.h>    /* functions sprintf, printf and BUFSIZ defined there */
#include <string.h>   /* function strcpy defined there */
#include <stdlib.h>   /* function exit defined there */
#include <mpi.h>      /* all MPI-2 functions defined there */

#define TRUE 1        
#define FALSE 0
#define MASTER_RANK 0 /* It is traditional to make process 0 the master. */

main(argc, argv)
     int argc;
     char *argv[];
{
  int count, pool_size, my_rank, my_name_length, i_am_the_master = FALSE;
  char my_name[BUFSIZ], master_name[BUFSIZ], send_buffer[BUFSIZ],
    recv_buffer[BUFSIZ];
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Get_processor_name(my_name, &my_name_length);

  if (my_rank == MASTER_RANK) {
    i_am_the_master = TRUE;
    strcpy (master_name, my_name);
  }

  MPI_Bcast(master_name, BUFSIZ, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);

  if (i_am_the_master) 
    for (count = 1; count < pool_size; count++) {
      MPI_Recv (recv_buffer, BUFSIZ, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG,
                MPI_COMM_WORLD, &status);
      printf ("%s\n", recv_buffer);
    }
  else {
    sprintf(send_buffer, "hello %s, greetings from %s, rank = %d",
            master_name, my_name, my_rank);
    MPI_Send (send_buffer, strlen(send_buffer) + 1, MPI_CHAR,
              MASTER_RANK, 0, MPI_COMM_WORLD);
  }

  MPI_Finalize();
   
  exit(0);
}
