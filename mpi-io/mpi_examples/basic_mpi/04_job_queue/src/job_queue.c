/******************************************************************************
 *                                                                            *
 *  Basic MPI Example - Job Queue                                             *
 *                                                                            *
 *  Demonstrates job scheduling through performing the matrix vector product  *
 *  c = A*b. The master assigns rows of A to the slaves to compute            *
 *  each entry of c.                                                          *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *  The original code was written by Gustav at University of Indiana in 2003. *
 *                                                                            *
 *  The current version has been tested/updated by the HPC department at      *
 *  the Norwegian University of Science and Technology in 2011.               *
 *                                                                            *
 ******************************************************************************/
#include <stdio.h>     /* [fs]printf, fopen and fclose defined here */
#include <stdlib.h>    /* exit defined here */
#include <sys/types.h> /* chmod defined here */
#include <sys/stat.h>  /* chmod defined here */
#include <mpi.h>

#define COLS 100
#define ROWS 100
#define TRUE 1
#define FALSE 0
#define MASTER_RANK 0

int main ( int argc, char **argv )
{
   int pool_size, my_rank, destination;
   int i_am_the_master = FALSE; 
   int A[ROWS][COLS], b[COLS], c[ROWS], i, j;
   int int_buffer[BUFSIZ];
   char my_logfile_name[BUFSIZ];
   FILE *my_logfile;
   MPI_Status status;
   

   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
   MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

   if (my_rank == MASTER_RANK) i_am_the_master = TRUE;

   sprintf(my_logfile_name, "output/file.%d", my_rank);
   my_logfile = fopen(my_logfile_name, "w");
   (void) chmod(my_logfile_name, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);

   if (i_am_the_master) {

      int row, count, sender;

      for (j = 0; j < COLS; j++) {
         b[j] = 1;
         for (i = 0; i < ROWS; i++) A[i][j] = i;
      }

      MPI_Bcast(b, COLS, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

      count = 0;
      for (destination = 0; destination < pool_size; destination++) {
         if (destination != my_rank) {
            for (j = 0; j < COLS; j++) int_buffer[j] = A[count][j];
            MPI_Send(int_buffer, COLS, MPI_INT, destination, count,
                     MPI_COMM_WORLD);
            fprintf(my_logfile, "sent row %d to %d\n", count, destination);
            count = count + 1;
         }
      }

      for (i = 0; i < ROWS; i++) {
         MPI_Recv (int_buffer, BUFSIZ, MPI_INT, MPI_ANY_SOURCE, 
                   MPI_ANY_TAG, MPI_COMM_WORLD, &status);
         sender = status.MPI_SOURCE;
         row = status.MPI_TAG;
         c[row] = int_buffer[0];
         fprintf(my_logfile, "\treceived row %d from %d\n", row, sender);
         if (count < ROWS) {
            for (j = 0; j < COLS; j++) int_buffer[j] = A[count][j];
            MPI_Send(int_buffer, COLS, MPI_INT, sender, count,
                     MPI_COMM_WORLD);
            fprintf(my_logfile, "sent row %d to %d\n", count, sender);
            count = count + 1;
         }
         else {
            MPI_Send(NULL, 0, MPI_INT, sender, ROWS, MPI_COMM_WORLD);
            fprintf(my_logfile, "terminated process %d with tag %d\n", sender, ROWS);
         }
      }
      for (row = 0; row < ROWS; row++) printf("%d ", c[row]);
      printf("\n");

   }
   else { /* I am not the master */

      int sum, row;

      MPI_Bcast(b, COLS, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
      fprintf(my_logfile, "received broadcast from %d\n", MASTER_RANK);
      MPI_Recv(int_buffer, COLS, MPI_INT, MASTER_RANK, MPI_ANY_TAG,
                    MPI_COMM_WORLD, &status);
      fprintf(my_logfile, "received a message from %d, tag %d\n",
                   status.MPI_SOURCE, status.MPI_TAG);
      while (status.MPI_TAG != ROWS) { /* The job is not finished */
         row = status.MPI_TAG; sum = 0;
         for (i = 0; i < COLS; i++) sum = sum + int_buffer[i] * b[i];
         int_buffer[0] = sum;
         MPI_Send (int_buffer, 1, MPI_INT, MASTER_RANK, row, MPI_COMM_WORLD);
         fprintf(my_logfile, "sent row %d to %d\n", row, MASTER_RANK);
         MPI_Recv (int_buffer, COLS, MPI_INT, MASTER_RANK, MPI_ANY_TAG,
                   MPI_COMM_WORLD, &status);
         fprintf(my_logfile, "received a message from %d, tag %d\n",
                 status.MPI_SOURCE, status.MPI_TAG);
      }
      fprintf(my_logfile, "exiting on  tag %d\n", status.MPI_TAG);
   }

   fclose (my_logfile);

   MPI_Finalize ();

   exit (0);
}
