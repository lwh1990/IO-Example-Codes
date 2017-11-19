/******************************************************************************
 *                                                                            *
 *  Basic MPI Example - Diffusion                                             *
 *                                                                            *
 *  Uses a cartesian grid to exchange boundary data between processes.        *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *  The original code was written by Gustav at University of Indiana in 2003. *
 *                                                                            *
 *  The current version has been tested/updated by the HPC department at      *
 *  the Norwegian University of Science and Technology in 2011.               *
 *                                                                            *
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

#define FALSE 0
#define TRUE  1
#define MASTER_RANK 0

#define UPDOWN              0
#define SIDEWAYS            1
#define RIGHT               1
#define UP                  1

#define PROCESS_DIMENSIONS  2

#define PROCESS_ROWS        4
#define ROWS                4
#define DISPLAY_ROWS       16      /* must be PROCESS_ROWS * ROWS */

#define PROCESS_COLUMNS     3
#define COLUMNS             4
#define DISPLAY_COLUMNS    12      /* must be PROCESS_COLUMNS * COLUMNS */

int main ( int argc, char **argv )
{
  int pool_size, my_rank, destination, source;
  MPI_Status status;
  char char_buffer[BUFSIZ];
  int i_am_the_master = FALSE; 

  int divisions[PROCESS_DIMENSIONS] = {PROCESS_ROWS, PROCESS_COLUMNS};
  int periods[PROCESS_DIMENSIONS] = {0, 0};
  int reorder = 1;
  MPI_Comm cartesian_communicator;
  int my_cartesian_rank, my_coordinates[PROCESS_DIMENSIONS];
  int left_neighbour, right_neighbour, bottom_neighbour, top_neighbour;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == MASTER_RANK) i_am_the_master = TRUE;

  MPI_Cart_create ( MPI_COMM_WORLD, PROCESS_DIMENSIONS, divisions, 
		    periods, reorder, &cartesian_communicator );

  if (cartesian_communicator != MPI_COMM_NULL) {

    int matrix [ROWS][COLUMNS];
    int i, j;
    MPI_Datatype column_type;

    MPI_Comm_rank ( cartesian_communicator, &my_cartesian_rank );
    MPI_Cart_coords ( cartesian_communicator, my_cartesian_rank,
		      PROCESS_DIMENSIONS, my_coordinates );
    MPI_Cart_shift ( cartesian_communicator, SIDEWAYS, RIGHT, 
		     &left_neighbour, &right_neighbour );
    MPI_Cart_shift ( cartesian_communicator, UPDOWN, UP,
		     &bottom_neighbour, &top_neighbour );

    if (! i_am_the_master ) {
      sprintf(char_buffer, "process %2d, cartesian %2d, \
coords (%2d,%2d), left %2d, right %2d, top %2d, bottom %2d",
	      my_rank, my_cartesian_rank, my_coordinates[0],
	      my_coordinates[1], left_neighbour, right_neighbour,
	      top_neighbour, bottom_neighbour);
      MPI_Send(char_buffer, strlen(char_buffer) + 1, MPI_CHAR,
	       MASTER_RANK, 3003, MPI_COMM_WORLD);
    } 
    else {

      int number_of_c_procs, count;

      number_of_c_procs = divisions[0] * divisions[1]; 
      for (count = 0; count < number_of_c_procs - 1; count++) {
	MPI_Recv(char_buffer, BUFSIZ, MPI_CHAR, MPI_ANY_SOURCE, 3003,
		 MPI_COMM_WORLD, &status);
	printf ("%s\n", char_buffer);
      }
      printf( "process %2d, cartesian %2d, \
coords (%2d,%2d), left %2d, right %2d, top %2d, bottom %2d\n",
	      my_rank, my_cartesian_rank, my_coordinates[0],
	      my_coordinates[1], left_neighbour, right_neighbour,
	      top_neighbour, bottom_neighbour);
    }

    for ( i = 0; i < ROWS; i++ ) {
      for ( j = 0; j < COLUMNS; j++ ) {
	matrix [i][j] = my_cartesian_rank;
      }
    }

    if (my_cartesian_rank != MASTER_RANK ) 
      MPI_Send ( matrix, COLUMNS * ROWS, MPI_INT, MASTER_RANK, 3003,
		 cartesian_communicator );
    else
      collect_matrices ( cartesian_communicator, my_cartesian_rank,
			 matrix, 3003 );

    MPI_Sendrecv ( &matrix[ROWS - 2][0], COLUMNS, MPI_INT, 
		   top_neighbour, 4004,
		   &matrix[0][0], COLUMNS, MPI_INT, bottom_neighbour, 
		   4004, 
		   cartesian_communicator, &status );

    MPI_Sendrecv ( &matrix[1][0], COLUMNS, MPI_INT, bottom_neighbour,
		   5005,
		   &matrix[ROWS - 1][0], COLUMNS, MPI_INT, 
		   top_neighbour, 5005,
		   cartesian_communicator, &status );

    if (my_cartesian_rank != MASTER_RANK ) 
      MPI_Send ( matrix, COLUMNS * ROWS, MPI_INT, MASTER_RANK, 6006,
		 cartesian_communicator );
    else 
      collect_matrices ( cartesian_communicator, my_cartesian_rank,
			 matrix, 6006 );

    MPI_Type_vector (ROWS, 1, COLUMNS, MPI_INT, &column_type);
    MPI_Type_commit (&column_type);

    MPI_Sendrecv ( &matrix[0][1], 1, column_type, left_neighbour, 7007,
		   &matrix[0][COLUMNS - 1], 1, column_type, 
		   right_neighbour, 7007,
		   cartesian_communicator, &status );
 
    MPI_Sendrecv ( &matrix[0][COLUMNS - 2], 1, column_type, 
		   right_neighbour, 8008,
		   &matrix[0][0], 1, column_type, left_neighbour, 8008,
		   cartesian_communicator, &status );

    if (my_cartesian_rank != MASTER_RANK )
      MPI_Send ( matrix, COLUMNS * ROWS, MPI_INT, MASTER_RANK, 9009,
		 cartesian_communicator );
    else
      collect_matrices ( cartesian_communicator, my_cartesian_rank,
			 matrix, 9009 );
  }

  MPI_Finalize ();
  exit(0);
}

int print_array (int array [DISPLAY_ROWS] [DISPLAY_COLUMNS], 
		 int vertical_break,
		 int horizontal_break)
{
  int k, l;

  printf ("\n");
  for (k = DISPLAY_ROWS - 1; k >= 0; k -- ) {
    for (l = 0; l < DISPLAY_COLUMNS; l ++ ) {
      if (l % horizontal_break == 0) printf (" ");
      printf ( "%2d ", array [k][l] );
    }
    printf ( "\n" );
    if (k % vertical_break == 0) printf ( "\n" );
  }
}

int collect_matrices (MPI_Comm cartesian_communicator, 
                      int my_cartesian_rank,
                      int matrix[ROWS][COLUMNS], 
                      int tag)
{
  int coordinates[PROCESS_DIMENSIONS];
  int client_matrix[ROWS][COLUMNS];
  int display[DISPLAY_ROWS][DISPLAY_COLUMNS];
  int i, j, k, l, source;
  MPI_Status status;
 
  for ( i = PROCESS_ROWS - 1; i >= 0; i -- ) {
    for ( j = 0; j < PROCESS_COLUMNS; j ++ ) {
      coordinates[0] = i;
      coordinates[1] = j;
      MPI_Cart_rank ( cartesian_communicator, coordinates,
		      &source );
      if (source != my_cartesian_rank) {
	MPI_Recv ( client_matrix, BUFSIZ, MPI_INT, source, tag,
		   cartesian_communicator, &status );
	for ( k = ROWS - 1; k >= 0; k -- ) {
	  for ( l = 0; l < COLUMNS; l++ ) {
	    display [i * ROWS + k] [j * COLUMNS + l] =
	      client_matrix[k][l];
	  }
	}
      }
      else {
	for ( k = ROWS - 1; k >= 0; k -- ) {
	  for ( l = 0; l < COLUMNS; l ++ ) {
	    display [i * ROWS + k] [j * COLUMNS + l]
	      = matrix[k][l];
	  }
	}
      }
    }
  }
 
  print_array (display, ROWS, COLUMNS);
}
