/******************************************************************************
 *                                                                            *
 *  Basic MPI Example - Handling Errors                                       *
 *                                                                            *
 *  Various mechanisms to handle errors.                                      *
 *                                                                            *
 ******************************************************************************
 *                                                                            *
 *  The original code was written by Gustav at University of Indiana in 2003. *
 *                                                                            *
 *  The current version has been tested/updated by the HPC department at      *
 *  the Norwegian University of Science and Technology in 2011.               *
 *                                                                            *
 ******************************************************************************/
#include <stdio.h>   /* all IO stuff lives here */
#include <stdlib.h>  /* exit lives here */
#include <unistd.h>  /* getopt lives here */
#include <errno.h>   /* UNIX error handling lives here */
#include <string.h>  /* strcpy lives here */
#include <mpi.h>     /* MPI and MPI-IO live here */

#define MASTER_RANK 0
#define TRUE 1
#define FALSE 0
#define BOOLEAN int
#define BLOCK_SIZE 1048576
#define SYNOPSIS printf ("synopsis: %s -f <file> -l <blocks>\n", argv[0])

int main(argc, argv)
     int argc;
     char *argv[];
{
  /* my variables */

  int my_rank, pool_size, number_of_blocks = 0, block, i;
  BOOLEAN i_am_the_master = FALSE, input_error = FALSE, 
    my_file_open_error = FALSE, file_open_error = FALSE,
    my_write_error = FALSE, write_error = FALSE;
  char *basename = NULL, file_name[BUFSIZ], message[BUFSIZ];
  int basename_length, junk[BLOCK_SIZE];
  FILE *fp;
  double start, finish, io_time = 0.0;

  /* getopt variables */

  extern char *optarg;
  int c;

  /* error handling variables */

  extern int errno;

  /* ACTION */

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
  if (my_rank == MASTER_RANK) i_am_the_master = TRUE;

  if (i_am_the_master) {

    /* read the command line */

    while ((c = getopt(argc, argv, "f:l:h")) != EOF)
      switch(c) {
      case 'f': 
        basename = optarg;
	break;
      case 'l': 
	if ((sscanf (optarg, "%d", &number_of_blocks) != 1) ||
	    (number_of_blocks < 1)) 
	  input_error = TRUE;
	break;
      case 'h':
	input_error = TRUE;
	break;
      case '?':
	input_error = TRUE;
	break;
      }

    /* Check if the command line has initialized basename and
     * number_of_blocks.
     */

    if ((basename == NULL) || (number_of_blocks == 0)) input_error = TRUE;

    if (input_error)
       SYNOPSIS;
    else {
      basename_length = strlen(basename) + 1;
#ifdef DEBUG
      printf("basename         = %s\n", basename);
      printf("basename_length  = %d\n", basename_length);
      printf("number_of_blocks = %d\n", number_of_blocks);
#endif
    }
  }

  /* Transmit the effect of reading the command line to other
     processes. */

  MPI_Bcast(&input_error, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);

  if (! input_error) {
    MPI_Bcast(&number_of_blocks, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    MPI_Bcast(&basename_length, 1, MPI_INT, MASTER_RANK, MPI_COMM_WORLD);
    if (! i_am_the_master) basename = (char*) malloc(basename_length);
    MPI_Bcast(basename, basename_length, MPI_CHAR, MASTER_RANK, MPI_COMM_WORLD);

#ifdef DEBUG
    printf("%3d: basename = %s, number_of_blocks = %d\n", 
	   my_rank, basename, number_of_blocks);
#endif

    /* Now every process creates its own file name and attempts
       to open the file. */

    sprintf(file_name, "%s.%d", basename, my_rank);

#ifdef DEBUG
    printf("%3d: opening file %s\n", my_rank, file_name);
#endif

    if (! (fp = fopen(file_name, "w"))) {
      sprintf(message, "%3d: %s", my_rank, file_name);
      perror(message);
      my_file_open_error = TRUE;
    }

    /* Now we must ALL check that NOBODY had problems
       with opening the file. */

    MPI_Allreduce (&my_file_open_error, &file_open_error, 1, MPI_INT, 
		   MPI_LOR, MPI_COMM_WORLD);

#ifdef DEBUG
    if (i_am_the_master)
      if (file_open_error)
	fprintf(stderr, "problem opening output files\n");
#endif

    /* If all files are open for writing, write to them */

    if (! file_open_error) {
      srand(28 + my_rank);
      for (block = 0; (block < number_of_blocks) || my_write_error; 
	   block++) {
	for (i = 0; i < BLOCK_SIZE; junk[i++] = rand());
        start = MPI_Wtime();
	if (fwrite(junk, sizeof(int), BLOCK_SIZE, fp) != BLOCK_SIZE) {
	  sprintf(message, "%3d: %s", my_rank, file_name);
	  perror(message);
	  my_write_error = TRUE;
	}
        finish = MPI_Wtime();
        io_time += finish - start;
      }
        
      /* Check if anybody had problems writing on the file */

      MPI_Allreduce (&my_write_error, &write_error, 1, MPI_INT,
		     MPI_LOR, MPI_COMM_WORLD);

#ifdef DEBUG
      if (i_am_the_master)
	if (write_error)
	  fprintf(stderr, "problem writing on files\n");
#endif
      if (i_am_the_master)
        if (!write_error)
          printf("io_time = %f\n", io_time);

    }

    /* Only processes that were successful opening the files
       need do close them here */

    if (!my_file_open_error) {
      fclose(fp);
#ifdef DEBUG
      printf ("%3d: closed %s\n", my_rank, file_name);
#endif
    }

    /* If we have either write errors or file open errors,
       then processes that managed to open their files
       are requested to throw them away */

    if ((write_error || file_open_error) && !my_file_open_error) {
      unlink(file_name);
#ifdef DEBUG
      printf("%3d: unlinked %s\n", my_rank, file_name);
#endif
    }

    /* We don't try to capture unlink or fclose errors here,
       because there is little we could do about them. */
       
  } 

  MPI_Finalize();
  exit(0);
}
