/******************************************************************************
 *                                                                            *
 *  Basic MPI Example - Random Pie                                            *
 *                                                                            *
 *  Approximating pi using a Monte Carlo approach.                            *
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
#include <stdlib.h> /* function random is defined here */
#include <limits.h> /* INT_MAX is defined here */
#include <mpi.h>
#include <math.h>

#define CHUNKSIZE 1000
#define REQUEST      1
#define REPLY        2
#define TOTAL     5000000

main(argc, argv)
     int argc;
     char *argv[];
{
  int iter;
  int in, out, i, iters, max, ix, iy, ranks[1], done, temp;
  double x, y, Pi, error, epsilon;
  int numprocs, myid, server, totalin, totalout, workerid;
  int rands[CHUNKSIZE], request;
  MPI_Comm world, workers;
  MPI_Group world_group, worker_group;
  MPI_Status stat;

  MPI_Init(&argc, &argv);
  world = MPI_COMM_WORLD;
  MPI_Comm_size(world, &numprocs);
  MPI_Comm_rank(world, &myid);
  server = numprocs - 1;
  if (myid == 0)
    sscanf( argv[1], "%lf", &epsilon);
  MPI_Bcast(&epsilon, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Comm_group( world, &world_group);
  ranks[0] = server;
  MPI_Group_excl(world_group, 1, ranks, &worker_group);
  MPI_Comm_create(world, worker_group, &workers);
  MPI_Group_free(&worker_group);

  if(myid == server) {
    do {
      MPI_Recv(&request, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST, world, &stat);
      if (request) {
        for (i = 0; i < CHUNKSIZE; i++)
          rands[i] = random();
        MPI_Send(rands, CHUNKSIZE, MPI_INT, stat.MPI_SOURCE, REPLY, world);
      }
    }
    while (request> 0);
  }
  else {
    request = 1;
    done = in = out = 0;
    max = INT_MAX;
    MPI_Send(&request, 1, MPI_INT, server, REQUEST, world);
    MPI_Comm_rank(workers, &workerid);
    iter = 0;
    while(!done) {
      iter++;
      request = 1;
      MPI_Recv(rands, CHUNKSIZE, MPI_INT, server, REPLY, world, &stat);
      for (i=0; i < CHUNKSIZE; ) {
        x = (((double) rands[i++])/max) * 2 - 1;
        y = (((double) rands[i++])/max) * 2 - 1;
        if (x*x + y*y < 1.0)
          in++;
        else
          out++;
      }
      MPI_Allreduce(&in, &totalin, 1, MPI_INT, MPI_SUM, workers);
      MPI_Allreduce(&out, &totalout, 1, MPI_INT, MPI_SUM, workers);
      Pi = (4.0*totalin)/(totalin + totalout);
      error = fabs( Pi-3.141592653589793238462643);
      done = ((error < epsilon) || ((totalin+totalout) > TOTAL));
      request = (done) ? 0 : 1;
      if (myid == 0) {
        printf( "\rpi = %23.20lf", Pi);
        MPI_Send(&request, 1, MPI_INT, server, REQUEST, world);
      }
      else {
        if (request)
          MPI_Send(&request, 1, MPI_INT, server, REQUEST, world);
      }
    }
    MPI_Comm_free(&workers);
  }
  if (myid == 0) 
    printf( "\npoints: %d\nin: %d, out: %d\n", 
            totalin+totalout, totalin, totalout);
  MPI_Finalize();
  exit (0);
}
