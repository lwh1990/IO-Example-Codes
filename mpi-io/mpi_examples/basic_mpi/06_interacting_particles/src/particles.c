/******************************************************************************
 *                                                                            *
 *  Basic MPI Example - Interacting Particles                                 *
 *                                                                            *
 *  Simulate the interaction between 10k electrically charged particles.      *
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
#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#define FALSE            0
#define TRUE             1
#define MASTER_RANK      0

#define MAX_PARTICLES 10000
#define MAX_PROCS       128
#define EPSILON           1.0E-10
#define DT               0.01
#define N_OF_ITERATIONS 20


int main ( int argc, char **argv )
{
  int pool_size, my_rank;
  int i_am_the_master = FALSE; 
  extern double drand48();
  extern void srand48();

  typedef struct {
    double x, y, z, vx, vy, vz, ax, ay, az, mass, charge;
  } Particle;

  Particle  particles[MAX_PARTICLES];  /* Particles on all nodes */
  int       counts[MAX_PROCS];         /* Number of ptcls on each proc */
  int       displacements[MAX_PROCS];  /* Offsets into particles */
  int       offsets[MAX_PROCS];        /* Offsets used by the master */
  int       particle_number, i, j, my_offset, true_i;
  int       total_particles;           /* Total number of particles */
  int       count;                     /* Count time steps */

  MPI_Datatype particle_type;

  double    dt = DT;                   /* Integration time step */
  double    comm_time, start_comm_time, end_comm_time, start_time, end_time;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &pool_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank == MASTER_RANK) i_am_the_master = TRUE;

  particle_number = MAX_PARTICLES / pool_size;

  if (i_am_the_master)
    printf ("%d particles per processor\n", particle_number);

  MPI_Type_contiguous ( 11, MPI_DOUBLE, &particle_type );
  MPI_Type_commit ( &particle_type );

  MPI_Allgather ( &particle_number, 1, MPI_INT, counts, 1, MPI_INT,
                  MPI_COMM_WORLD );

  displacements[0] = 0;
  for (i = 1; i < pool_size; i++)
    displacements[i] = displacements[i-1] + counts[i-1];
  total_particles = displacements[pool_size - 1]
    + counts[pool_size - 1];
       
  if (i_am_the_master)
    printf ("total number of particles = %d\n", total_particles);

  my_offset = displacements[my_rank];

  MPI_Gather ( &my_offset, 1, MPI_INT, offsets, 1, MPI_INT, MASTER_RANK,
               MPI_COMM_WORLD );

  if (i_am_the_master) {
    printf ("offsets: ");
    for (i = 0; i < pool_size; i++)
      printf ("%d ", offsets[i]);
    printf("\n");
  }

  srand48((long) (my_rank + 28));

  /* Here each process initializes its own particles. */

  for (i = 0; i < particle_number; i++) {
    particles[my_offset + i].x = drand48();
    particles[my_offset + i].y = drand48();
    particles[my_offset + i].z = drand48();
    particles[my_offset + i].vx = 0.0;
    particles[my_offset + i].vy = 0.0;
    particles[my_offset + i].vz = 0.0;
    particles[my_offset + i].ax = 0.0;
    particles[my_offset + i].ay = 0.0;
    particles[my_offset + i].az = 0.0;
    particles[my_offset + i].mass = 1.0;
    particles[my_offset + i].charge = 1.0 - 2.0 * (i % 2);
  }

  start_time = MPI_Wtime();
  comm_time = 0.0;

  for (count = 0; count < N_OF_ITERATIONS; count++) {

    if (i_am_the_master) printf("Iteration %d.\n", count + 1);

    /* Here processes exchange their particles with each other. */

    start_comm_time = MPI_Wtime();

    MPI_Allgatherv ( particles + my_offset, particle_number,
                     particle_type,
                     particles, counts, displacements, particle_type,
                     MPI_COMM_WORLD );

    end_comm_time = MPI_Wtime();
    comm_time += end_comm_time - start_comm_time;

    for (i = 0; i < particle_number; i++) {

      true_i = i + my_offset;

      /* initialize accelerations to zero */

      particles[true_i].ax = 0.0;
      particles[true_i].ay = 0.0;
      particles[true_i].az = 0.0;

      for (j = 0; j < total_particles; j++) {

        /* Do not evaluate interaction with yourself. */

        if (j != true_i) {

          /* Evaluate forces that j-particles exert on the i-particle. */

          double dx, dy, dz, r2, r, qj_by_r3;

          /* Here we absorb the minus sign by changing the order
             of i and j. */

          dx = particles[true_i].x - particles[j].x;
          dy = particles[true_i].y - particles[j].y;
          dz = particles[true_i].z - particles[j].z;

          r2 = dx * dx + dy * dy + dz * dz; r = sqrt(r2);

          /* Quench the force if the particles are too close. */

          if (r < EPSILON) qj_by_r3 = 0.0;
          else qj_by_r3 = particles[j].charge / (r2 * r);

          /* accumulate the contribution from particle j */

          particles[true_i].ax += qj_by_r3 * dx;
          particles[true_i].ay += qj_by_r3 * dy;
          particles[true_i].az += qj_by_r3 * dz;
        }
      }
    }

    /*
     * We advance particle positions and velocities only *after* 
     * we have evaluated all accelerations using the *old* positions.
     */

    for (i = 0; i < particle_number; i++) {

      double qidt_by_m, dt_by_2, vx0, vy0, vz0;

      true_i = i + my_offset;

      /* Save old velocities */

      vx0 = particles[true_i].vx;
      vy0 = particles[true_i].vy;
      vz0 = particles[true_i].vz;

      /* Now advance the velocity of particle i */

      qidt_by_m = particles[true_i].charge * dt / particles[true_i].mass;
      particles[true_i].vx += particles[true_i].ax * qidt_by_m;
      particles[true_i].vy += particles[true_i].ay * qidt_by_m;
      particles[true_i].vz += particles[true_i].az * qidt_by_m;

      /* Use average velocity in the interval to advance the particles' 
         positions */

      dt_by_2 = 0.5 * dt; 
      particles[true_i].x += (vx0 + particles[true_i].vx) * dt_by_2;
      particles[true_i].y += (vy0 + particles[true_i].vy) * dt_by_2;
      particles[true_i].z += (vz0 + particles[true_i].vz) * dt_by_2;
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  end_time = MPI_Wtime();

  if (i_am_the_master) {
    printf ("Communication time %8.5f seconds\n", comm_time);
    printf ("Computation time   %8.5f seconds\n", 
            end_time - start_time - comm_time);
    printf ("\tEvaluated %d interactions\n",
            N_OF_ITERATIONS * total_particles * (total_particles - 1));
  }

  MPI_Finalize ();

  exit(0);
}
