#!/bin/bash
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ___________________________________________________________________
#  |                                                                   | 
#  | Set initial information for the Queuing system                    | 
#  | ==============================================                    | 
#  |                                                                   | 
#  | All PBS directives (the lines starting with #PBS) below           |
#  | are optional and can be omitted, resulting in the use of          | 
#  | the system defaults.                                              | 
#  |___________________________________________________________________|
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#PBS -l select=2:ncpus=32:mpiprocs=16
#
#    2 nodes and 16 MPI processes per node. Each node has 16 cores, but
#    with Intel hyper-threading 2 threads can run on each core, so PBS 
#    therefore see a total of 32 cores per node. You must always specify
#    ncpus=32, so that you reserve the entire nodes that you use.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#PBS -lwalltime=00:05:00
#
#    Expecting to run no more than 5 minutes.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
###PBS -m abe
#
#    The queuing system will send an email on job Abort, Begin, End
#    Remember to specify a .forward file in your home directory with
#    your actual e-mail address. Otherwise use the following option:
###PBS -M user@example.com
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#PBS -o output
#PBS -e output
#
#    Standard output and error paths respectively.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#PBS -N cpi
#
#    Name of the process in the queue.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#PBS -A ACCOUNT_MISSING
#
#    Account to use - one must be provided.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#   ______________________________________________________________
#  |                                                              | 
#  |                                                              |
#  |  Setting up and running your job on vilje                    |
#  |  ========================================                    |
#  |                                                              |
#  |  We are now ready to begin running commands.                 |
#  |                                                              |
#  |  This job script is run as a regular shell script on the     |
#  |  first node assigned to this job, hereafter called Mother    |
#  |  Superior.                                                   |
#  |______________________________________________________________|

module load mpt

cd $PBS_O_WORKDIR

mpiexec_mpt bin/cpi

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
#    Running your MPI job by launching it with "mpirun". This works for
#    when using OpenMPI. I you are using some other MPI distribution you may
#    have to use a different launcher.
#
#    Don't put the job commands in the background, e.g. by adding & at
#    the end.  Doing so will make the job escape the queuing system, 
#    create havoc with the scheduling and result in you getting angry
#    email from the sysadmins.
#
#    The job ends when the script exit.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

# This patchscript is copied from stallo.uit.no and modified for use on vilje.hpc.ntnu.no
