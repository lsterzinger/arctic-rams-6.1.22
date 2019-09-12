#!/bin/bash
#BSUB -P JobID
#BSUB -J Job_Name
#BSUB -o /glade/scratch/yourname/%J.out
#BSUB -e /glade/scratch/youname/%J.err
#BSUB -W 12:00
#BSUB -q economy
#BSUB -n 320
#BSUB -R span[ptile=16]
#BSUB -B
#BSUB -N
#
# Run an MPI job with the "mpirun.lsf" MPI job starter.
mpirun.lsf ./rams-6.1.12_dm
