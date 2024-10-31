#!/bin/bash
# The job name
#SBATCH --job-name=dwarf-p-cloudsc-tests
# Set the error and output files to [job name]-[JobID].out
# Note that some cloudsc output is sent in the error stream, so we use the same file for 
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.out
# Set the initial working directory
#SBATCH --chdir=/scratch/user
# Choose the queue
#SBATCH --qos=dg
# Wall clock time limit
#SBATCH --time=00:05:00
# Send an email on failure
#SBATCH --mail-type=FAIL
# This is the job

