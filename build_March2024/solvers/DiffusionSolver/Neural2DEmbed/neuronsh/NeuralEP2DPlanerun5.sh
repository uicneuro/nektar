#!/bin/sh
#PBS -V
#PBS -N N_run5
#PBS -q long
#PBS -A etc
#PBS -l select=6:ncpus=64:mpiprocs=64
#PBS -l walltime=120:00:00

cd /scratch/e1460a01/nektar/build/solvers/DiffusionSolver/Neural2DEmbed

module purge
module load craype-mic-knl gcc/7.2.0 openmpi/3.1.0

mpirun ../MMFNeuralEP Fiber2DEphapticPlaneL6rho4.xml > runFiber2DEphapticPlaneL6rho4.out