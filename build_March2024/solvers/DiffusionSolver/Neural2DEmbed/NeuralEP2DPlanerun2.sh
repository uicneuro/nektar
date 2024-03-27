#!/bin/sh
#PBS -V
#PBS -N N_run2
#PBS -q normal
#PBS -A etc
#PBS -l select=8:ncpus=48:mpiprocs=48
#PBS -l walltime=48:00:00

cd /scratch/e1460a01/nektar/build/solvers/DiffusionSolver/Neural2DEmbed

module purge
module load craype-mic-knl gcc/7.2.0 openmpi/3.1.0

mpirun ../MMFNeuralEP Fiber2DEphapticPlaneL6rho1.xml > runFiber2DEphapticPlaneL6rho1.out