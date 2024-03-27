#!/bin/sh
#PBS -V
#PBS -N NeuralEP2DPlane
#PBS -q normal
#PBS -A etc
#PBS -l select=2:ncpus=16:mpiprocs=16
#PBS -l walltime=48:00:00

cd /scratch/e1460a01/nektar/build/solvers/DiffusionSolver/Neural2DEmbed

module purge
module load craype-mic-knl gcc/7.2.0 openmpi/3.1.0

mpirun ../MMFNeuralEP Fiber2DEphapticPlaneL6.xml > runFiber2DEphapticPlaneL6.out