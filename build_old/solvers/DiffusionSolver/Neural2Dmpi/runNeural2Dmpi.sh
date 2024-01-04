nohup mpirun -np 4 ../MMFNeuralEP NeuralEPL2NoExtmpi_v1.xml > runNeuralEPL2NoExtmpi_v1.out 2>&1 &
nohup mpirun -np 4 ../MMFNeuralEP NeuralEPL2NoExtmpi_v2.xml > runNeuralEPL2NoExtmpi_v2.out 2>&1 &
nohup mpirun -np 4 ../MMFNeuralEP NeuralEPL2NoExtmpi_v3.xml > runNeuralEPL2NoExtmpi_v3.out 2>&1 &

nohup mpirun -np 4 ../MMFNeuralEP NeuralEPL2mpi_v1.xml > runNeuralEPL2mpi_v1.out 2>&1 &
nohup mpirun -np 4 ../MMFNeuralEP NeuralEPL2mpi_v2.xml > runNeuralEPL2mpi_v2.out 2>&1 &
nohup mpirun -np 4 ../MMFNeuralEP NeuralEPL2mpi_v3.xml > runNeuralEPL2mpi_v3.out 2>&1 &
