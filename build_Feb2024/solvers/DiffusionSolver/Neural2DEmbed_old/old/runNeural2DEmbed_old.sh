

nohup mpirun -np 8 ../MMFNeuralEP NeuralEP2DmonoDL6mpi8.xml > runNeuralEP2DmonoDL6mpi8.out 2>&1 &
nohup mpirun -np 16 ../MMFNeuralEP NeuralEP2DmonoDL6mpi16.xml > runNeuralEP2DmonoDL6mpi16.out 2>&1 &

nohup mpirun -np 8 ../MMFNeuralEP NeuralEP2DmonoDNoPhieL6mpi8.xml > runNeuralEP2DmonoDNoPhieL6mpi8.out 2>&1 &
nohup mpirun -np 16 ../MMFNeuralEP NeuralEP2DmonoDNoPhieL6mpi16.xml > runNeuralEP2DmonoDNoPhieL6mpi16.out 2>&1 &

nohup ../MMFNeuralEP NeuralEP2DmonoDWCL2_D.xml > runNeuralEP2DmonoDWCL2_D.out 2>&1 &

nohup ../MMFNeuralEP NeuralEP2DmonoDWCL2_Np3.xml > runNeuralEP2DmonoDWCL2_Np3.out 2>&1 &
nohup ../MMFNeuralEP NeuralEP2DmonoDWCL2_Np4.xml > runNeuralEP2DmonoDWCL2_Np4.out 2>&1 &
nohup ../MMFNeuralEP NeuralEP2DmonoDWCL2_Np5.xml > runNeuralEP2DmonoDWCL2_Np5.out 2>&1 &
nohup ../MMFNeuralEP NeuralEP2DmonoDWCL2_Np6.xml > runNeuralEP2DmonoDWCL2_Np6.out 2>&1 &
nohup ../MMFNeuralEP NeuralEP2DmonoDWCL2_Np7.xml > runNeuralEP2DmonoDWCL2_Np7.out 2>&1 &
