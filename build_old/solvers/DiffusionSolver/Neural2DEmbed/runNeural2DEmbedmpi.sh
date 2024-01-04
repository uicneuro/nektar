

nohup mpirun -np 4 ../MMFNeuralEP NeuralEP2DExtL2NoExtmpi.xml > runNeuralEP2DExtL2NoExtmpi.out 2>&1 &
nohup mpirun -np 4 ../MMFNeuralEP NeuralEP2DExtL2Phiempi.xml > runNeuralEP2DExtL2Phiempi.out 2>&1 &
nohup mpirun -np 4 ../MMFNeuralEP NeuralEP2DExtL2Approxmpi.xml > runNeuralEP2DExtL2Approxmpi.out 2>&1 &

