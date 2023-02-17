nohup ./MMFNeuralEP TestMMFNeuralEP2Dbi.xml > runTestMMFNeuralEP2Dbi.out 2>&1 &
nohup mpirun -np 10 ./MMFNeuralEP TestMMFNeuralEP2Dbi.xml > runTestMMFNeuralEP2Dbimpi.out 2>&1 &