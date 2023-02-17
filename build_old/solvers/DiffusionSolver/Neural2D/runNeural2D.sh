nohup mpirun -np 8 ../MMFNeuralEP Fiber2DL06bimpi.xml > runFiber2DL06bimpi.out 2>&1 &
nohup ../MMFNeuralEP Fiber2DL06bi.xml > runFiber2DL06bi.out 2>&1 &

nohup ../MMFNeuralEP TestMMFNeuralEP2DbiS10.xml > runTestMMFNeuralEP2DbiS10.out 2>&1 &
nohup ../MMFNeuralEP TestMMFNeuralEP2DbiS20.xml > runTestMMFNeuralEP2DbiS20.out 2>&1 &
nohup ../MMFNeuralEP TestMMFNeuralEP2DbiS50.xml > runTestMMFNeuralEP2DbiS50.out 2>&1 &
