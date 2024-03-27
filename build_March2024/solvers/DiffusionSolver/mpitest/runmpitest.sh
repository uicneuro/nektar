

nohup ../MMFNeuralEP NeuralEP2DOneL2mono.xml > runNeuralEP2DOneL2mono.out 2>&1 &
nohup ../MMFNeuralEP NeuralEP2DbiL2.xml > runNeuralEP2DbiL2.out 2>&1 &
nohup ../MMFNeuralEP NeuralEP2DbiL2NoExt.xml > runNeuralEP2DbiL2NoExt.out 2>&1 &

echo "Start sleep"
sleep 2
echo "End sleep"

nohup mpirun -np 4 ../MMFNeuralEP NeuralEP2DOneL2monompi4.xml > runNeuralEP2DOneL2monompi4.out 2>&1 &
nohup mpirun -np 4 ../MMFNeuralEP NeuralEP2DbiL2mpi4.xml > runNeuralEP2DbiL2mpi4.out 2>&1 &
nohup mpirun -np 4 ../MMFNeuralEP NeuralEP2DbiL2NoExtmpi4.xml > runNeuralEP2DbiL2NoExtmpi4.out 2>&1 &
