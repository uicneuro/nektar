nohup mpirun -np 8 ./MMFNeuralEP FiberLineL20Duo.xml > runFiberLineL20Duo.out 2>&1 &
nohup mpirun -np 8 ./MMFNeuralEP FiberLineL20DuoIsol.xml > runFiberLineL20DuoIsol.out 2>&1 &
nohup mpirun -np 8 ./MMFNeuralEP FiberLineL20Single.xml > runFiberLineL20Single.out 2>&1 &
nohup mpirun -np 8 ./MMFNeuralEP FiberLineL20SingleIsol.xml > runFiberLineL20SingleIsol.out 2>&1 &

nohup ./MMFNeuralEP FiberLineL20DuoC1.xml > runFiberLineL20DuoC1.out 2>&1 &
nohup ./MMFNeuralEP FiberLineL20SingleC1.xml > runFiberLineL20SingleC1.out 2>&1 &