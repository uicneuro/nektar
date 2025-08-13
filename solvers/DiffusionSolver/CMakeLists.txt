CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_DIFFUSION
    "Build the Diffusion Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF (NOT NEKTAR_BUILD_SOLVERS)
    SET(NEKTAR_SOLVER_DIFFUSION OFF CACHE INTERNAL "")
ENDIF()

IF(NEKTAR_SOLVER_DIFFUSION)
    ADD_SOLVER_EXECUTABLE(DiffusionSolver SOURCES DiffusionSolver.cpp)
    ADD_SOLVER_EXECUTABLE(DiffusionSolverTimeInt SOURCES DiffusionSolverTimeInt.cpp COMPONENT diffusionsolver)

    ADD_SOLVER_EXECUTABLE(MMFNeuralEP SOURCES MMFNeuralEP.cpp
        LIBRARY_SOURCES
        NeuronModels/NeuronModel.cpp
        NeuronModels/FitzhughNagumo.cpp
        NeuronModels/HodgkinHuxley.cpp
        NeuronModels/FrankenHuxley.cpp
        NeuronModels/NeuralStimuli/NeuralStimulus.cpp
        NeuronModels/NeuralStimuli/NeuralStimulusRegion.cpp
        NeuronModels/NeuralStimuli/NeuralStimulusZone.cpp
        NeuronModels/NeuralStimuli/Protocol.cpp
        NeuronModels/NeuralStimuli/ProtocolSingle.cpp
        NeuronModels/NeuralStimuli/ProtocolS1.cpp
        NeuronModels/NeuralStimuli/ProtocolS1S2.cpp
    )

    ADD_NEKTAR_TEST(ImDiffusion_m6)
    ADD_NEKTAR_TEST(ImDiffusion_m6_time_int)

    # MMF Diffusion
    ADD_NEKTAR_TEST(TestMMFDiffPlane)
    ADD_NEKTAR_TEST(TestMMFDiffSphere)
    ADD_NEKTAR_TEST(TestMMFDiffCube)
ENDIF(NEKTAR_SOLVER_DIFFUSION)
