CMAKE_DEPENDENT_OPTION(NEKTAR_SOLVER_MMF
    "Build the ADR Solver." ON
    "NEKTAR_BUILD_SOLVERS" OFF)

IF(NEKTAR_SOLVER_MMF)
    ADD_SOLVER_EXECUTABLE(MMFSolver SOURCES MMFSolver.cpp
        LIBRARY_SOURCES
        ./EquationSystems/MMFDiffusion.cpp
        ./EquationSystems/MMFMaxwell.cpp
    )

    ADD_SOLVER_EXECUTABLE(MMFCardiacEP SOURCES MMFCardiacEP.cpp
        LIBRARY_SOURCES
        ../CardiacEPSolver/CellModels/CellModel.cpp
        ../CardiacEPSolver/CellModels/FitzhughNagumo.cpp
        ../CardiacEPSolver/CellModels/AlievPanfilov.cpp
        ../CardiacEPSolver/CellModels/CourtemancheRamirezNattel98.cpp
        ../CardiacEPSolver/CellModels/FentonKarma.cpp
        ../CardiacEPSolver/CellModels/LuoRudy91.cpp
        ../CardiacEPSolver/CellModels/Fox02.cpp
        ../CardiacEPSolver/CellModels/Winslow99.cpp
        ../CardiacEPSolver/CellModels/TenTusscher06.cpp
        ../CardiacEPSolver/Stimuli/Stimulus.cpp
        ../CardiacEPSolver/Stimuli/StimulusCircle.cpp
        ../CardiacEPSolver/Stimuli/StimulusRect.cpp
        ../CardiacEPSolver/Stimuli/StimulusPoint.cpp
        ../CardiacEPSolver/Stimuli/Protocol.cpp
        ../CardiacEPSolver/Stimuli/ProtocolSingle.cpp
        ../CardiacEPSolver/Stimuli/ProtocolS1.cpp
        ../CardiacEPSolver/Stimuli/ProtocolS1S2.cpp
    )

    ADD_SOLVER_EXECUTABLE(MMFNeuralEP SOURCES MMFNeuralEP.cpp
        LIBRARY_SOURCES
        ./EquationSystems/NeuronModels/NeuronModel.cpp
        ./EquationSystems/NeuronModels/NeuronModelmulti.cpp
        ./EquationSystems/NeuronModels/FitzhughNagumo.cpp
        ./EquationSystems/NeuronModels/HodgkinHuxley.cpp
        ./EquationSystems/NeuronModels/FrankenHuxley.cpp
        ./EquationSystems/NeuronModels/NeuralStimuli/NeuralStimulus.cpp
        ./EquationSystems/NeuronModels/NeuralStimuli/NeuralStimulusRegion.cpp
        ./EquationSystems/NeuronModels/NeuralStimuli/NeuralStimulusZone.cpp
        ./EquationSystems/NeuronModels/NeuralStimuli/Protocol.cpp
        ./EquationSystems/NeuronModels/NeuralStimuli/ProtocolSingle.cpp
        ./EquationSystems/NeuronModels/NeuralStimuli/ProtocolS1.cpp
        ./EquationSystems/NeuronModels/NeuralStimuli/ProtocolS1S2.cpp
    )

    IF(NEKTAR_BUILD_SOLVER_LIBS)
        MESSAGE(WARNING, "Solver library build not set up for MMFSolver - skipping.")
    ENDIF()

    ADD_NEKTAR_TEST        (TestMaxwellSphere)
ENDIF( NEKTAR_SOLVER_MMF )
