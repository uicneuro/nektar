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

    IF(NEKTAR_BUILD_SOLVER_LIBS)
        MESSAGE(WARNING, "Solver library build not set up for MMFSolver - skipping.")
    ENDIF()

    ADD_NEKTAR_TEST        (TestMaxwellSphere)
ENDIF( NEKTAR_SOLVER_MMF )
