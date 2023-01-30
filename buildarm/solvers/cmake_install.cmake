# Install script for directory: /Users/schun/nektar/solvers

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/schun/nektar/buildarm/dist")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/schun/nektar/buildarm/solvers/ADRSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/AcousticSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/CardiacEPSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/CompressibleFlowSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/DiffusionSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/DummySolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/ImageWarpingSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/IncNavierStokesSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/LinearElasticSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/MMFSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/PulseWaveSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/ShallowWaterSolver/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/solvers/VortexWaveInteraction/cmake_install.cmake")

endif()

