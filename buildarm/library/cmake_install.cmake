# Install script for directory: /Users/schun/nektar/library

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
  include("/Users/schun/nektar/buildarm/library/LibUtilities/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/LocalRegions/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/SpatialDomains/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/StdRegions/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/Collections/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/MultiRegions/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/MatrixFreeOps/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/SolverUtils/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/GlobalMapping/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/FieldUtils/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/NekMesh/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/UnitTests/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/UnitTests/SIMD/cmake_install.cmake")
  include("/Users/schun/nektar/buildarm/library/Demos/cmake_install.cmake")

endif()

