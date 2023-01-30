# Install script for directory: /Users/schun/nektar/solvers/CompressibleFlowSolver

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/schun/nektar/build/dist")
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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcompressibleflow-solverx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/solvers/CompressibleFlowSolver/CompressibleFlowSolver")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleFlowSolver")
    endif()
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/schun/nektar/build/solvers/CompressibleFlowSolver/RiemannSolvers/UnitTests/cmake_install.cmake")
  include("/Users/schun/nektar/build/solvers/CompressibleFlowSolver/Utilities/cmake_install.cmake")

endif()

