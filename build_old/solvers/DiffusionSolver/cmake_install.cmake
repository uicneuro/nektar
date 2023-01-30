# Install script for directory: /Users/schun/nektar/solvers/DiffusionSolver

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdiffusion-solverx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/solvers/DiffusionSolver/DiffusionSolver")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolver")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdiffusionsolverx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/solvers/DiffusionSolver/DiffusionSolverTimeInt")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/DiffusionSolverTimeInt")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmmfdiffusionx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/solvers/DiffusionSolver/MMFDiffusion")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFDiffusion")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmmfcardiacepx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/solvers/DiffusionSolver/MMFCardiacEP")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFCardiacEP")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmmfneuralepx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/solvers/DiffusionSolver/MMFNeuralEP")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MMFNeuralEP")
    endif()
  endif()
endif()

