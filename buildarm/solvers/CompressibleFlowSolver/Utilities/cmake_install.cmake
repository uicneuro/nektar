# Install script for directory: /Users/schun/nektar/solvers/CompressibleFlowSolver/Utilities

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcompressibleflow-solverx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buildarm/solvers/CompressibleFlowSolver/Utilities/ExtractSurface2DCFS")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface2DCFS")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcompressibleflow-solverx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buildarm/solvers/CompressibleFlowSolver/Utilities/ExtractSurface3DCFS")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/ExtractSurface3DCFS")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xcompressibleflow-solverx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buildarm/solvers/CompressibleFlowSolver/Utilities/CompressibleBL")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SolverUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/CompressibleBL")
    endif()
  endif()
endif()

