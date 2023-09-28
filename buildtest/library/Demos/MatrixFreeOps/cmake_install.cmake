# Install script for directory: /Users/schun/nektar/library/Demos/MatrixFreeOps

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/schun/nektar/buildtest/dist")
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

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buildtest/library/Demos/MatrixFreeOps/MFO_Timing")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MFO_Timing" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MFO_Timing")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buildtest/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/Collections"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/buildtest/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MFO_Timing")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MFO_Timing")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buildtest/library/Demos/MatrixFreeOps/MFO_Timing_3D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MFO_Timing_3D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MFO_Timing_3D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buildtest/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/Collections"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/buildtest/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MFO_Timing_3D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/MFO_Timing_3D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buildtest/library/Demos/MatrixFreeOps/Helmholtz2D_MFO")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D_MFO" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D_MFO")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buildtest/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/Collections"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/buildtest/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D_MFO")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D_MFO")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buildtest/library/Demos/MatrixFreeOps/Helmholtz3D_MFO")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D_MFO" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D_MFO")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buildtest/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/Collections"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/buildtest/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/buildtest/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D_MFO")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D_MFO")
    endif()
  endif()
endif()

