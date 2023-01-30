# Install script for directory: /Users/schun/nektar/library/UnitTests/LocalRegions

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xunit-testx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buildarm/library/UnitTests/LocalRegions/LocalRegionsUnitTests")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/LocalRegionsUnitTests")
    endif()
  endif()
endif()

