# Install script for directory: /Users/schun/nektar/utilities/FieldConvert

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/schun/nektar/buidtest/dist")
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

if(CMAKE_INSTALL_COMPONENT STREQUAL "fieldconvert" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/buidtest/utilities/FieldConvert/FieldConvert")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FieldConvert" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FieldConvert")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buidtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buidtest/library/FieldUtils"
      -delete_rpath "/Users/schun/nektar/buidtest/library/GlobalMapping"
      -delete_rpath "/Users/schun/nektar/buidtest/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/buidtest/library/Collections"
      -delete_rpath "/Users/schun/nektar/buidtest/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/buidtest/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/buidtest/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/buidtest/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/buidtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buidtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buidtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FieldConvert")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/FieldConvert")
    endif()
  endif()
endif()

