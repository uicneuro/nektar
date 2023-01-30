# Install script for directory: /Users/schun/nektar/library/StdRegions

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xstdregionsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/buildarm/library/StdRegions/libStdRegions.5.3.0.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.5.3.0.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.5.3.0.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.5.3.0.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xstdregionsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/buildarm/library/StdRegions/libStdRegions.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libStdRegions.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdExpansion0D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdExpansion1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdExpansion2D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdExpansion3D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdExpansion.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdHexExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdLinSysKey.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdMatrixKey.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdNodalPrismExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdNodalTetExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdNodalTriExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdPrismExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdPyrExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdQuadExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdRegionsDeclspec.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdSegExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdPointExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdTetExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdTriExp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE FILE FILES "/Users/schun/nektar/library/StdRegions/StdRegions.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/StdRegions" TYPE DIRECTORY FILES "/Users/schun/nektar/library/StdRegions/./" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\.hpp$")
endif()

