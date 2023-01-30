# Install script for directory: /Users/schun/nektar/library/MultiRegions

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmultiregionsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/build/library/MultiRegions/libMultiRegions.5.3.0.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.5.3.0.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xmultiregionsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/build/library/MultiRegions/libMultiRegions.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMultiRegions.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ContField.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ContField3DHomogeneous1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ContField3DHomogeneous2D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/DisContField.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/DisContField3DHomogeneous1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/DisContField3DHomogeneous2D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ExpList.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ExpListHomogeneous1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ExpListHomogeneous2D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ExpList2DHomogeneous1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ExpList1DHomogeneous2D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ExpList3DHomogeneous1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/ExpList3DHomogeneous2D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GJPStabilisation.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSys.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSysKey.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSysDirect.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSysDirectFull.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSysDirectStaticCond.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSysIterative.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSysIterativeFull.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSysIterativeStaticCond.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalLinSysStaticCond.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalMatrix.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/GlobalMatrixKey.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/MultiRegions.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/MultiRegionsDeclspec.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/Preconditioner.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/PreconditionerDiagonal.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/PreconditionerLowEnergy.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/PreconditionerBlock.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE FILE FILES "/Users/schun/nektar/library/MultiRegions/SubStructuredGraph.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MultiRegions" TYPE DIRECTORY FILES "/Users/schun/nektar/library/MultiRegions/./" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\.hpp$")
endif()

