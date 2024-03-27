# Install script for directory: /Users/schun/nektar/library/Demos/MultiRegions

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

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Helmholtz1D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz1D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz1D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz1D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz1D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Helmholtz2D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz2D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Helmholtz3D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Helmholtz3DHomo1D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo1D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo1D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo1D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo1D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Helmholtz3DHomo2D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo2D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo2D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo2D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Helmholtz3DHomo2D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/HDGHelmholtz1D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz1D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz1D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz1D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz1D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/HDGHelmholtz2D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz2D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz2D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz2D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz2D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/HDGHelmholtz3D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/HDGHelmholtz3DHomo1D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3DHomo1D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3DHomo1D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3DHomo1D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/HDGHelmholtz3DHomo1D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/PostProcHDG2D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG2D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG2D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG2D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG2D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/PostProcHDG3D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG3D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG3D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG3D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/PostProcHDG3D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Deriv3DHomo1D_SingleMode")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D_SingleMode" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D_SingleMode")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D_SingleMode")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D_SingleMode")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Deriv3DHomo1D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo1D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Deriv3DHomo2D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo2D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo2D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo2D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Deriv3DHomo2D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/Int3DHomo1D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Int3DHomo1D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Int3DHomo1D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Int3DHomo1D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/Int3DHomo1D")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "demos" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/bin" TYPE EXECUTABLE OPTIONAL FILES "/Users/schun/nektar/build/library/Demos/MultiRegions/SteadyAdvectionDiffusionReaction2D")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/SteadyAdvectionDiffusionReaction2D" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/SteadyAdvectionDiffusionReaction2D")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/SteadyAdvectionDiffusionReaction2D")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -u -r "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/bin/SteadyAdvectionDiffusionReaction2D")
    endif()
  endif()
endif()

