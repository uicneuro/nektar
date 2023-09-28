# Install script for directory: /Users/schun/nektar/library/MatrixFreeOps

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "matrixfreeops" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/buildtest/library/MatrixFreeOps/libMatrixFreeOps.5.3.0.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMatrixFreeOps.5.3.0.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMatrixFreeOps.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMatrixFreeOps.5.3.0.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMatrixFreeOps.5.3.0.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "matrixfreeops" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/buildtest/library/MatrixFreeOps/libMatrixFreeOps.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMatrixFreeOps.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMatrixFreeOps.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buildtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buildtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMatrixFreeOps.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libMatrixFreeOps.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/Operator.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTrans.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTransKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTrans.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTransKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTrans.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTransKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTrans.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTransKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTrans.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTransKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTrans.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTransKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTrans.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/BwdTransKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDeriv.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDerivKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDeriv.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDerivKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDeriv.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDerivKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDeriv.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDerivKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDeriv.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDerivKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDeriv.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDerivKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDeriv.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/PhysDerivKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/Helmholtz.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/HelmholtzKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/Helmholtz.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/HelmholtzKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/Helmholtz.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/HelmholtzKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/Helmholtz.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/HelmholtzKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/Helmholtz.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/HelmholtzKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/Helmholtz.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/HelmholtzKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProduct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProduct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProduct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProduct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProduct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProduct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProduct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBase.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBaseKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBase.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBaseKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBase.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBaseKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBase.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBaseKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBase.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBaseKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBase.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBaseKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBase.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE FILE FILES "/Users/schun/nektar/library/MatrixFreeOps/IProductWRTDerivBaseKernels.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/MatrixFreeOps" TYPE DIRECTORY FILES "/Users/schun/nektar/library/MatrixFreeOps/./" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\.hpp$")
endif()

