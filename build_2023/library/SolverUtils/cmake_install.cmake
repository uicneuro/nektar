# Install script for directory: /Users/schun/nektar/library/SolverUtils

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xsolverutilsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/build/library/SolverUtils/libSolverUtils.5.3.0.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.5.3.0.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xsolverutilsx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/build/library/SolverUtils/libSolverUtils.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/FieldUtils"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/GlobalMapping"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MultiRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/Collections"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/MatrixFreeOps"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/build/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/build/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libSolverUtils.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Core" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Core/Coupling.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Core" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Core/CouplingFile.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Core" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Core/Misc.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Core" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Core/SessionFunction.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/AdvectionSystem.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Advection" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Advection/Advection.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Advection" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Advection/AdvectionFR.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Advection" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Advection/Advection3DHomogeneous1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Advection" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Advection/AdvectionNonConservative.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Advection" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Advection/AdvectionWeakDG.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Diffusion" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Diffusion/Diffusion.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Diffusion" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Diffusion/Diffusion3DHomogeneous1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Diffusion" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Diffusion/DiffusionLDG.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Diffusion" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Diffusion/DiffusionLFR.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Diffusion" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Diffusion/DiffusionLFRNS.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Diffusion" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Diffusion/DiffusionIP.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Driver.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/DriverAdaptive.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/DriverArnoldi.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/DriverModifiedArnoldi.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/DriverParareal.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/DriverStandard.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/DriverSteadyState.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/EquationSystem.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/Filter.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterAeroForces.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterAverageFields.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterMaxMinFields.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterCheckpoint.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterEnergy1D.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterEnergy.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterError.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterHistoryPoints.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterInterfaces.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterIntegral.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterMean.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterModalEnergy.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterMovingAverage.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterFieldConvert.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterThresholdMax.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterThresholdMin.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Filters" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Filters/FilterBodyFittedVelocity.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/RiemannSolvers" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/RiemannSolvers/RiemannSolver.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/RiemannSolvers" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/RiemannSolvers/UpwindSolver.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/RiemannSolvers" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/RiemannSolvers/UpwindLDGSolver.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/SolverUtils.hpp")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/SolverUtilsDeclspec.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/UnsteadySystem.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Forcing" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Forcing/Forcing.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Forcing" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Forcing/ForcingAbsorption.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Forcing" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Forcing/ForcingBody.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Forcing" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Forcing/ForcingMovingReferenceFrame.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Forcing" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Forcing/ForcingNoise.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/SolverUtils/Forcing" TYPE FILE FILES "/Users/schun/nektar/library/SolverUtils/Forcing/ForcingProgrammatic.h")
endif()

