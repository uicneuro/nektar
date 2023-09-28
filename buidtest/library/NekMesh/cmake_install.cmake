# Install script for directory: /Users/schun/nektar/library/NekMesh

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

if(CMAKE_INSTALL_COMPONENT STREQUAL "libnekmesh" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/buidtest/library/NekMesh/libNekMesh.5.3.0.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buidtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buidtest/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/buidtest/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/buidtest/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/buidtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buidtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buidtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "libnekmesh" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/buidtest/library/NekMesh/libNekMesh.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buidtest/dist/lib/nektar++"
      -delete_rpath "/Users/schun/nektar/buidtest/library/LocalRegions"
      -delete_rpath "/Users/schun/nektar/buidtest/library/SpatialDomains"
      -delete_rpath "/Users/schun/nektar/buidtest/library/StdRegions"
      -delete_rpath "/Users/schun/nektar/buidtest/library/LibUtilities"
      -add_rpath "/Users/schun/nektar/buidtest/dist/lib"
      -add_rpath "/Users/schun/nektar/buidtest/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/NekMeshDeclspec.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/Module.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputGmsh.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputNek.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputNek5000.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputNekpp.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputPly.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputSem.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputSwan.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputStarTec.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/OutputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/OutputModules/OutputGmsh.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/OutputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/OutputModules/OutputNekpp.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/OutputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/OutputModules/OutputSTL.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/OutputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/OutputModules/OutputStdOut.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessBL.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessCurve.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessCurvedEdges.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessCyl.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessDetectSurf.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessExtractSurf.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessExtractTetPrismInterface.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessJac.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessLinkCheck.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessLinear.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessPerAlign.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessScalar.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessSpherigon.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessTetSplit.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessOptiExtract.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessInsertSurface.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessExtrude.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules/ProcessVarOpti" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessVarOpti/ProcessVarOpti.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules/ProcessVarOpti" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessVarOpti/NodeOpti.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules/ProcessVarOpti" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessVarOpti/ElUtil.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Node.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Edge.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Face.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Element.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Composite.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Mesh.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Point.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Line.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Triangle.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Quadrilateral.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Tetrahedron.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Pyramid.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Prism.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Hexahedron.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/HOAlignment.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/ElementConfig.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADObject.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADSystem.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADVert.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADCurve.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADSurf.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/ProcessLoadCAD.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "dev" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/NekMesh" TYPE DIRECTORY FILES "/Users/schun/nektar/library/NekMesh/./" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\.hpp$")
endif()

