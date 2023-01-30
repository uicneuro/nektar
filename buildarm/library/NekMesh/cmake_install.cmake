# Install script for directory: /Users/schun/nektar/library/NekMesh

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

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibnekmeshx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/buildarm/library/NekMesh/libNekMesh.5.3.0.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.5.3.0.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xlibnekmeshx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY OPTIONAL FILES "/Users/schun/nektar/buildarm/library/NekMesh/libNekMesh.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LocalRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/SpatialDomains"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/StdRegions"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -delete_rpath "/Users/schun/nektar/buildarm/library/LibUtilities"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    execute_process(COMMAND /usr/bin/install_name_tool
      -add_rpath "/Users/schun/nektar/buildarm/dist/lib/nektar++/"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libNekMesh.dylib")
    endif()
  endif()
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/NekMeshDeclspec.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/Module.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputGmsh.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputNek.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputNek5000.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputNekpp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputPly.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputSem.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputSwan.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/InputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/InputModules/InputStarTec.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/OutputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/OutputModules/OutputGmsh.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/OutputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/OutputModules/OutputNekpp.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/OutputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/OutputModules/OutputSTL.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/OutputModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/OutputModules/OutputStdOut.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessBL.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessCurve.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessCurvedEdges.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessCyl.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessDetectSurf.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessExtractSurf.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessExtractTetPrismInterface.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessJac.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessLinkCheck.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessLinear.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessPerAlign.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessScalar.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessSpherigon.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessTetSplit.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessOptiExtract.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessInsertSurface.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessExtrude.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules/ProcessVarOpti" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessVarOpti/ProcessVarOpti.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules/ProcessVarOpti" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessVarOpti/NodeOpti.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/Module/ProcessModules/ProcessVarOpti" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/Module/ProcessModules/ProcessVarOpti/ElUtil.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Node.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Edge.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Face.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Element.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Composite.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Mesh.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Point.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Line.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Triangle.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Quadrilateral.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Tetrahedron.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Pyramid.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Prism.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/Hexahedron.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/HOAlignment.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/MeshElements" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/MeshElements/ElementConfig.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADObject.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADSystem.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADVert.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADCurve.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/CADSurf.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/libNekMesh/CADSystem" TYPE FILE FILES "/Users/schun/nektar/library/NekMesh/CADSystem/ProcessLoadCAD.h")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xdevx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/nektar++/NekMesh" TYPE DIRECTORY FILES "/Users/schun/nektar/library/NekMesh/./" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\.hpp$")
endif()

