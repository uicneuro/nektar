# - Nektar++ Config File
#
# Use this module by invoking find_package with the form:
#  FIND_PACKAGE(Nektar++
#    [version] [EXACT]      # Minimum or EXACT version e.g. 1.36.0
#    [REQUIRED]             # Fail with error if Nektar++ is not found
#    )                      #
#
# This sets the following variables:
#  Nektar++_FOUND             - True if headers and requested libraries were found
#  Nektar++_VERSION           - Nektar++_VERSION
#  Nektar++_INCLUDE_DIRS      - Nektar++ include directories
#  Nektar++_LIBRARY_DIRS      - Link directories for Nektar++ libraries
#  Nektar++_DEFINITIONS       - Nektar++ build flags
#  Nektar++_LIBRARIES         - Nektar++ component libraries to be linked
#
#  Nektar++_TP_INCLUDE_DIRS   - Nektar++ ThirdParty include directories
#  Nektar++_TP_LIBRARY_DIRS   - Link directories for Nektar++ ThirdParty libraries
#  Nektar++_TP_LIBRARIES      - Nektar++ ThirdParty libraries to be linked
#
# Example Use:
#  FIND_PACKAGE(Nektar++ REQUIRED)
#  ADD_DEFINITIONS(${NEKTAR++_DEFINITIONS})
#  INCLUDE_DIRECTORIES(${NEKTAR++_INCLUDE_DIRS} ${NEKTAR++_TP_INCLUDE_DIRS})
#  LINK_DIRECTORIES(${NEKTAR++_LIBRARY_DIRS} ${NEKTAR++_TP_LIBRARY_DIRS})
#  TARGET_LINK_LIBRARIES(${ProjectName} ${NEKTAR++_LIBRARIES} ${NEKTAR++_TP_LIBRARIES})
#

# set basic variables
SET(NEKTAR++_FOUND "ON")
SET(NEKTAR++_VERSION "5.3.0")
SET(NEKTAR++_ROOT_DIR "/Users/schun/nektar/build/dist")
SET(NEKTAR++_INCLUDE_DIRS "/Users/schun/nektar/build/dist/include/nektar++")
SET(NEKTAR++_LIBRARY_DIRS "${NEKTAR++_ROOT_DIR}/lib")
SET(NEKTAR++_DEFINITIONS -DNEKTAR_USE_SCOTCH -DTIXML_USE_STL -DNEKTAR_MEMORY_POOL_ENABLED)
SET(NEKTAR++_LIBRARIES LibUtilities;LocalRegions;SpatialDomains;StdRegions;Collections;MultiRegions;MatrixFreeOps;SolverUtils;GlobalMapping;FieldUtils;libNekMesh)

SET(NEKTAR++_TP_INCLUDE_DIRS "include/nektar++")
SET(NEKTAR++_TP_LIBRARIES "")
SET(NEKTAR++_TP_LIBRARY_DIRS "")

# Set up Boost and other ThirdParty include directories.
SET(Boost_INCLUDE_DIRS "/opt/local/include")
SET(NEKTAR++_TP_INCLUDE_DIRS
    ${NEKTAR++_TP_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})

SET(Boost_LIBRARIES /opt/local/lib/libboost_thread-mt.dylib;/opt/local/lib/libboost_iostreams-mt.dylib;/opt/local/lib/libboost_filesystem-mt.dylib;/opt/local/lib/libboost_system-mt.dylib;/opt/local/lib/libboost_program_options-mt.dylib;/opt/local/lib/libboost_regex-mt.dylib;/opt/local/lib/libboost_chrono-mt.dylib;/opt/local/lib/libboost_date_time-mt.dylib;/opt/local/lib/libboost_atomic-mt.dylib)
SET(Boost_LIBRARY_DIRS "/opt/local/lib")
SET(NEKTAR++_TP_LIBRARIES
    ${NEKTAR++_TP_LIBRARIES} ${Boost_LIBRARIES})

SET(TINYXML_INCLUDE_DIRS "")
SET(NEKTAR++_TP_INCLUDE_DIRS
    ${NEKTAR++_TP_INCLUDE_DIRS} ${TINYXML_INCLUDE_DIRS})

SET(ZLIB_INCLUDE_DIRS "/opt/local/include")
SET(NEKTAR++_TP_INCLUDE_DIRS
    ${NEKTAR++_TP_INCLUDE_DIRS} ${ZLIB_INCLUDE_DIRS})

SET(NEKTAR_USE_MPI "OFF")
IF( NEKTAR_USE_MPI )
    SET(MPI_C_LIBRARIES "")
    SET(MPI_C_INCLUDE_PATH "")
    SET(NEKTAR++_TP_INCLUDE_DIRS
        ${NEKTAR++_TP_INCLUDE_DIRS} ${MPI_C_INCLUDE_PATH})
    SET(NEKTAR++_TP_LIBRARIES
        ${NEKTAR++_TP_LIBRARIES} ${MPI_C_LIBRARIES})
    SET(NEKTAR++_DEFINITIONS
        "${NEKTAR++_DEFINITIONS} -DNEKTAR_USE_MPI")
ENDIF( NEKTAR_USE_MPI )

SET(NEKTAR_USE_SCOTCH "ON")
IF( NEKTAR_USE_SCOTCH )
    SET(SCOTCH_INCLUDE_DIRS "/opt/local/include")
    SET(NEKTAR++_TP_INCLUDE_DIRS
        ${NEKTAR++_TP_INCLUDE_DIRS} ${SCOTCH_INCLUDE_DIRS})
ENDIF( NEKTAR_USE_SCOTCH )

SET(NEKTAR_USE_METIS "OFF")
IF( NEKTAR_USE_METIS )
    SET(METIS_INCLUDE_DIRS "")
    SET(NEKTAR++_TP_INCLUDE_DIRS
        ${NEKTAR++_TP_INCLUDE_DIRS} ${METIS_INCLUDE_DIRS})
ENDIF( NEKTAR_USE_METIS )

SET(NEKTAR_USE_FFTW "OFF")
IF( NEKTAR_USE_FFTW )
    SET(FFTW_INCLUDE_DIR "")
    SET(NEKTAR++_TP_INCLUDE_DIRS
        ${NEKTAR++_TP_INCLUDE_DIRS} ${FFTW_INCLUDE_DIR})
ENDIF( NEKTAR_USE_FFTW )

SET(NEKTAR_USE_ARPACK "OFF")
IF( NEKTAR_USE_ARPACK )
    SET(ARPACK_INCLUDE_DIR "")
    SET(NEKTAR++_TP_INCLUDE_DIRS
        ${NEKTAR++_TP_INCLUDE_DIRS} ${ARPACK_INCLUDE_DIR})
ENDIF( NEKTAR_USE_ARPACK )

SET(NEKTAR_USE_VTK "OFF")
IF( NEKTAR_USE_VTK )
    IF (OFF)
        SET(NEKTAR++_TP_INCLUDE_DIRS
            ${NEKTAR++_TP_INCLUDE_DIRS} ${VTK_INCLUDE_DIRS})
    ELSE()
        # Make a FIND_DEPENDENCY call to find the relevant version of
        # VTK and the required modules.
        INCLUDE(CMakeFindDependencyMacro)
        SET(VTK_VERSION "")
        IF (VTK_VERSION EQUAL 9)
            FIND_DEPENDENCY(VTK 9 QUIET COMPONENTS
                FiltersCore IOLegacy IOXML IOImage RenderingCore)
        ELSEIF (VTK_VERSION LESS 9)
            FIND_DEPENDENCY(VTK ${VTK_VERSION} QUIET COMPONENTS
                vtkFiltersGeometry vtkIOLegacy vtkIOXML vtkIOImage vtkRenderingCore)
        ENDIF()
    ENDIF()
ENDIF( NEKTAR_USE_VTK )

SET(NEKTAR_USE_HDF5 "OFF")
IF( NEKTAR_USE_HDF5 )
    SET(HDF5_INCLUDE_DIRS "")
    SET(NEKTAR++_TP_INCLUDE_DIRS
        ${NEKTAR++_TP_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})
ENDIF( NEKTAR_USE_HDF5 )

# find and add Nektar++ libraries
INCLUDE(/Users/schun/nektar/build/dist/lib/nektar++/cmake/Nektar++Libraries.cmake)

# platform dependent options
IF(${CMAKE_SYSTEM} MATCHES "Linux.*")
    set(NEKTAR++_TP_LIBRARIES ${NEKTAR++_TP_LIBRARIES} rt)
    SET(NEKTAR++_DEFINITIONS "${NEKTAR++_DEFINITIONS} -pthread")
ENDIF(${CMAKE_SYSTEM} MATCHES "Linux.*")
