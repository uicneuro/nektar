# CMake generated Testfile for 
# Source directory: /Users/schun/nektar/library/UnitTests/LibUtilities
# Build directory: /Users/schun/nektar/buidtest/library/UnitTests/LibUtilities
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(LibUtilitiesUnitTests "/Users/schun/nektar/buidtest/library/UnitTests/LibUtilities/LibUtilitiesUnitTests" "--detect_memory_leaks=0")
set_tests_properties(LibUtilitiesUnitTests PROPERTIES  _BACKTRACE_TRIPLES "/Users/schun/nektar/library/UnitTests/LibUtilities/CMakeLists.txt;23;ADD_TEST;/Users/schun/nektar/library/UnitTests/LibUtilities/CMakeLists.txt;0;")
subdirs("LinearAlgebra")
subdirs("VmathTimer")
