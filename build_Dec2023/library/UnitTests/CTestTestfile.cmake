# CMake generated Testfile for 
# Source directory: /Users/schun/nektar/library/UnitTests
# Build directory: /Users/schun/nektar/build/library/UnitTests
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(UnitTests "/Users/schun/nektar/build/library/UnitTests/UnitTests" "--detect_memory_leaks=0")
set_tests_properties(UnitTests PROPERTIES  _BACKTRACE_TRIPLES "/Users/schun/nektar/library/UnitTests/CMakeLists.txt;20;ADD_TEST;/Users/schun/nektar/library/UnitTests/CMakeLists.txt;0;")
subdirs("LibUtilities")
subdirs("LocalRegions")
subdirs("Collections")
