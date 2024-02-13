# CMake generated Testfile for 
# Source directory: /Users/schun/nektar/solvers/LinearElasticSolver
# Build directory: /Users/schun/nektar/build/solvers/LinearElasticSolver
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(LinearElasticSolver_L-domain "/Users/schun/nektar/build/tests/Tester" "/Users/schun/nektar/solvers/LinearElasticSolver/Tests/L-domain.tst")
set_tests_properties(LinearElasticSolver_L-domain PROPERTIES  _BACKTRACE_TRIPLES "/Users/schun/nektar/cmake/NektarCommon.cmake;263;ADD_TEST;/Users/schun/nektar/solvers/LinearElasticSolver/CMakeLists.txt;18;ADD_NEKTAR_TEST;/Users/schun/nektar/solvers/LinearElasticSolver/CMakeLists.txt;0;")
add_test(LinearElasticSolver_L-domain-par "/Users/schun/nektar/build/tests/Tester" "/Users/schun/nektar/solvers/LinearElasticSolver/Tests/L-domain-par.tst")
set_tests_properties(LinearElasticSolver_L-domain-par PROPERTIES  _BACKTRACE_TRIPLES "/Users/schun/nektar/cmake/NektarCommon.cmake;263;ADD_TEST;/Users/schun/nektar/solvers/LinearElasticSolver/CMakeLists.txt;23;ADD_NEKTAR_TEST;/Users/schun/nektar/solvers/LinearElasticSolver/CMakeLists.txt;0;")
