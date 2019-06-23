find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)

foreach(PROG
	test_atomcgconverter
  test_bead
  test_beadtriple
  test_basebead
  test_beadmap
  test_beadmotif_algorithms
  test_beadmotif_base
  test_beadmotif_connector
  test_beadstructure_base
  test_beadstructure_algorithms
  test_bondedstatistics
	test_cgmoleculestencil
	test_csg_topology
  test_interaction
  test_lammpsdatareader 
  test_lammpsdumpreaderwriter
	test_nblist_3body
  test_nblistgrid_3body
  test_pdbreader
  test_tabulatedpotential
  test_topologytypecontainer
  test_triplelist
  test_xyzreaderwriter)

  file(GLOB ${PROG}_SOURCES ${PROG}.cc)
  add_executable(unit_${PROG} ${${PROG}_SOURCES})
  target_link_libraries(unit_${PROG} votca_csg_boltzmann votca_csg ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  target_compile_definitions(unit_${PROG} PRIVATE BOOST_TEST_DYN_LINK)
  add_test(unit_${PROG} unit_${PROG})
  # run tests for csg (for coverage) as well
  set_tests_properties(unit_${PROG} PROPERTIES LABELS "csg;votca")

  if(VALGRIND_EXECUTABLE)
    # Do not create memory test for this specific unit test because it takes 
    # too long
    if(NOT ${PROG} MATCHES "test_tabulatedpotential")
    add_test(NAME memory_${PROG} COMMAND "${VALGRIND_EXECUTABLE}" --tool=memcheck --leak-check=yes --show-reachable=yes --num-callers=20 --track-fds=yes --error-exitcode=1 $<TARGET_FILE:unit_${PROG}>)
    set_tests_properties(memory_${PROG} PROPERTIES LABELS "csg;votca")
  endif(NOT ${PROG} MATCHES "test_tabulatedpotential")
  endif(VALGRIND_EXECUTABLE)

endforeach(PROG)
