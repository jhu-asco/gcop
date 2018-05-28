link_directories(${GCOP_BINARY_DIR}/lib) 
# Use these dependency libs when linking with gcop
# EST_LIBS VIEWER_LIBS SYS_LIBS  UTIL_LIBS  ALGO_LIBS ALL_LIBS
message(${SYS_LIBS})
if (BUILD_TESTS)
  add_executable(test_load_eigen_matrix test_load_eigen_matrix.cpp)
  target_compile_definitions(test_load_eigen_matrix PUBLIC DATA_PATH="${CMAKE_SOURCE_DIR}/data")
  target_link_libraries(test_load_eigen_matrix ${SYS_LIBS} ${GTEST_BOTH_LIBRARIES})
  add_test(test_load_eigen_matrix test_load_eigen_matrix)
  if (casadi_FOUND)
    add_executable(test_casadi_system test_casadi_system.cc)
    target_link_libraries(test_casadi_system gcop_systems ${SYS_LIBS} ${UTIL_LIBS} ${GTEST_BOTH_LIBRARIES})
    add_test(test_casadi_system test_casadi_system)

    add_executable(test_quad_casadi_system test_quad_casadi_system.cpp)
    target_link_libraries(test_quad_casadi_system gcop_systems ${SYS_LIBS} ${UTIL_LIBS} ${GTEST_BOTH_LIBRARIES})
    add_test(test_quad_casadi_system test_quad_casadi_system)

    add_executable(test_airm_ff_system test_airm_ff_system.cpp)
    target_compile_definitions(test_airm_ff_system PUBLIC DATA_PATH="${CMAKE_SOURCE_DIR}/data")
    target_link_libraries(test_airm_ff_system gcop_systems ${SYS_LIBS} ${UTIL_LIBS} ${GTEST_BOTH_LIBRARIES})
    add_test(test_airm_ff_system test_airm_ff_system)

    add_executable(test_fully_connected_layer test_fully_connected_layer.cpp)
    target_compile_definitions(test_fully_connected_layer PUBLIC DATA_PATH="${CMAKE_SOURCE_DIR}/data")
    target_link_libraries(test_fully_connected_layer gcop_systems ${SYS_LIBS} ${UTIL_LIBS} ${GTEST_BOTH_LIBRARIES})
    add_test(test_fully_connected_layer test_fully_connected_layer)

    add_executable(test_airm_residual_network_model test_airm_residual_network_model.cpp)
    target_compile_definitions(test_airm_residual_network_model PUBLIC DATA_PATH="${CMAKE_SOURCE_DIR}/data")
    target_link_libraries(test_airm_residual_network_model gcop_systems ${SYS_LIBS} ${UTIL_LIBS} ${GTEST_BOTH_LIBRARIES})
    add_test(test_airm_residual_network_model test_airm_residual_network_model)
  endif (casadi_FOUND)
endif (BUILD_TESTS)


