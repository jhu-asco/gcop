link_directories(${RPC_BINARY_DIR}/lib) 
# Use these dependency libs when linking with gcop
# EST_LIBS VIEWER_LIBS SYS_LIBS  UTIL_LIBS  ALGO_LIBS ALL_LIBS
message(" EST_LIBS: ${EST_LIBS}\n")
message(" VIEWER_LIBS: ${VIEWER_LIBS}\n")
message(" SYS_LIBS: ${SYS_LIBS}\n")
message(" UTIL_LIBS: ${UTIL_LIBS}\n")
message(" ALGO_LIBS: ${ALGO_LIBS}\n")
if (BUILD_TESTS)
	add_executable(particle2dtest particle2dtest.cc)
	target_link_libraries(particle2dtest  gcop_systems gcop_views gcop_algos ${ALL_LIBS})# We could just use VIEWER,ALGO,SYS,UTIL

	add_executable(unicycletest unicycletest.cc)
	target_link_libraries(unicycletest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(gunicycletest gunicycletest.cc)
	target_link_libraries(gunicycletest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(body2dtest body2dtest.cc)
	target_link_libraries(body2dtest  gcop_systems gcop_views gcop_algos ${ALL_LIBS})

	add_executable(body3dtest body3dtest.cc)
	target_link_libraries(body3dtest  gcop_systems gcop_views gcop_algos ${ALL_LIBS})

	add_executable(qrotortest qrotortest.cc)
	target_link_libraries(qrotortest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(hrotortest hrotortest.cc)
	target_link_libraries(hrotortest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(helitest helitest.cc)
	target_link_libraries(helitest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(helicartest helicartest.cc)
	target_link_libraries(helicartest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(batest batest.cc)
	target_link_libraries(batest gcop  ${ALL_LIBS})

	add_executable(body2dslamtest body2dslamtest.cc)
	target_link_libraries(body2dslamtest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(body2dtracktest body2dtracktest.cc)
	target_link_libraries(body2dtracktest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(body2dtracktest2 body2dtracktest2.cc)
	target_link_libraries(body2dtracktest2  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})


	add_executable(extforceslamtest extforceslamtest.cc)
	target_link_libraries(extforceslamtest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})


	add_executable(cartest cartest.cc)
	target_link_libraries(cartest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(rccartest rccartest.cc)
	target_link_libraries(rccartest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(gnrccartest gnrccartest.cc)
	target_link_libraries(gnrccartest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})


	add_executable(cecartest cecartest.cc)
	target_link_libraries(cecartest  gcop_systems gcop_views gcop_algos gcop_est ${ALL_LIBS})

	add_executable(cechaintest cechaintest.cc)
	target_link_libraries(cechaintest  gcop_systems gcop_views gcop_algos gcop_est ${ALL_LIBS})

	add_executable(spsacartest spsacartest.cc)
	target_link_libraries(spsacartest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(aspsacartest aspsacartest.cc)
	target_link_libraries(aspsacartest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(chaintest chaintest.cc)
	target_link_libraries(chaintest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(chain1test chain1test.cc)
	target_link_libraries(chain1test  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(airmtest airmtest.cc)
	target_link_libraries(airmtest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(airmopttest airmopttest.cc)
	target_link_libraries(airmopttest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(chainopttest chainopttest.cc)
	target_link_libraries(chainopttest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(fixedchaintest fixedchaintest.cc)
	target_link_libraries(fixedchaintest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(pendtest pendtest.cc)
	target_link_libraries(pendtest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(chaineltest chaineltest.cc)
	target_link_libraries(chaineltest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

#add_executable(snakeopttest snakeopttest.cc)
#target_link_libraries(snakeopttest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})

	add_executable(airbottest airbottest.cc)
	target_link_libraries(airbottest  gcop_systems gcop_views gcop_algos  ${ALL_LIBS})


# Estimation-related tests

	add_executable(ukftest ukftest.cc)
	target_link_libraries(ukftest gcop_est gcop_systems ${ALL_LIBS})

	add_executable(imuukftest imuukftest.cc)
	target_link_libraries(imuukftest gcop_est gcop_systems ${ALL_LIBS})

#add_executable(gptest gptest.cc)
#target_link_libraries(gptest gcop_est gcop_systems gcop_views gcop_algos ${ALL_LIBS})

#add_executable(gpadd gpadd.cc)
#target_link_libraries(gpadd gcop ${ALL_LIBS})


# Cross-entropy simple example

	add_executable(cetest cetest.cc)
	target_link_libraries(cetest gcop_est gcop_systems gcop_est ${ALL_LIBS})

	add_executable(airplanetest airplanetest.cc)
	target_link_libraries(airplanetest gcop_est gcop_systems gcop_views gcop_algos ${ALL_LIBS})

# Camera testing 6dof pose example
#add_executable(camtest camtest.cc)
#target_link_libraries(camtest gcop_est gcop_systems gcop_views gcop_algos ${ALL_LIBS})
#INSTALL(TARGETS gcoptest
#        RUNTIME
#        DESTINATION bin)
endif (BUILD_TESTS)
