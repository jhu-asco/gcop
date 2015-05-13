link_directories(${GCOP_BINARY_DIR}/lib) 
# Use these dependency libs when linking with gcop
# EST_LIBS VIEWER_LIBS SYS_LIBS  UTIL_LIBS  ALGO_LIBS ALL_LIBS
message(" EST_LIBS: ${EST_LIBS}\n")
message(" VIEWER_LIBS: ${VIEWER_LIBS}\n")
message(" SYS_LIBS: ${SYS_LIBS}\n")
message(" UTIL_LIBS: ${UTIL_LIBS}\n")
message(" ALGO_LIBS: ${ALGO_LIBS}\n")
if (BUILD_TESTS)
	add_executable(particle2dtest particle2dtest.cc)
	target_link_libraries(particle2dtest  gcop_algos gcop_views gcop_systems ${ALL_LIBS})# We could just use VIEWER,ALGO,SYS,UTIL

	add_executable(unicycletest unicycletest.cc)
	target_link_libraries(unicycletest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(gunicycletest gunicycletest.cc)
	target_link_libraries(gunicycletest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(body2dtest body2dtest.cc)
	target_link_libraries(body2dtest  gcop_algos gcop_views gcop_systems ${ALL_LIBS})

	add_executable(body3dtest body3dtest.cc)
	target_link_libraries(body3dtest  gcop_algos gcop_views gcop_systems ${ALL_LIBS})

	add_executable(body3dstab body3dstab.cc)
	target_link_libraries(body3dstab  gcop_algos gcop_views gcop_systems ${ALL_LIBS})


# use qrotor.cfg as input
	add_executable(qrotortest qrotortest.cc)
	target_link_libraries(qrotortest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(hrotortest hrotortest.cc)
	target_link_libraries(hrotortest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})


IF (PQP_FOUND)
#	add_executable(hrotordem hrotordem.cc)
#	target_link_libraries(hrotordem  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

#	add_executable(body3ddem body3ddem.cc)
#	target_link_libraries(body3ddem  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(body3drhc body3drhc.cc)
	target_link_libraries(body3drhc  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(body3ddemstab body3ddemstab.cc)
	target_link_libraries(body3ddemstab  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(body3dcestab body3dcestab.cc)
	target_link_libraries(body3dcestab  gcop_algos gcop_views gcop_systems ${ALL_LIBS})

IF (DSL_FOUND)

	add_executable(body3dcedemstab body3dcedemstab.cc)
	target_link_libraries(body3dcedemstab gcop_algos gcop_views gcop_systems ${ALL_LIBS})

ENDIF (DSL_FOUND)


#	add_executable(body3dstops body3dstops.cc)
#	target_link_libraries(body3dstops  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})


	add_executable(body3dunc body3dunc.cc)
	target_link_libraries(body3dunc  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

#	add_executable(body3davoids body3davoids.cc)
#	target_link_libraries(body3davoids  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

ENDIF (PQP_FOUND)


	add_executable(hrotorshell hrotorshell.cc)
	target_link_libraries(hrotorshell  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(body3dshell body3dshell.cc)
	target_link_libraries(body3dshell  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(helitest helitest.cc)
	target_link_libraries(helitest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(helicartest helicartest.cc)
	target_link_libraries(helicartest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

#	add_executable(batest batest.cc)
#	target_link_libraries(batest gcop_algos gcop_views gcop_systems ${ALL_LIBS})

	#add_executable(body2dslamtest body2dslamtest.cc)
	#target_link_libraries(body2dslamtest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

#	add_executable(body2dtracktest body2dtracktest.cc)
#	target_link_libraries(body2dtracktest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(body2dtracktest2 body2dtracktest2.cc)
	target_link_libraries(body2dtracktest2 gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(body3dtracktest body3dtracktest.cc)
	target_link_libraries(body3dtracktest gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(kinbody3dtracktest kinbody3dtracktest.cc)
	target_link_libraries(kinbody3dtracktest gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(kinrccartracktest kinrccartracktest.cc)
	target_link_libraries(kinrccartracktest gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(kinbodyprojtracktest kinbodyprojtracktest.cc)
	target_link_libraries(kinbodyprojtracktest gcop_algos gcop_views gcop_systems  ${ALL_LIBS})
	
        #add_executable(extforceslamtest extforceslamtest.cc)
	#target_link_libraries(extforceslamtest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(cartest cartest.cc)
	target_link_libraries(cartest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(rccartest rccartest.cc)
	target_link_libraries(rccartest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(rccar_flatparamverify rccar_flatparamverify.cc)
	target_link_libraries(rccar_flatparamverify gcop_algos gcop_views gcop_systems ${ALL_LIBS})

	add_executable(gnrccartest gnrccartest.cc)
	target_link_libraries(gnrccartest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})
if (USE_BULLET)
	add_executable(ddpbulletrccartest ddpbulletrccartest.cc)
	target_link_libraries(ddpbulletrccartest  gcop_algos gcop_views gcop_bulletsystems gcop_systems ${BULLET_LIBRARIES} ${ALL_LIBS})

	add_executable(sddpbulletrccartest sddpbulletrccartest.cc)
	target_link_libraries(sddpbulletrccartest  gcop_algos gcop_views gcop_bulletsystems gcop_systems ${BULLET_LIBRARIES} ${ALL_LIBS})

	add_executable(bulletrccartest bulletrccartest.cc)
	target_link_libraries(bulletrccartest  gcop_algos gcop_views gcop_bulletsystems gcop_systems ${BULLET_LIBRARIES} ${ALL_LIBS})

	add_executable(gnbulletrccartest gnbulletrccartest.cc)
	target_link_libraries(gnbulletrccartest  gcop_algos gcop_views gcop_bulletsystems gcop_systems ${BULLET_LIBRARIES} ${ALL_LIBS})

	add_executable(cebulletcartest cebulletcartest.cc)
	target_link_libraries(cebulletcartest  gcop_algos gcop_views gcop_bulletsystems gcop_systems ${BULLET_LIBRARIES} ${ALL_LIBS})

	add_executable(spsabulletcartest spsabulletcartest.cc)
	target_link_libraries(spsabulletcartest  gcop_algos gcop_views gcop_bulletsystems gcop_systems ${BULLET_LIBRARIES} ${ALL_LIBS})
endif (USE_BULLET)

	add_executable(gnrccarestimation gnrccarestimation.cc)
	target_link_libraries(gnrccarestimation  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(cecartest cecartest.cc)
	target_link_libraries(cecartest  gcop_algos gcop_views gcop_systems gcop_est ${ALL_LIBS})

	#add_executable(cechaintest cechaintest.cc)
	#target_link_libraries(cechaintest  gcop_algos gcop_views gcop_systems gcop_est ${ALL_LIBS})

	add_executable(spsacartest spsacartest.cc)
	target_link_libraries(spsacartest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(aspsacartest aspsacartest.cc)
	target_link_libraries(aspsacartest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(chaintest chaintest.cc)
	target_link_libraries(chaintest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	#add_executable(chain1test chain1test.cc)
	#target_link_libraries(chain1test  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})


	#add_executable(airmtest airmtest.cc)
	#target_link_libraries(airmtest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(airmopttest airmopttest.cc)
	target_link_libraries(airmopttest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	#add_executable(chainopttest chainopttest.cc)
	#target_link_libraries(chainopttest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})




# to be used with bin/fixedchain1.cfg
	add_executable(fixedchaintest fixedchaintest.cc)
	target_link_libraries(fixedchaintest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

#	add_executable(pendtest pendtest.cc)
#	target_link_libraries(pendtest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

#	add_executable(chaineltest chaineltest.cc)
#	target_link_libraries(chaineltest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

#add_executable(snakeopttest snakeopttest.cc)
#target_link_libraries(snakeopttest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

	add_executable(airbottest airbottest.cc)
	target_link_libraries(airbottest  gcop_algos gcop_views gcop_systems  ${ALL_LIBS})


# Estimation-related tests
	#add_executable(ukftest ukftest.cc)
	#target_link_libraries(ukftest gcop_est gcop_systems ${ALL_LIBS})

	#add_executable(imuukftest imuukftest.cc)
	#target_link_libraries(imuukftest gcop_est gcop_systems ${ALL_LIBS})

add_executable(gptest gptest.cc)
target_link_libraries(gptest gcop_est gcop_algos gcop_views gcop_systems ${ALL_LIBS})

add_executable(gpadd gpadd.cc)
target_link_libraries(gpadd gcop_est gcop_algos gcop_views gcop_systems ${ALL_LIBS})

add_executable(imuekftest imuekftest.cc)
target_link_libraries(imuekftest gcop_est gcop_algos gcop_views gcop_systems ${ALL_LIBS})

add_executable(insekftest insekftest.cc)
target_link_libraries(insekftest gcop_est gcop_algos gcop_views gcop_systems ${ALL_LIBS})

add_executable(point3dkftest point3dkftest.cc)
target_link_libraries(point3dkftest gcop_est gcop_algos gcop_views gcop_systems ${ALL_LIBS})

add_executable(point3dstab point3dstab.cc)
target_link_libraries(point3dstab gcop_est gcop_algos gcop_views gcop_systems ${ALL_LIBS})


add_executable(helirhc helirhc.cc)
target_link_libraries(helirhc gcop_est gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

# with bin/jhurov.cfg
add_executable(uuvtest uuvtest.cc)
target_link_libraries(uuvtest gcop_est gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

#add_executable(armtest armtest.cc)
#target_link_libraries(armtest gcop_est gcop_algos gcop_views gcop_systems  ${ALL_LIBS})

# Cross-entropy simple example

add_executable(cetest cetest.cc)
target_link_libraries(cetest gcop_est gcop_systems gcop_est ${ALL_LIBS})

#add_executable(airplanetest airplanetest.cc)
#target_link_libraries(airplanetest gcop_est gcop_algos gcop_views gcop_systems ${ALL_LIBS})

# Camera testing 6dof pose example
if(OPENCV_FOUND)
  add_executable(camtest camtest.cc)
  target_link_libraries(camtest gcop_est gcop_algos gcop_views gcop_systems ${cam_LIBS} ${ALL_LIBS})
endif(OPENCV_FOUND)

#INSTALL(TARGETS gcoptest
#        RUNTIME
#        DESTINATION bin)
endif (BUILD_TESTS)

if (USE_BULLET)
  if (BULLET_FOUND)
    add_subdirectory (SimpleOpenGL3)
    add_subdirectory (BulletVehicleEstimation)
    add_subdirectory (BulletVehicleEstimation1)
    add_subdirectory (BulletSTLDemo)
    add_subdirectory (BulletZVehicleTest)
  endif (BULLET_FOUND)
endif (USE_BULLET)
