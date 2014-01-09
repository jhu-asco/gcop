link_directories(${RPC_BINARY_DIR}/lib) 

add_executable(particle2dtest particle2dtest.cc)
target_link_libraries(particle2dtest gcop  ${LIBS})

add_executable(unicycletest unicycletest.cc)
target_link_libraries(unicycletest gcop  ${LIBS})

add_executable(gunicycletest gunicycletest.cc)
target_link_libraries(gunicycletest gcop  ${LIBS})

add_executable(body3dtest body3dtest.cc)
target_link_libraries(body3dtest gcop  ${LIBS})

add_executable(qrotortest qrotortest.cc)
target_link_libraries(qrotortest gcop  ${LIBS})

add_executable(hrotortest hrotortest.cc)
target_link_libraries(hrotortest gcop  ${LIBS})


add_executable(helitest helitest.cc)
target_link_libraries(helitest gcop  ${LIBS})

add_executable(helicartest helicartest.cc)
target_link_libraries(helicartest gcop  ${LIBS})

add_executable(batest batest.cc)
target_link_libraries(batest gcop  ${LIBS})

add_executable(cartest cartest.cc)
target_link_libraries(cartest gcop  ${LIBS})

add_executable(rccartest rccartest.cc)
target_link_libraries(rccartest gcop  ${LIBS})

add_executable(cecartest cecartest.cc)
target_link_libraries(cecartest gcop  ${LIBS})

add_executable(cechaintest cechaintest.cc)
target_link_libraries(cechaintest gcop  ${LIBS})


add_executable(chaintest chaintest.cc)
target_link_libraries(chaintest gcop  ${LIBS})

add_executable(chain1test chain1test.cc)
target_link_libraries(chain1test gcop  ${LIBS})

add_executable(airmtest airmtest.cc)
target_link_libraries(airmtest gcop  ${LIBS})

add_executable(airmopttest airmopttest.cc)
target_link_libraries(airmopttest gcop  ${LIBS})

add_executable(chainopttest chainopttest.cc)
target_link_libraries(chainopttest gcop  ${LIBS})

add_executable(fixedchaintest fixedchaintest.cc)
target_link_libraries(fixedchaintest gcop  ${LIBS})

add_executable(pendtest pendtest.cc)
target_link_libraries(pendtest gcop  ${LIBS})


add_executable(chaineltest chaineltest.cc)
target_link_libraries(chaineltest gcop  ${LIBS})

#add_executable(snakeopttest snakeopttest.cc)
#target_link_libraries(snakeopttest gcop  ${LIBS})

add_executable(airbottest airbottest.cc)
target_link_libraries(airbottest gcop  ${LIBS})


# Estimation-related tests

add_executable(ukftest ukftest.cc)
target_link_libraries(ukftest gcop ${LIBS})

add_executable(imuukftest imuukftest.cc)
target_link_libraries(imuukftest gcop ${LIBS})

#add_executable(gptest gptest.cc)
#target_link_libraries(gptest gcop ${LIBS})

#add_executable(gpadd gpadd.cc)
#target_link_libraries(gpadd gcop ${LIBS})


# Cross-entropy simple example

add_executable(cetest cetest.cc)
target_link_libraries(cetest gcop ${LIBS})

add_executable(airplanetest airplanetest.cc)
target_link_libraries(airplanetest gcop ${LIBS})


#INSTALL(TARGETS gcoptest
#        RUNTIME
#        DESTINATION bin)

