# This is basically the overall name of the project in Visual Studio this is the name of the Solution File


# For every executable you have with a main method you should have an add_executable line below.
# For every add executable line you should list every .cpp and .h file you have associated with that executable.  



# You shouldn't have to modify anything below this line 
########################################################

message("BULLET INCLUDE: ${BULLET_INCLUDE_DIR}")
message("Bullet LIBS: ${BULLET_LIBRARIES}")

INCLUDE_DIRECTORIES(
${BULLET_PHYSICS_SOURCE_DIR}/Demos/OpenGL 
${GLUT_INCLUDE_DIR}
)
#BulletDynamics BulletCollision LinearMath 
LINK_LIBRARIES(
 OpenGLSupport ${GLUT_glut_LIBRARY} ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY}
 gcop_algos gcop_bulletsystems gcop_systems ${BULLET_LIBRARIES} ${ALGO_LIBS} ${SYS_LIBS} ${UTIL_LIBS}
)

	#heightfield128x128.cpp
ADD_EXECUTABLE(AppVehicleDemo
	VehicleDemo.cpp
	main.cpp
)

IF (INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
			SET_TARGET_PROPERTIES(AppVehicleDemo PROPERTIES  DEBUG_POSTFIX "_Debug")
			SET_TARGET_PROPERTIES(AppVehicleDemo PROPERTIES  MINSIZEREL_POSTFIX "_MinsizeRel")
			SET_TARGET_PROPERTIES(AppVehicleDemo PROPERTIES  RELWITHDEBINFO_POSTFIX "_RelWithDebugInfo")
ENDIF(INTERNAL_ADD_POSTFIX_EXECUTABLE_NAMES)
