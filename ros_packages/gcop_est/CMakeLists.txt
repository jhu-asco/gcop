cmake_minimum_required(VERSION 2.4.6)
include($ENV{ROS_ROOT}/core/rosbuild/rosbuild.cmake)

ADD_DEFINITIONS(-g  -pthread -fPIC)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC") 

SET(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# Set the build type.  Options are:
#  Coverage       : w/ debug symbols, w/o optimization, w/ code-coverage
#  Debug          : w/ debug symbols, w/o optimization
#  Release        : w/o debug symbols, w/ optimization
#  RelWithDebInfo : w/ debug symbols, w/ optimization
#  MinSizeRel     : w/o debug symbols, w/ optimization, stripped binaries
#set(ROS_BUILD_TYPE RelWithDebInfo)

rosbuild_init()

#set the default path for built executables to the "bin" directory
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
#set the default path for built libraries to the "lib" directory
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

#uncomment if you have defined messages
#rosbuild_genmsg()
#uncomment if you have defined services
#rosbuild_gensrv()



find_package(GCOP REQUIRED)
include_directories(${GCOP_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GCOP_LIBRARIES})

find_package(OpenCV REQUIRED core highgui imgproc calib3d)
include_directories( ${OpenCV_INCLUDE_DIRS} )
set(LIBS ${LIBS} ${OpenCV_LIBS})
message(STATUS "Hello ${OpenCV_LIBS}")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(GLU REQUIRED)
include_directories(${GLU_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GLU_LIBRARIES} ${GLU_LIBRARY})

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GLUT_LIBRARY})

# Nodelet library
rosbuild_add_library(${PROJECT_NAME} src/nodelets/estimator1.cpp)
target_link_libraries(${PROJECT_NAME} ${LIBS})
#src/nodelets/debayer.cpp
#			src/nodelets/rectify.cpp
#			src/nodelets/crop_decimate.cpp
#			src/libimage_proc/advertisement_checker.cpp
#			src/nodelets/edge_aware.cpp
#			src/nodelets/yuv422.cpp

#common commands for building c++ executables and libraries
#rosbuild_add_library(${PROJECT_NAME} src/example.cpp)
#target_link_libraries(${PROJECT_NAME} another_library)
#rosbuild_add_boost_directories()
#rosbuild_link_boost(${PROJECT_NAME} thread)
#rosbuild_add_executable(example examples/example.cpp)
#target_link_libraries(example ${PROJECT_NAME})
