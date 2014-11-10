cmake_minimum_required (VERSION 2.6) 
project (gcop)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

#ADD_DEFINITIONS(-g -O -std=c++0x -pthread)
ADD_DEFINITIONS(-g -O -std=c++0x -pthread -Wl,--no-undefined  -I/usr/include/eigen3)
#ADD_DEFINITIONS(-g -O -pthread)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Robotic Planning and Control")
SET(CPACK_PACKAGE_VENDOR "Marin Kobilarov")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "2")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

## Options for building this Project:
option(BUILD_TESTS "BUILD_TESTS" ON)
# Will add more like disabling est if not needed etc
#option(USE_VIEWER "USE_VIEWER" ON)
#option(DISABLE_EST "DISABLE_ESTIMATION" OFF)

#SET(CPACK_STRIP_FILES "bin/MyExecutable")
#SET(CPACK_SOURCE_STRIP_FILES "")

#find_package(OpenGL REQUIRED)
#find_package(GLUT REQUIRED)
#include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
#target_link_libraries( ${OpenGL_LIBRARIES} ${GLUT_LIBRARY} )

# We break this project gcop into smaller libraries as:
# gcop_systems.so gcop_views.so gcop_algos.so gcop_est.so
# Utils folder is linked into all the above libraries

# Instead of a common variable called LIBS we make separate LIBS for each the submodule library
# 	SYS_LIBS	ALGO_LIBS	VIEWER_LIBS	UTIL_LIBS EST_LIBS
# Any Libs needed specifically by testing binaries is named based on the testing bin name
# Example: Camtest bin needs ABC library then there will be a variable Cam_LIBS with the Libraries needed by the Camtest (It may be empty if nothing is needed)

find_package(PkgConfig REQUIRED)


find_package( PkgConfig )
pkg_check_modules( EIGEN3 REQUIRED eigen3 )
include_directories( ${EIGEN3_INCLUDE_DIRS} )

#find_package(Eigen REQUIRED)
#include_directories(${Eigen_INCLUDE_DIRS})

find_package( OpenCV 2.4 COMPONENTS core highgui imgproc calib3d)
include_directories( ${OpenCV_INCLUDE_DIRS} )
set(EST_LIBS ${EST_LIBS} ${OpenCV_LIBS})
message(STATUS "Hello ${OpenCV_LIBS}")

find_package(DC1394)
include_directories(${DC1394_INCLUDE_DIR})
set(cam_LIBS ${cam_LIBS} ${DC1394_LIBRARIES})
message(STATUS "Hello2 ${LIBS}")

#find_package(LibUSB REQUIRED)
#include_directories(${LIBUSB_INCLUDE_DIR})
#set(LIBS ${LIBS} ${LIBUSB_LIBRARIES})
set(UTIL_LIBS ${UTIL_LIBS} gcop_utils)#Adding gcop_utils to the util libraries

find_package(GLU REQUIRED)
include_directories(${GLU_INCLUDE_DIRS})
set(VIEWER_LIBS ${VIEWER_LIBS} ${GLU_LIBRARIES} ${GLU_LIBRARY})

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
set(VIEWER_LIBS ${VIEWER_LIBS} ${GLUT_LIBRARY})

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
set(VIEWER_LIBS ${VIEWER_LIBS} ${OpenGL_LIBRARIES})

#find_package(Adolc)
#include_directories(${Adolc_INCLUDE_DIRS})
#set(SYS_LIBS ${SYS_LIBS} ${Adolc_LIBRARIES})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
set(UTIL_LIBS ${UTIL_LIBS} ${Boost_LIBRARY})

find_package (Threads)
set (VIEWER_LIBS ${VIEWER_LIBS} ${CMAKE_THREAD_LIBS_INIT})

find_package(TinyXML)
include_directories(${TINYXML_INCLUDE_DIR})
set(UTIL_LIBS ${UTIL_LIBS} ${TINYXML_LIBRARIES})


#FIND_PACKAGE(Ceres)
#INCLUDE_DIRECTORIES(${CERES_INCLUDE_DIRS})
#set(ALGO_LIBS ${ALGO_LIBS})


set(ALL_LIBS ${ALGO_LIBS} ${SYS_LIBS} ${VIEWER_LIBS} ${EST_LIBS} ${UTIL_LIBS})

include_directories(
  lib
  lib/algos
  lib/est
  lib/systems
  lib/utils
  lib/views)

add_subdirectory (lib) 
add_subdirectory (bin)
#add_subdirectory (docs)

#SET(CPACK_PACKAGE_EXECUTABLES test "An Example Trajectory Planner")

#set(CPACK_SOURCE_INSTALLED_DIRECTORIES "${CMAKE_SOURCE_DIR}/lib;./lib")
set(CPACK_SOURCE_IGNORE_FILES
"/docs/*;/logs/*;/ros_packages/;/build/;/.git/;/CSV/;/\\\\.svn/;/\\\\.gz/;/\\\\.eps/;/\\\\.pdf/;/\\\\.jpg/;/\\\\.ppm/;/\\\\.log/;/\\\\.avi/;/\\\\.zip/;/\\\\.so/;/\\\\.a/;/\\\\.o/;\\\\.swp$;~$;\\\\.\\\\#;/\\\\#")

INCLUDE(CPack)

INSTALL(FILES ${headers} DESTINATION include/gcop)

# uncomment the following line to get CMake variables to print to screen
# include(CMakePrintSystemInformation)o
