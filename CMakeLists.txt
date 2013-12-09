cmake_minimum_required (VERSION 2.6) 
project (gcop)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)

ADD_DEFINITIONS(-g -O -std=c++0x -pthread)

SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Robotic Planning and Control")
SET(CPACK_PACKAGE_VENDOR "Marin Kobilarov")
#SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/ReadMe.txt")
#SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "CMake ${CMake_VERSION_MAJOR}.${CMake_VERSION_MINOR}")

#SET(CPACK_STRIP_FILES "bin/MyExecutable")
#SET(CPACK_SOURCE_STRIP_FILES "")

#find_package(OpenGL REQUIRED)
#find_package(GLUT REQUIRED)
#include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
#target_link_libraries( ${OpenGL_LIBRARIES} ${GLUT_LIBRARY} )


find_package(GLU REQUIRED)
include_directories(${GLU_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GLU_LIBRARIES} ${GLU_LIBRARY})

find_package(GLUT REQUIRED)
include_directories(${GLUT_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GLUT_LIBRARY})

find_package(OpenGL REQUIRED)
include_directories(${OpenGL_INCLUDE_DIRS})
set(LIBS ${LIBS} ${OpenGL_LIBRARIES})

find_package(Adolc REQUIRED)
include_directories(${Adolc_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Adolc_LIBRARIES})



find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
set(LIBS ${LIBS} ${Boost_LIBRARY})


find_package (Threads)
set (LIBS ${LIBS} ${CMAKE_THREAD_LIBS_INIT})

find_package(Eigen REQUIRED)
include_directories(${Eigen_INCLUDE_DIRS})

find_package(TinyXML REQUIRED)
include_directories(${TINYXML_INCLUDE_DIR})
set(LIBS ${LIBS} ${TINYXML_LIBRARIES})

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
"/CSV/;/\\\\.svn/;/\\\\.gz/;/\\\\.eps/;/\\\\.pdf/;/\\\\.jpg/;/\\\\.ppm/;/\\\\.log/;/\\\\.avi/;/\\\\.zip/;/\\\\.so/;/\\\\.a/;/\\\\.o/;\\\\.swp$;~$;\\\\.\\\\#;/\\\\#")

INCLUDE(CPack)

INSTALL(FILES ${headers} DESTINATION include/gcop)

# uncomment the following line to get CMake variables to print to screen
# include(CMakePrintSystemInformation)o
