macro(append_subdir_files variable variable2 dirname)
  get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable2})
  foreach(depfile ${holder})
    list(APPEND ${variable} "${dirname}/${depfile}")
  endforeach()
endmacro()

add_subdirectory (algos)
#add_subdirectory (devices)
#add_subdirectory (obstacles) 
add_subdirectory (systems) 
add_subdirectory (utils)
add_subdirectory (views)
add_subdirectory (est)

#SOURCES and HEADERS are divided for 
# systems viewers algos utils est

append_subdir_files(sys_sources sources systems)
append_subdir_files(util_sources sources utils)
append_subdir_files(algo_sources sources algos)
append_subdir_files(viewer_sources sources views)
append_subdir_files(est_sources sources est)

#append_subdir_files(sources sources devices)
#append_subdir_files(headers headers devices)
#append_subdir_files(sources obstacles)
#append_subdir_files(headers obstacles)

#append_subdir_files(sys_headers headers systems)
#append_subdir_files(util_headers headers utils)
#append_subdir_files(est_headers headers est)
#append_subdir_files(viewer_headers headers views)
#append_subdir_files(algo_headers headers algos)

append_subdir_files(headers headers systems)
append_subdir_files(headers headers utils)
append_subdir_files(headers headers views)
append_subdir_files(headers headers est)
append_subdir_files(headers headers algos)

# Creating separate libraries for each of the subdirectories 
# they will be treated as separate modules
#Util Headers:
#INSTALL(FILES ${util_headers} DESTINATION include/gcop/utils)              
INSTALL(FILES ${headers} DESTINATION include/gcop)#This will be changed later TODO              
#Utils
add_library(gcop_utils  SHARED ${util_sources})
INSTALL(TARGETS gcop_utils DESTINATION lib)

SET_TARGET_PROPERTIES(gcop_utils
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 0)

#Systems:
add_library(gcop_systems SHARED ${sys_sources})
#target_link_libraries(gcop_systems gcop_utils)#Dependency of modules

INSTALL(TARGETS gcop_systems DESTINATION lib)
#INSTALL(FILES ${sys_headers} DESTINATION include/gcop/systems)              

SET_TARGET_PROPERTIES(gcop_systems
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 0)
# Viewers:
add_library(gcop_views SHARED ${viewer_sources})
#target_link_libraries(gcop_views gcop_utils)#Dependency of modules

INSTALL(TARGETS gcop_views DESTINATION lib)
#INSTALL(FILES ${viewer_headers} DESTINATION include/gcop/views)              

SET_TARGET_PROPERTIES(gcop_views
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 0)
#Algorithms
add_library(gcop_algos SHARED ${algo_sources})
#target_link_libraries(gcop_algos gcop_utils)#Dependency of modules

INSTALL(TARGETS gcop_algos DESTINATION lib)
#INSTALL(FILES ${algo_headers} DESTINATION include/gcop/algos)              

SET_TARGET_PROPERTIES(gcop_algos
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 0)
# Estimation:
add_library(gcop_est SHARED ${est_sources})
#target_link_libraries(gcop_est gcop_utils)#Dependency of modules

INSTALL(TARGETS gcop_est DESTINATION lib)
#INSTALL(FILES ${est_headers} DESTINATION include/gcop/systems)              

SET_TARGET_PROPERTIES(gcop_est
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 0)
