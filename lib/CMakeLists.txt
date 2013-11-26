macro(append_subdir_files variable dirname)
  get_directory_property(holder DIRECTORY ${dirname} DEFINITION ${variable})
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

append_subdir_files(sources systems)
append_subdir_files(headers systems)
append_subdir_files(sources utils)
append_subdir_files(headers utils)
#append_subdir_files(sources devices)
#append_subdir_files(headers devices)
append_subdir_files(sources algos)
append_subdir_files(headers algos)
#append_subdir_files(sources obstacles)
#append_subdir_files(headers obstacles)
append_subdir_files(sources views)
append_subdir_files(headers views)
append_subdir_files(sources est)
append_subdir_files(headers est)

add_library(gcop SHARED ${sources})

INSTALL(TARGETS gcop DESTINATION lib)
INSTALL(FILES ${headers} DESTINATION include/gcop)              

SET_TARGET_PROPERTIES(gcop
                      PROPERTIES
                      VERSION 1.0.0
                      SOVERSION 0)
