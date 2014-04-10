include(LibFindMacros)

# Dependencies
#libfind_package(eigen)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(EIGEN_PKGCONF eigen3)

# Include dir
find_path(EIGEN_INCLUDE_DIR
  NAMES EigenBase.h
  PATHS ${EIGEN_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
#find_library(EIGEN_LIBRARY
#  NAMES eigen
#  PATHS ${EIGEN_PKGCONF_LIBRARY_DIRS}
#)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(EIGEN_PROCESS_INCLUDES EIGEN_INCLUDE_DIR)
#set(EIGEN_PROCESS_LIBS _LIBRARY)
#libfind_process(GCOP)