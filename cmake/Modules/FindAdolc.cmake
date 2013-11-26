include(LibFindMacros)

# Dependencies
#libfind_package(adolc)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(ADOLC_PKGCONF adolc)

# Include dir
find_path(ADOLC_INCLUDE_DIR
  NAMES adouble.h
  PATHS ${ADOLC_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(ADOLC_LIBRARY
  NAMES adolc
#  PATHS ${ADOLC_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(ADOLC_PROCESS_INCLUDES ADOLC_INCLUDE_DIR)
set(ADOLC_PROCESS_LIBS ADOLC_LIBRARY)
libfind_process(ADOLC)