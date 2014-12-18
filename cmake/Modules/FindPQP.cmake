include(LibFindMacros)

# Dependencies
#libfind_package(dsl)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(PQP_PKGCONF PQP)

# Include dir
find_path(PQP_INCLUDE_DIR
  NAMES PQP/PQP.h
  PATHS ${PQP_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(PQP_LIBRARY
  NAMES PQP
  PATHS ${PQP_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(PQP_PROCESS_INCLUDES PQP_INCLUDE_DIR)
set(PQP_PROCESS_LIBS PQP_LIBRARY)
libfind_process(PQP)