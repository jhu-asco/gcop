include(LibFindMacros)

# Dependencies
#libfind_package(gcop)

# Use pkg-config to get hints about paths
libfind_pkg_check_modules(GCOP_PKGCONF gcop)

# Include dir
find_path(GCOP_INCLUDE_DIR
  NAMES gcop/system.h
  PATHS ${GCOP_PKGCONF_INCLUDE_DIRS}
)

# Finally the library itself
find_library(GCOP_LIBRARY
  NAMES gcop
#  PATHS ${GCOP_PKGCONF_LIBRARY_DIRS}
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
set(GCOP_PROCESS_INCLUDES GCOP_INCLUDE_DIR)
set(GCOP_PROCESS_LIBS GCOP_LIBRARY)
libfind_process(GCOP)

