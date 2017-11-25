
# Find the NLOPT library
# Usage:
#   find_package( NLOPT [REQUIRED] [QUIET] )
#

if( NLOPT_ROOT )

    find_library(  NLOPT_LIBRARY
                   NAMES "nlopt" 
                   PATHS ${NLOPT_ROOT}
                   PATH_SUFFIXES "lib" "lib64"
                   NO_DEFAULT_PATH
    )

    find_path(  NLOPT_INCLUDE_DIR
                NAMES "nlopt.h"
                PATHS ${NLOPT_ROOT}
                PATH_SUFFIXES "include"
                NO_DEFAULT_PATH
    )

else()

    find_library(   NLOPT_LIBRARY
                    NAMES "nlopt"  
                    PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH 
    )

    get_filename_component( TMP ${NLOPT_LIBRARY} PATH )
    get_filename_component( TMP ${TMP} PATH )
    set( NLOPT_INCLUDE_DIR ${TMP}/include CACHE STRING INTERNAL )

endif()

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(NLOPT
      REQUIRED_VARS NLOPT_INCLUDE_DIR NLOPT_LIBRARY
      HANDLE_COMPONENTS
      )

mark_as_advanced(
      NLOPT_LIBRARY
      NLOPT_INCLUDE_DIR
      )
