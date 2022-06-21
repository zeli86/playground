cmake_minimum_required(VERSION 3.10)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS ON)
set(CONAN_DISABLE_CHECK_COMPILER ON)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

#execute_process(COMMAND sh -c "module list")

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

find_package(deal.II HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR})
if(NOT ${deal.II_FOUND})
  message(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
endif ()

#
# Are all dependencies fulfilled?
#
if(NOT DEAL_II_WITH_PETSC OR DEAL_II_PETSC_WITH_COMPLEX) # keep in one line
  message(FATAL_ERROR "
Error! This tutorial requires a deal.II library that was configured with the following options:
    DEAL_II_WITH_PETSC = ON
    DEAL_II_PETSC_WITH_COMPLEX = OFF
However, the deal.II library found at ${DEAL_II_PATH} was configured with these options
    DEAL_II_WITH_PETSC = ${DEAL_II_WITH_PETSC}
    DEAL_II_PETSC_WITH_COMPLEX = ${DEAL_II_PETSC_WITH_COMPLEX}
which conflict with the requirements."
    )
endif()

DEAL_II_INITIALIZE_CACHED_VARIABLES()

include(${DEAL_II_FEATURE_CONFIG})

message("${DEAL_II_FEATURE_CONFIG}")

project(atus-playground)

### nlopt
find_package(NLopt REQUIRED)

### doxygen
find_package(Doxygen)

### Boost
set(Boost_VERBOSE ON)
find_package(Boost REQUIRED COMPONENTS log log_setup program_options filesystem HINTS ${DEAL_II_BOOST_DIR})
#include_directories( ${DEAL_II_BOOST_INCLUDE_DIRS} )

include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
add_definitions("-DGIT_SHA1=${GIT_SHA1}")

option( BUILD_3D "3D - in the future this is going to be removed" OFF )
option( BUILD_HTRAP "3D" OFF )
option( BUILD_NEHARI "Disables the initial point in function space being set to the Nehari value" ON )
option( BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

set( HOME $ENV{HOME} CACHE STRING INTERNAL )
set( DIR_INC ${PROJECT_SOURCE_DIR}/include CACHE STRING INTERNAL )
set( DIR_INC_2 ${PROJECT_SOURCE_DIR}/source/lib CACHE STRING INTERNAL )

set( EXECUTABLE_OUTPUT_PATH  ${HOME}/out/bin )
set( LIBRARY_OUTPUT_PATH  ${HOME}/out/bin/lib )
#SET( EXECUTABLE_INSTALL_PATH  ${HOME}/bin )

# set include dirs
include_directories( ${DIR_INC} ${DIR_INC_2} ${NLOPT_INCLUDE_DIR} ${DEAL_II_MUPARSER_INCLUDE_DIRS} ${DEAL_II_GSL_INCLUDE_DIRS} ${DEAL_II_INCLUDE_DIRS} )

if( BUILD_HTRAP )
  add_definitions( -DPOTENTIAL=1 )
  message(STATUS "*** Building harmonic trap" )
else()
add_definitions( -DPOTENTIAL=2 )
  message(STATUS "*** Building gravito optical trap" )
endif()

if( BUILD_3D )
  add_definitions( -DSPATIAL_DIM=3 )
  message(STATUS "*** Building 3D" )
else()
  add_definitions( -DSPATIAL_DIM=2 )
  message(STATUS "*** Building 2D" )
endif()

if( BUILD_NEHARI )
  add_definitions( -DNEHARI )
endif()

add_subdirectory( source/lib )
add_subdirectory( source/case_gen )
add_subdirectory( source/stationary )
#add_subdirectory( source/groundstates )
#add_subdirectory( source/realtime )
#add_subdirectory( source/oct )
#add_subdirectory( source/utils )
#add_subdirectory( source/edu )

if(BUILD_DOCUMENTATION)
if(NOT DOXYGEN_FOUND)
    message(FATAL_ERROR "Doxygen is needed to build the documentation.")
endif()

set(doxyfile ${PROJECT_SOURCE_DIR}/Doxyfile)

add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${doxyfile}
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Generating API documentation with Doxygen"
    VERBATIM)
endif()

# enable_testing()
# serial progs
#add_test( NAME CASEGENTEST1 COMMAND ${EXECUTABLE_OUTPUT_PATH}/batch_case_generation_xml --dim 1 -f test1d -N 2 --Ndmu 2)
#add_test( NAME CASEGENTEST2 COMMAND ${EXECUTABLE_OUTPUT_PATH}/batch_case_generation_xml --dim 2 -f test2d -N 1 --Ndmu 2)
#add_test( NAME BREED1_0 COMMAND ${EXECUTABLE_OUTPUT_PATH}/breed_1 WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test1d/0/)
#add_test( NAME BREED1_1 COMMAND ${EXECUTABLE_OUTPUT_PATH}/breed_1 WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test1d/1/)
#add_test( NAME BREEDSOB1_0 COMMAND ${EXECUTABLE_OUTPUT_PATH}/breed_sob_1 WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test1d/0/0001)

#execute_process(COMMAND sh -c "mkdir -p rt_prop_1_demo")
#execute_process(COMMAND sh -c "cp ${CMAKE_CURRENT_SOURCE_DIR}/xml/rt_prop_1_demo.xml ${CMAKE_CACHEFILE_DIR}/rt_prop_1_demo/")
#add_test( NAME RTPROP0 COMMAND ${EXECUTABLE_OUTPUT_PATH}/gen_Psi_1 -p rt_prop_1_demo.xml WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/rt_prop_1_demo)
#add_test( NAME RTPROP1 COMMAND ${EXECUTABLE_OUTPUT_PATH}/rt_prop -p rt_prop_1_demo.xml -i Psi0.1d.bin WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/rt_prop_1_demo)

#execute_process(COMMAND sh -c "mkdir -p rt_prop_2_demo")
#execute_process(COMMAND sh -c "cp ${CMAKE_CURRENT_SOURCE_DIR}/xml/rt_prop_1_demo.xml ${CMAKE_CACHEFILE_DIR}/rt_prop_2_demo/")
#add_test( NAME RTPROP2 COMMAND ${EXECUTABLE_OUTPUT_PATH}/gen_Psi_1 -p rt_prop_1_demo.xml WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/rt_prop_2_demo)
#add_test( NAME RTPROP3 COMMAND ${EXECUTABLE_OUTPUT_PATH}/rt_prop_2 -p rt_prop_1_demo.xml -i Psi0.1d.bin WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/rt_prop_2_demo)

# mpi progs
#set (test_parameters -np ${MPIEXEC_MAX_NUMPROCS} ${EXECUTABLE_OUTPUT_PATH}/breed_mpi)
#add_test( NAME BREEDMPI_0 COMMAND ${MPIEXEC} ${test_parameters} WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test2d/0_0/)
#add_test( NAME BREEDMPI_1 COMMAND ${MPIEXEC} ${test_parameters} WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test2d/1_0/)

FILE(GLOB bash_sh "${PROJECT_SOURCE_DIR}/bash/*")
FOREACH( file_i ${bash_sh})
    MESSAGE(STATUS ${file_i} )
    INSTALL(FILES ${file_i} PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION "${HOME}/bin" )
ENDFOREACH( file_i )
