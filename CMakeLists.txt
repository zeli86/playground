PROJECT (atus-playground)
cmake_minimum_required(VERSION 3.3)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

message("**********************************************************************")
execute_process(COMMAND sh -c "module list")
message("**********************************************************************")

find_package(MPI REQUIRED)
find_package(GSL REQUIRED)
find_package(NLOPT REQUIRED)
find_package(MUPARSER REQUIRED)
find_package(Doxygen)
#find_package(Boost REQUIRED)
message("**********************************************************************")


include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
add_definitions("-DGIT_SHA1=${GIT_SHA1}")

#include(ExternalProject)

#ExternalProject_add( CXXOPT
#                     PREFIX ${CMAKE_CACHEFILE_DIR}/ThirdParty
#                     GIT_REPOSITORY https://github.com/jarro2783/cxxopts
#                     INSTALL_COMMAND ""
#)

#file(SHA256 ${PROJECT_SOURCE_DIR}/include/cxxopts.hpp CXXOPTSHA)
#file(SHA256 ${CMAKE_CACHEFILE_DIR}/ThirdParty/src/CXXOPT/include/cxxopts.hpp CXXOPTGITSHA)

#if(NOT CXXOPTSHA EQUAL CXXOPTGITSHA )
#  file( COPY ${CMAKE_CACHEFILE_DIR}/ThirdParty/src/CXXOPT/include/cxxopts.hpp DESTINATION ${PROJECT_SOURCE_DIR}/include/)
#endif()

#ExternalProject_add( PUGIXML
#                     PREFIX ${CMAKE_CACHEFILE_DIR}/ThirdParty
#                     GIT_REPOSITORY https://github.com/zeux/pugixml
#                     INSTALL_COMMAND ""
#)


option( BUILD_3D "3D - in the future this is going to be removed" OFF )
option( BUILD_HTRAP "3D" OFF )
option( BUILD_NEHARI "Disables the initial point in function space being set to the Nehari value" ON )
option( BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

set(CMAKE_CXX_FLAGS_RELEASE "-std=gnu++11 -march=native -O3 -fopenmp -mieee-fp -funroll-loops -ftree-vectorize -w -s")
set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++11 -g -march=native -fopenmp -mieee-fp -fsanitize=address")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" )
ENDIF(NOT CMAKE_BUILD_TYPE)

set( HOME $ENV{HOME} CACHE STRING INTERNAL )
set( DIR_INC ${PROJECT_SOURCE_DIR}/include CACHE STRING INTERNAL )
set( DIR_INC_2 ${PROJECT_SOURCE_DIR}/source/lib CACHE STRING INTERNAL )

### p4est library paths
find_library( P4EST_LIBRARY NAMES p4est PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT P4EST_LIBRARY )
  message( FATAL_ERROR "Could not find the p4est library." )
endif()
get_filename_component( TMP ${P4EST_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( P4EST_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

### PETSc library paths
find_library( PETSC_LIBRARY NAMES petsc  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT PETSC_LIBRARY )
  message( FATAL_ERROR "Could not find the PETSc library." )
endif()
get_filename_component( TMP ${PETSC_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( PETSC_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

### deal.ii library paths
find_library( DEAL_II_LIBRARY NAMES deal_II  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT DEAL_II_LIBRARY )
  message( FATAL_ERROR "Could not find the deal.ii library." )
endif()
get_filename_component( TMP ${DEAL_II_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( DEAL_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

SET( EXECUTABLE_OUTPUT_PATH  ${HOME}/bin )
#SET( LIBRARY_OUTPUT_PATH  ${HOME}/bin/lib )
#SET( EXECUTABLE_INSTALL_PATH  ${HOME}/bin )

# set include dirs
include_directories(  ${MPI_CXX_INCLUDE_PATH} 
                      ${GSL_INCLUDE_DIR} 
                      ${NLOPT_INCLUDE_DIR} 
                      ${MUPARSER_INCLUDE_DIR}
                      ${DIR_INC} 
                      ${DIR_INC_2} 
                      ${PETSC_INC_DIR} 
                      ${P4EST_INC_DIR} 
                      ${DEAL_INC_DIR} 
                      ${DEAL_INC_DIR}/deal.II/bundled 
)

IF( BUILD_HTRAP )
  ADD_DEFINITIONS( -DPOTENTIAL=1 )
  MESSAGE(STATUS "*** Building harmonic trap" )
ELSE()
  ADD_DEFINITIONS( -DPOTENTIAL=2 )
  MESSAGE(STATUS "*** Building gravito optical trap" )
ENDIF()

IF( BUILD_3D )
  ADD_DEFINITIONS( -DDIMENSION=3 )
  MESSAGE(STATUS "*** Building 3D" )
ELSE()
  ADD_DEFINITIONS( -DDIMENSION=2 )
  MESSAGE(STATUS "*** Building 2D" )
ENDIF()

if( BUILD_NEHARI )
  add_definitions( -DNEHARI )
endif()

ADD_SUBDIRECTORY( source/lib )
ADD_SUBDIRECTORY( source/groundstates )
ADD_SUBDIRECTORY( source/realtime )
ADD_SUBDIRECTORY( source/oct )
ADD_SUBDIRECTORY( source/utils )
ADD_SUBDIRECTORY( source/case_gen )
ADD_SUBDIRECTORY( source/stationary )
#ADD_SUBDIRECTORY( source/resonances )
ADD_SUBDIRECTORY( source/trash )
#ADD_SUBDIRECTORY( source/belen )
ADD_SUBDIRECTORY( source/edu )
#ADD_SUBDIRECTORY( source/cont )

#set(cpuinfo_file "/proc/cpuinfo")
#if(EXISTS "${cpuinfo_file}")
#  file(STRINGS "${cpuinfo_file}" procs REGEX "^processor.: [0-9]+$")
#  list(LENGTH procs PROCESSOR_COUNT)
#endif()

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

message(STATUS "number of cores: ${MPIEXEC_MAX_NUMPROCS}\n")

enable_testing()
# serial progs
add_test( NAME CASEGENTEST1 COMMAND ${EXECUTABLE_OUTPUT_PATH}/batch_case_generation_xml --dim 1 -f test1d -N 2 --Ndmu 2)
add_test( NAME CASEGENTEST2 COMMAND ${EXECUTABLE_OUTPUT_PATH}/batch_case_generation_xml --dim 2 -f test2d -N 1 --Ndmu 2)
add_test( NAME BREED1_0 COMMAND ${EXECUTABLE_OUTPUT_PATH}/breed_1 WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test1d/0/)
add_test( NAME BREED1_1 COMMAND ${EXECUTABLE_OUTPUT_PATH}/breed_1 WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test1d/1/)
add_test( NAME BREEDSOB1_0 COMMAND ${EXECUTABLE_OUTPUT_PATH}/breed_sob_1 WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test1d/0/0001)

execute_process(COMMAND sh -c "mkdir -p rt_prop_1_demo")
execute_process(COMMAND sh -c "cp ${CMAKE_CURRENT_SOURCE_DIR}/xml/rt_prop_1_demo.xml ${CMAKE_CACHEFILE_DIR}/rt_prop_1_demo")
add_test( NAME RTPROP0 COMMAND ${EXECUTABLE_OUTPUT_PATH}/gen_Psi_1 -p rt_prop_1_demo.xml WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/rt_prop_1_demo)
add_test( NAME RTPROP1 COMMAND ${EXECUTABLE_OUTPUT_PATH}/rt_prop -p rt_prop_1_demo.xml -i Psi0.1d.bin WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/rt_prop_1_demo)
 
# mpi progs
set (test_parameters -np ${MPIEXEC_MAX_NUMPROCS} ${EXECUTABLE_OUTPUT_PATH}/breed_mpi)
add_test( NAME BREEDMPI_0 COMMAND ${MPIEXEC} ${test_parameters} WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test2d/0_0/)
add_test( NAME BREEDMPI_1 COMMAND ${MPIEXEC} ${test_parameters} WORKING_DIRECTORY ${CMAKE_CACHEFILE_DIR}/test2d/1_0/)

FILE(GLOB bash_sh "${PROJECT_SOURCE_DIR}/bash/*")
FOREACH( file_i ${bash_sh})
    MESSAGE(STATUS ${file_i} )
    INSTALL(FILES ${file_i} PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION "${HOME}/bin" )
ENDFOREACH( file_i )
