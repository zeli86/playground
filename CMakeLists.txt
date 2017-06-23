PROJECT (doktorarbeit)

# specify compiler
# cmake -DCMAKE_CXX_COMPILER=g++-6 -DCMAKE_C_COMPILER=gcc-6 ..

cmake_minimum_required(VERSION 2.8)
find_package(MPI REQUIRED)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

option( BUILD_3D "3D" OFF )
option( BUILD_HTRAP "3D" OFF )
option( BUILD_VARIANT2 "Variant 2" ON )
option( BUILD_NEHARI "Disables the initial point in function space being set to the Nehari value" ON )

set(CMAKE_CXX_FLAGS_RELEASE "-std=gnu++11 -march=native -O3 -fopenmp -mieee-fp -funroll-loops -ftree-vectorize -w -s")
set(CMAKE_CXX_FLAGS_DEBUG "-std=gnu++11 -g -march=native -fopenmp -mieee-fp ")

IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Debug Release." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Release" "Debug" )
ENDIF(NOT CMAKE_BUILD_TYPE)

set( HOME $ENV{HOME} CACHE STRING INTERNAL )
set( DIR_INC ${PROJECT_SOURCE_DIR}/include CACHE STRING INTERNAL )
set( DIR_INC_2 ${PROJECT_SOURCE_DIR}/source/lib CACHE STRING INTERNAL )

### NLOPT library paths
#find_library( NLOPT_LIBRARY NAMES nlopt PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
#if( NOT NLOPT_LIBRARY )
#  message( FATAL_ERROR "Could not find the nlopt library." )
#endif()
#get_filename_component( TMP ${NLOPT_LIBRARY} PATH )
#get_filename_component( TMP ${TMP} PATH )
#set( NLOPT_INC_DIR ${TMP}/include CACHE STRING INTERNAL )
#if( NLOPT_LIBRARY )
#  include_directories( ${NLOPT_INC_DIR} )
#endif()

### GSL library paths
find_library( GSL_LIBRARY_1 NAMES gsl  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT GSL_LIBRARY_1 )
  message( FATAL_ERROR "Could not find the gsl library." )
endif()
find_library( GSL_LIBRARY_2 NAMES gslcblas PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT GSL_LIBRARY_2 )
  message( FATAL_ERROR "Could not find the gsl library." )
endif()
get_filename_component( TMP ${GSL_LIBRARY_1} PATH )
get_filename_component( TMP ${TMP} PATH )
set( GSL_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

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

### SLEPc library paths
#find_library( SLEPC_LIBRARY NAMES slepc  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
#if( NOT SLEPC_LIBRARY )
#  message( FATAL_ERROR "Could not find the SLEPc library." )
#endif()
#get_filename_component( TMP ${SLEPC_LIBRARY} PATH )
#get_filename_component( TMP ${TMP} PATH )
#set( SLEPC_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

# muParser library
find_library( MUPARSER_LIBRARY NAMES muparser  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT MUPARSER_LIBRARY )
  message( FATAL_ERROR "Could not find libmuparser." )
endif()
get_filename_component( TMP ${MUPARSER_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( MUPARSER_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

### deal.ii library paths
find_library( DEAL_II_LIBRARY NAMES deal_II  PATHS ENV LD_LIBRARY_PATH NO_DEFAULT_PATH )
if( NOT DEAL_II_LIBRARY )
  message( FATAL_ERROR "Could not find the deal.ii library." )
endif()
get_filename_component( TMP ${DEAL_II_LIBRARY} PATH )
get_filename_component( TMP ${TMP} PATH )
set( DEAL_INC_DIR ${TMP}/include CACHE STRING INTERNAL )

#SET( BOOST_INC_DIR $ENV{MY_PATH_BOOST}/include CACHE STRING INTERNAL )
#SET( BOOST_LIB_DIR $ENV{MY_PATH_BOOST}/lib CACHE STRING INTERNAL )

SET( EXECUTABLE_OUTPUT_PATH  ${HOME}/bin )
#SET( LIBRARY_OUTPUT_PATH  ${HOME}/bin/lib )
#SET( EXECUTABLE_INSTALL_PATH  ${HOME}/bin )

# Zuerst die Pfade angeben
#=========================
# ${DEAL_INC_DIR}/deal.II/bundled ${BOOST_INC_DIR}
include_directories( ${MPI_CXX_INCLUDE_PATH} ${DIR_INC} ${DIR_INC_2} ${DEAL_INC_DIR} ${DEAL_INC_DIR}/deal.II/bundled ${GSL_INC_DIR} ${P4EST_INC_DIR} ${PETSC_INC_DIR} ${MUPARSER_INC_DIR} )

#LINK_DIRECTORIES( ${DEAL_LIB_DIR} ${GSL_LIB_DIR} ${P4EST_LIB_DIR} ${PETSC_LIB_DIR} ${BOOST_LIB_DIR} )

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")
include(GetGitRevisionDescription)
get_git_head_revision(GIT_REFSPEC GIT_SHA1)
add_definitions("-DGIT_SHA1=${GIT_SHA1}")

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

IF( BUILD_VARIANT2 )
  ADD_DEFINITIONS( -D__variant_2__ )
  MESSAGE(STATUS "*** Building __variant_2__" )
ELSE()
  ADD_DEFINITIONS( -D__variant_1__ )
  MESSAGE(STATUS "*** Building __variant_1__" )  
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
#ADD_SUBDIRECTORY( source/trash )
#ADD_SUBDIRECTORY( source/belen )
ADD_SUBDIRECTORY( source/cont )

FILE(GLOB bash_sh "${PROJECT_SOURCE_DIR}/bash/*")

FOREACH( file_i ${bash_sh})
    MESSAGE(STATUS ${file_i} )
    INSTALL(FILES ${file_i} PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ DESTINATION "${HOME}/bin" )
ENDFOREACH( file_i )
