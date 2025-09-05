cmake_minimum_required(VERSION 3.24)

project(atus-playground)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_FIND_DEBUG_MODE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option( ENABLE_TESTS "3D - in the future this is going to be removed" OFF)
option( BUILD_3D "3D - in the future this is going to be removed" OFF )
option( BUILD_HTRAP "3D" OFF )
option( BUILD_NEHARI "Disables the initial point in function space being set to the Nehari value" ON )
option( BUILD_DOCUMENTATION "Create and install the HTML based API documentation (requires Doxygen)" ${DOXYGEN_FOUND})

if(ENABLE_TESTS)
  enable_testing()
  find_package(Catch2 3 REQUIRED)
  include(CTest)
  include(Catch)
endif()

find_package(deal.II REQUIRED CONFIG)

### nlopt
find_package(NLopt REQUIRED CONFIG)

### Boost
set(Boost_VERBOSE ON)
find_package(Boost 1.82.0 EXACT REQUIRED COMPONENTS program_options)

set( DIR_INC ${PROJECT_SOURCE_DIR}/include CACHE STRING INTERNAL )
set( DIR_INC_2 ${PROJECT_SOURCE_DIR}/source/lib CACHE STRING INTERNAL )

### doxygen
#find_package(Doxygen)

include_directories( ${DIR_INC} ${DIR_INC_2} ) # ${NLOPT_INCLUDE_DIR} ${DEAL_II_MUPARSER_INCLUDE_DIRS} ${DEAL_II_GSL_INCLUDE_DIRS} ${DEAL_II_INCLUDE_DIRS} )

if( BUILD_HTRAP )
  add_definitions( -DPOTENTIAL=1 )
  message(STATUS "Building harmonic trap" )
else()
add_definitions( -DPOTENTIAL=2 )
  message(STATUS "Building gravito optical trap" )
endif()

if( BUILD_3D )
  add_definitions( -DSPATIAL_DIM=3 )
  message(STATUS "Building 3D" )
else()
  add_definitions( -DSPATIAL_DIM=2 )
  message(STATUS "Building 2D" )
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
#add_subdirectory( tests/lib )


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

# serial progs
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
