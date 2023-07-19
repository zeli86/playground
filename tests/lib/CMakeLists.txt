###################################################################################################
set(TARGET TestUtilsRealWavefunction)

add_executable(${TARGET}
  ${TARGET}.cpp 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/source/lib )
target_link_libraries(${TARGET} myatus Catch2::Catch2WithMain)

catch_discover_tests(${TARGET})

unset(TARGET)

###################################################################################################
set(TARGET TestMPIUtilsRealWavefunction)

add_executable(${TARGET}
  ${TARGET}.cpp 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/source/lib )
target_link_libraries(${TARGET} myatus Catch2::Catch2)

set_target_properties(${TARGET} PROPERTIES CROSSCOMPILING_EMULATOR "mpirun;-np;2;-v")

catch_discover_tests(${TARGET})

unset(TARGET)

###################################################################################################
set(TARGET TestGPUtilsRealWavefunction)

add_executable(${TARGET}
  ${TARGET}.cpp 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/source/lib )
target_link_libraries(${TARGET} myatus Catch2::Catch2WithMain)

catch_discover_tests(${TARGET})

unset(TARGET)

###################################################################################################
set(TARGET TestUtilsComplexWavefunction)

add_executable(${TARGET}
  ${TARGET}.cpp 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/source/lib )
target_link_libraries(${TARGET} myatus Catch2::Catch2WithMain)

catch_discover_tests(${TARGET})

unset(TARGET)

###################################################################################################
set(TARGET TestGPUtilsComplexWavefunction)

add_executable(${TARGET}
  ${TARGET}.cpp 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/source/lib )
target_link_libraries(${TARGET} myatus Catch2::Catch2WithMain)

catch_discover_tests(${TARGET})

unset(TARGET)