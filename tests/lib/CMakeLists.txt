###################################################################################################
set(TARGET TestUtilsRealwavefunction)

add_executable(${TARGET}
  ${TARGET}.cpp 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/source/lib )
target_link_libraries(${TARGET} myatus Catch2::Catch2WithMain)

catch_discover_tests(${TARGET})

unset(TARGET)

###################################################################################################
set(TARGET TestGPUtilsRealwavefunction)

add_executable(${TARGET}
  ${TARGET}.cpp 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/source/lib )
target_link_libraries(${TARGET} myatus Catch2::Catch2WithMain)

catch_discover_tests(${TARGET})

unset(TARGET)