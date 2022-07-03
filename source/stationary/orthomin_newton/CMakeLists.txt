
set(TARGET orthomin_newton_mpi )

add_executable(${TARGET}
  ${TARGET}_main.cpp 
  ${TARGET}.hpp 
  ${TARGET}.cpp 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

target_include_directories(${TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include )
target_link_libraries(${TARGET} myatus NLopt::nlopt)

unset(${TARGET})