
ADD_EXECUTABLE( quantum_mirror quantum_mirror.cpp )
TARGET_LINK_LIBRARIES( quantum_mirror ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )
