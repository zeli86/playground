


ADD_EXECUTABLE( prep_oct_gs prep_oct_gs_from_es.cpp )
TARGET_LINK_LIBRARIES( prep_oct_gs ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( prep_oct_transport prep_oct_transport.cpp )
#TARGET_LINK_LIBRARIES( prep_oct_transport ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( breed_4_oct_1 breed_4_oct_1.cpp )
#TARGET_LINK_LIBRARIES( breed_4_oct_1 ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )
