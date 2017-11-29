
ADD_EXECUTABLE( breed_minj_mpi breed_minj_mpi.cpp )
TARGET_LINK_LIBRARIES( breed_minj_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breed_minj_cs_mpi breed_minj_cs_mpi.cpp )
TARGET_LINK_LIBRARIES( breed_minj_cs_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breed_sob_mpi breed_sob_mpi.cpp )
TARGET_LINK_LIBRARIES( breed_sob_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( my_newton_2_mpi my_newton_2_mpi.cpp )
TARGET_LINK_LIBRARIES( my_newton_2_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( my_schroedinger_newton_v0_mpi my_schroedinger_newton_v0_mpi.cpp )
#TARGET_LINK_LIBRARIES( my_schroedinger_newton_v0_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( my_schroedinger_newton_cs_mpi my_schroedinger_newton_cs_mpi.cpp )
#TARGET_LINK_LIBRARIES( my_schroedinger_newton_cs_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )
