ADD_EXECUTABLE( breed_mpi my_newton_mpi.cpp )
TARGET_LINK_LIBRARIES( breed_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} ${NLOPT_LIBRARY} ${MUPARSER_LIBRARY})

ADD_EXECUTABLE( breed_cs_mpi my_newton_cs_mpi.cpp )
TARGET_LINK_LIBRARIES( breed_cs_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} ${NLOPT_LIBRARY} )

ADD_EXECUTABLE( breed_minj_mpi breed_minj_mpi.cpp )
TARGET_LINK_LIBRARIES( breed_minj_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breed_minj_cs_mpi breed_minj_cs_mpi.cpp )
TARGET_LINK_LIBRARIES( breed_minj_cs_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breed_sob_mpi breed_sob_mpi.cpp )
TARGET_LINK_LIBRARIES( breed_sob_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( breed_1 breed_1.cpp )
TARGET_LINK_LIBRARIES( breed_1 ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} ${NLOPT_LIBRARY} )

ADD_EXECUTABLE( breed_sob_1 breed_sob_1.cpp )
TARGET_LINK_LIBRARIES( breed_sob_1 ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} ${NLOPT_LIBRARY} )
