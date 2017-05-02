
ADD_EXECUTABLE( rt_prop_mpi rt_prop_mpi.cpp )
TARGET_LINK_LIBRARIES( rt_prop_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( rt_prop_mpi_2 rt_prop_mpi_2.cpp )
#TARGET_LINK_LIBRARIES( rt_prop_mpi_2 ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( rt_prop_mpi_cs rt_prop_mpi_cs.cpp )
TARGET_LINK_LIBRARIES( rt_prop_mpi_cs ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( rt_prop_mpi_cs_cmplx_pot rt_prop_mpi_cs_cmplx_pot.cpp )
#TARGET_LINK_LIBRARIES( rt_prop_mpi_cs_cmplx_pot ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( rt_prop rtprop_1.cpp )
TARGET_LINK_LIBRARIES( rt_prop ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( rt_prop_err_mpi_cs rt_prop_err_mpi_cs.cpp )
#TARGET_LINK_LIBRARIES( rt_prop_err_mpi_cs ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )
