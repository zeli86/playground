
ADD_EXECUTABLE( gen_Psi_1 gen_Psi_1.cpp )
DEAL_II_SETUP_TARGET( gen_Psi_1 )
TARGET_LINK_LIBRARIES( gen_Psi_1 myatus )

ADD_EXECUTABLE( binR_to_atus2 binR_to_atus2.cpp )
DEAL_II_SETUP_TARGET( binR_to_atus2 )
TARGET_LINK_LIBRARIES( binR_to_atus2 myatus )

ADD_EXECUTABLE( binC_to_atus2 binC_to_atus2.cpp )
DEAL_II_SETUP_TARGET( binC_to_atus2 )
TARGET_LINK_LIBRARIES( binC_to_atus2 myatus )

ADD_EXECUTABLE( binR_to_binC binR_to_binC.cpp )
DEAL_II_SETUP_TARGET( binR_to_binC )
TARGET_LINK_LIBRARIES( binR_to_binC myatus )

ADD_EXECUTABLE( binR_to_binC_1 binR_to_binC_1.cpp )
DEAL_II_SETUP_TARGET( binR_to_binC_1 )
TARGET_LINK_LIBRARIES( binR_to_binC_1 myatus )

#ADD_EXECUTABLE( intp_gamma_mu intp_gamma_mu.cpp )
#TARGET_LINK_LIBRARIES( intp_gamma_mu m gsl gslcblas )

#ADD_EXECUTABLE( virial virial.cpp )
#TARGET_LINK_LIBRARIES( virial ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( virial_cs virial_cs.cpp )
#TARGET_LINK_LIBRARIES( virial_cs ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( compute_lifetimes compute_lifetimes.cpp )
#TARGET_LINK_LIBRARIES( compute_lifetimes ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( compute_lin_op compute_lin_op.cpp )
#TARGET_LINK_LIBRARIES( compute_lin_op ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( test_ef test_ef.cpp )
#TARGET_LINK_LIBRARIES( test_ef ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE(gen_pvd gen_pvd.cpp)
