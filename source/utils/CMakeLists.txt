
ADD_EXECUTABLE( binR_to_atus2 binR_to_atus2.cpp )
TARGET_LINK_LIBRARIES( binR_to_atus2 ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( binC_to_atus2 binC_to_atus2.cpp )
TARGET_LINK_LIBRARIES( binC_to_atus2 ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( gen_qsubs gen_qsubs.cpp )
#TARGET_LINK_LIBRARIES( gen_qsubs ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( gen_hlrn gen_hlrn.cpp )
#TARGET_LINK_LIBRARIES( gen_hlrn ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( gen_bash gen_bash.cpp )
#TARGET_LINK_LIBRARIES( gen_bash ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( gen_slurm gen_slurm.cpp )
#TARGET_LINK_LIBRARIES( gen_slurm ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( intp_gamma_mu intp_gamma_mu.cpp )
#TARGET_LINK_LIBRARIES( intp_gamma_mu m gsl gslcblas )

#ADD_EXECUTABLE( bin_2_vtu bin_2_vtu.cpp )
#TARGET_LINK_LIBRARIES( bin_2_vtu ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( virial virial.cpp )
#TARGET_LINK_LIBRARIES( virial ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( virial_cs virial_cs.cpp )
#TARGET_LINK_LIBRARIES( virial_cs ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( compute_lifetimes compute_lifetimes.cpp )
#TARGET_LINK_LIBRARIES( compute_lifetimes ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( compute_lin_op compute_lin_op.cpp )
#TARGET_LINK_LIBRARIES( compute_lin_op ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( compute_inf_norm_atus2_wf compute_inf_norm_atus2_wf.cpp )
#TARGET_LINK_LIBRARIES( compute_inf_norm_atus2_wf ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

#ADD_EXECUTABLE( test_ef test_ef.cpp )
#TARGET_LINK_LIBRARIES( test_ef ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE(gen_pvd gen_pvd.cpp)
