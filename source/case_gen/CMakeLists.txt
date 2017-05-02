
#ADD_EXECUTABLE( batch_case_generation_elias batch_case_generation_elias.cpp )
#TARGET_LINK_LIBRARIES( batch_case_generation_elias ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( batch_case_generation_xml batch_case_generation_xml.cpp )
TARGET_LINK_LIBRARIES( batch_case_generation_xml ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )

ADD_EXECUTABLE( batch_case_generation_cs_xml batch_case_generation_cs_xml.cpp )
TARGET_LINK_LIBRARIES( batch_case_generation_cs_xml ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} )
