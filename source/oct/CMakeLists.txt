
ADD_EXECUTABLE( oct oct.cpp )
DEAL_II_SETUP_TARGET( oct )
TARGET_LINK_LIBRARIES( oct myatus )

#ADD_EXECUTABLE( oct_mpi oct_mpi.cpp )
#DEAL_II_SETUP_TARGET( oct )
#TARGET_LINK_LIBRARIES( oct_mpi ${MPI_CXX_LIBRARIES} myatus ${GSL_LIBRARY_1} ${GSL_LIBRARY_2} m ${DEAL_II_LIBRARY} ${PETSC_LIBRARY} ${MUPARSER_LIBRARY} )

ADD_EXECUTABLE( prep_testcase prep_testcase.cpp )
DEAL_II_SETUP_TARGET( prep_testcase )
TARGET_LINK_LIBRARIES( prep_testcase myatus )
