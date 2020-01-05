
ADD_EXECUTABLE( edu_mpi edu_mpi.cpp )
DEAL_II_SETUP_TARGET( edu_mpi )
TARGET_LINK_LIBRARIES( edu_mpi myatus )
