ADD_EXECUTABLE( breed_mpi my_newton_mpi.cpp )
DEAL_II_SETUP_TARGET( breed_mpi )
TARGET_LINK_LIBRARIES( breed_mpi myatus NLopt::nlopt)

# ADD_EXECUTABLE( breed_cs_mpi my_newton_cs_mpi.cpp )
# DEAL_II_SETUP_TARGET( breed_cs_mpi )
# TARGET_LINK_LIBRARIES( breed_cs_mpi myatus NLopt::nlopt)

#ADD_EXECUTABLE( breed_minj_mpi breed_minj_mpi.cpp )
#DEAL_II_SETUP_TARGET( breed_minj_mpi )
#TARGET_LINK_LIBRARIES( breed_minj_mpi myatus NLopt::nlopt)

# ADD_EXECUTABLE( breed_minj_cs_mpi breed_minj_cs_mpi.cpp )
# DEAL_II_SETUP_TARGET( breed_minj_cs_mpi )
# TARGET_LINK_LIBRARIES( breed_minj_cs_mpi myatus )

#ADD_EXECUTABLE( rbreed_mpi rbreed_minj_mpi.cpp )
#DEAL_II_SETUP_TARGET( rbreed_mpi )
#TARGET_LINK_LIBRARIES( rbreed_mpi myatus NLopt::nlopt)

#ADD_EXECUTABLE( breed_sob_mpi breed_sob_mpi.cpp )
#DEAL_II_SETUP_TARGET( breed_sob_mpi )
#TARGET_LINK_LIBRARIES( breed_sob_mpi myatus )

# ADD_EXECUTABLE( breed_1 breed_1.cpp )
# DEAL_II_SETUP_TARGET( breed_1 )
# TARGET_LINK_LIBRARIES( breed_1 myatus NLopt::nlopt )

# ADD_EXECUTABLE( breed_sob_1 breed_sob_1.cpp )
# DEAL_II_SETUP_TARGET( breed_sob_1 )
# TARGET_LINK_LIBRARIES( breed_sob_1 myatus )
