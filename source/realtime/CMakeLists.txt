
ADD_EXECUTABLE( rt_prop_mpi rt_prop_mpi.cpp )
DEAL_II_SETUP_TARGET( rt_prop_mpi )
TARGET_LINK_LIBRARIES( rt_prop_mpi myatus Boost::headers Boost::program_options )

ADD_EXECUTABLE( rt_prop_mpi_2 rt_prop_mpi_2.cpp )
DEAL_II_SETUP_TARGET( rt_prop_mpi_2 )
TARGET_LINK_LIBRARIES( rt_prop_mpi_2 myatus Boost::headers Boost::program_options )

ADD_EXECUTABLE( rt_prop_mpi_cs rt_prop_mpi_cs.cpp )
DEAL_II_SETUP_TARGET( rt_prop_mpi_cs )
TARGET_LINK_LIBRARIES( rt_prop_mpi_cs myatus Boost::headers Boost::program_options )

ADD_EXECUTABLE( rt_prop rtprop_1.cpp )
DEAL_II_SETUP_TARGET( rt_prop )
TARGET_LINK_LIBRARIES( rt_prop myatus Boost::headers Boost::program_options )

ADD_EXECUTABLE( rt_prop_2 rtprop_2.cpp )
DEAL_II_SETUP_TARGET( rt_prop_2 )
TARGET_LINK_LIBRARIES( rt_prop_2 myatus Boost::headers Boost::program_options )
