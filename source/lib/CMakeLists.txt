
set(TARGET myatus )

add_library( ${TARGET} SHARED
  UtilsRealWavefunction.hpp
  UtilsRealWavefunction.cpp
  GPUtilsRealWavefunction.hpp
  GPUtilsRealWavefunction.cpp
  UtilsComplexWavefunction.hpp
  UtilsComplexWavefunction.cpp
  GPUtilsComplexWavefunction.hpp
  GPUtilsComplexWavefunction.cpp
  MyParameterHandler.cpp 
  MyLogging.cpp 
  mmio.c 
)

DEAL_II_SETUP_TARGET( ${TARGET} )

TARGET_LINK_LIBRARIES( ${TARGET} m Boost::headers Boost::log Boost::log_setup)

