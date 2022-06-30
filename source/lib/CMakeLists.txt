
set(TARGET myatus )

add_library( ${TARGET} 
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

TARGET_LINK_LIBRARIES( ${TARGET} m gomp Boost::headers Boost::log Boost::log_setup)

