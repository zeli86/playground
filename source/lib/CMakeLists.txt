
set(TARGET myatus )

ADD_LIBRARY( ${TARGET} MyParameterHandler.cpp MyRealTools.cpp MyComplexTools.cpp MyLogging.cpp mmio.c )
TARGET_LINK_LIBRARIES( ${TARGET} m gomp Boost::headers Boost::log Boost::log_setup)

