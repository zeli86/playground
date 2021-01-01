
set(TARGET myatus )

ADD_LIBRARY( ${TARGET} MyParameterHandler.cpp MyTable.cpp mmio.c )
TARGET_LINK_LIBRARIES( ${TARGET} m gomp Boost::headers)

