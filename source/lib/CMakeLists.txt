
ADD_LIBRARY( myatus MyParameterHandler.cpp my_table.cpp mmio.c pugixml.cpp )
TARGET_LINK_LIBRARIES( myatus m gomp )

ADD_EXECUTABLE( testread_xml testread_xml.cpp )
TARGET_LINK_LIBRARIES( testread_xml myatus m gomp )

