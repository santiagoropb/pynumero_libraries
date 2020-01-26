
cmake_minimum_required(VERSION 3.2)

add_library( libpynumero_ASL SHARED "${CMAKE_CURRENT_SOURCE_DIR}/../src/AmplInterface.cpp")
add_dependencies(libpynumero_ASL mp)

IF(MSVC)
set_target_properties(libpynumero_ASL PROPERTIES LINK_FLAGS_RELEASE "/NODEFAULTLIB")
set_target_properties(libpynumero_ASL PROPERTIES LINK_FLAGS_DEBUG "/NODEFAULTLIB")
set_target_properties(libpynumero_ASL PROPERTIES STATIC_LIBRARY_FLAGS "/NODEFAULTLIB")	
ENDIF()

TARGET_INCLUDE_DIRECTORIES( libpynumero_ASL PUBLIC "${CMAKE_BINARY_DIR}/third_party/include/asl" )
TARGET_LINK_DIRECTORIES(libpynumero_ASL PUBLIC "${CMAKE_BINARY_DIR}/third_party/lib")
TARGET_LINK_LIBRARIES( libpynumero_ASL "asl" ${CMAKE_DL_LIBS})
INSTALL(TARGETS libpynumero_ASL
        DESTINATION lib
        )
INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../src/"
        DESTINATION include
	)
