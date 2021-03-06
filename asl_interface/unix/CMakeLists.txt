cmake_minimum_required(VERSION 3.2)

set(MP_PATH ""
    CACHE FILEPATH
    "Set the path to the ampl-mp package to link against ampl-mp asl library")

# set the default directory for ASL
set(DEFAULT_ASL_PATH "${PROJECT_SOURCE_DIR}/third_party/ASL/solvers")

set(ASL_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/AmplInterface.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/../src/AssertUtils.hpp"
   )

if (IPOPT_PATH)
   add_library( pynumero_ASL SHARED ${ASL_SOURCES})
   find_library(ASL_LIB
	              NAMES coinasl libcoinasl
                HINTS "${IPOPT_PATH}/lib/"
		           )
   set(ASL_INCLUDE_DIRS "${IPOPT_PATH}/include/coin/ThirdParty/")
   set( ASL_LIBRARIES "${ASL_LIB}" )
   set( WITH_AMPL true)

elseif (MP_PATH)
   add_library( pynumero_ASL SHARED ${ASL_SOURCES})

   find_library(ASL_LIB
	              NAMES asl libdasl
                HINTS "${MP_PATH}/lib/"
		           )

   set(ASL_INCLUDE_DIRS "${MP_PATH}/include/asl")
   set( ASL_LIBRARIES "${ASL_LIB}" )
   set( WITH_AMPL true)

elseif(DEFAULT_ASL_PATH)
   ADD_LIBRARY( pynumero_ASL SHARED ${ASL_SOURCES})
   set(ASL_INCLUDE_DIRS "${DEFAULT_ASL_PATH}")
   set( ASL_LIBRARIES "${DEFAULT_ASL_PATH}/amplsolver.a" )
   set( WITH_AMPL true)
else()
   MESSAGE( STATUS "*** ASL (AMPL) support not included - set USE_ASL_PATH to compile in ASL support." )
   set( ASL_INCLUDE_DIRS "" )
   set( ASL_LIBRARIES "" )
   set( WITH_AMPL false)
endif()


if (${WITH_AMPL})
   MESSAGE("-- ASL_INCLUDE_DIRS ${ASL_INCLUDE_DIRS}")
   MESSAGE("-- ASL_LIBRARIES ${ASL_LIBRARIES}")
   TARGET_INCLUDE_DIRECTORIES( pynumero_ASL PUBLIC ${ASL_INCLUDE_DIRS} )
   TARGET_LINK_LIBRARIES( pynumero_ASL ${ASL_LIBRARIES} ${CMAKE_DL_LIBS})
   INSTALL(TARGETS pynumero_ASL
           DESTINATION lib
          )
   INSTALL(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/../src/"
           DESTINATION include
	   )
endif()
