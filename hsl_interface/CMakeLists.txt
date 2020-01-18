enable_language(Fortran)
set(DEFAULT_COIN_HSL_PATH "${PROJECT_SOURCE_DIR}/third_party/HSL/coinhsl" )

if (USE_COIN_HSL_PATH)
  SET(HSL_FILES ${USE_COIN_HSL_PATH})
elseif(DEFAULT_COIN_HSL_PATH)
  SET(HSL_FILES ${DEFAULT_COIN_HSL_PATH})
endif()

if (${BUILD_MA27})

  ADD_LIBRARY(pynumero_MA27 SHARED src/ma27Interface.cpp)

  if (IPOPT_PATH)
    SET( MA27_INCLUDE_DIRS "")
    FIND_LIBRARY(HSL_IPOPT_LIB
                 NAMES coinhsl libcoinhsl
                 HINTS "${IPOPT_PATH}/lib/"
                 )
    SET(MA27_LIBRARY "")
  else()
    SET(HSL_IPOPT_LIB "")
    SET( MA27_INCLUDE_DIRS "${HSL_FILES}/ma27")
    ADD_LIBRARY( ma27 STATIC "${HSL_FILES}/ma27/ma27d.f")
    SET_PROPERTY(TARGET ma27 PROPERTY POSITION_INDEPENDENT_CODE ON)
    SET(MA27_LIBRARY ma27)
    SET_TARGET_PROPERTIES(${MA27_LIBRARY} PROPERTIES COMPILE_FLAGS "-fno-common")
  endif()

  SET_PROPERTY(TARGET pynumero_MA27 PROPERTY POSITION_INDEPENDENT_CODE ON)
  TARGET_INCLUDE_DIRECTORIES( pynumero_MA27 PUBLIC ${MA27_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES( pynumero_MA27
                         ${HSL_IPOPT_LIB}
                         ${MA27_LIBRARY}
                         ${CMAKE_DL_LIBS}
                         gfortran)

  INSTALL(TARGETS pynumero_MA27
          DESTINATION lib
         )
  INSTALL(DIRECTORY src
          DESTINATION include
  	     )
endif()

if (${BUILD_MA57})
  ADD_LIBRARY(pynumero_MA57 SHARED src/ma57Interface.cpp)

  if (IPOPT_PATH)
    SET( MA57_INCLUDE_DIRS "")
    FIND_LIBRARY(HSL_IPOPT_LIB
                 NAMES coinhsl libcoinhsl
                 HINTS "${IPOPT_PATH}/lib/"
                 )
  else()
    MESSAGE( STATUS "Need IPOPT_PATH to build ma57 library")
  endif()

  SET_PROPERTY(TARGET pynumero_MA57 PROPERTY POSITION_INDEPENDENT_CODE ON)
  TARGET_INCLUDE_DIRECTORIES( pynumero_MA57 PUBLIC ${MA57_INCLUDE_DIRS})
  TARGET_LINK_LIBRARIES( pynumero_MA57
                         ${HSL_IPOPT_LIB}
                         ${CMAKE_DL_LIBS}
                         gfortran)

  INSTALL(TARGETS pynumero_MA57
          DESTINATION lib
         )

  INSTALL(DIRECTORY src
          DESTINATION include
 	       )
endif()
