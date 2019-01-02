
enable_language(Fortran)
set(DEFAULT_COIN_HSL_PATH "${PROJECT_SOURCE_DIR}/third_party/HSL/coinhsl" )

# get the coinhsl library
if (USE_COIN_HSL_PATH)
    ADD_LIBRARY( coinhsl STATIC	"${USE_COIN_HSL_PATH}/ma27/ma27d.f")
    SET( COIN_HSL_INCLUDE_DIRS "${USE_COIN_HSL_PATH}/ma27")
    SET_PROPERTY(TARGET coinhsl PROPERTY POSITION_INDEPENDENT_CODE ON)
    SET(COIN_HSL_LIBRARIES coinhsl)
    SET( WITH_HSL true)
elseif(DEFAULT_COIN_HSL_PATH)
    ADD_LIBRARY( coinhsl STATIC "${DEFAULT_COIN_HSL_PATH}/ma27/ma27d.f")
    SET(COIN_HSL_INCLUDE_DIRS "${DEFAULT_COIN_HSL_PATH}/ma27")
    SET_PROPERTY(TARGET coinhsl PROPERTY POSITION_INDEPENDENT_CODE ON)
    SET(COIN_HSL_LIBRARIES coinhsl)
    SET( WITH_HSL true)
else()
    # for now throw an error - when we support more linear solvers, we can
    # allow this to proceed
    MESSAGE( FATAL_ERROR "This software currently requires MA27 provided by coinhsl. Please set USE_COIN_HSL_PATH to the location of coinhsl codes." )
    SET(COIN_HSL_LIRARIES "")
    MESSAGE( STATUS "*** HSL support not included - set USE_COIN_HSL_PATH to compile in HSL support." )
    SET( COIN_HSL_INCLUDE_DIRS "" )
    SET( WITH_HSL false)
endif()

if (${WITH_HSL})
    MESSAGE("-- COIN_HSL_INCLUDE_DIRS ${COIN_HSL_INCLUDE_DIRS}")
    SET(BLAS_LAPACK_INCLUDE_DIRS "")
    SET(BLAS_LAPACK_LIBRARIES "")
    ADD_LIBRARY(pynumero_HSL SHARED src/ma27Interface.cpp)
    SET_PROPERTY(TARGET pynumero_HSL PROPERTY POSITION_INDEPENDENT_CODE ON)
    TARGET_INCLUDE_DIRECTORIES( pynumero_HSL PUBLIC ${BLAS_LAPACK_INCLUDE_DIRS} PUBLIC ${COIN_HSL_INCLUDE_DIRS})
    TARGET_LINK_LIBRARIES( pynumero_HSL ${BLAS_LAPACK_LIBRARIES} ${COIN_HSL_LIBRARIES} ${CMAKE_DL_LIBS} gfortran)
    INSTALL(TARGETS pynumero_HSL		
            DESTINATION lib
           )
    INSTALL(DIRECTORY src
            DESTINATION include
	   )
endif()
