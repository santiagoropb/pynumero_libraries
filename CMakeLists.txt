
cmake_minimum_required(VERSION 3.2)

PROJECT( Pynumero )

SET(CMAKE_BUILD_TYPE release)

##################### Checks for compiler #####################
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)

IF(COMPILER_SUPPORTS_CXX11)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
ELSE()
  MESSAGE(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
  IF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 3.1)
      MESSAGE(FATAL_ERROR "CLANG version must be at least 3.1!")
    ENDIF()
  ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 4.7)
      MESSAGE(FATAL_ERROR "GCC version must be at least 4.7!")
    ENDIF()
  ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.0)
      MESSAGE(FATAL_ERROR "ICC version must be at least 12.0!")
    ENDIF()
  ELSEIF("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    IF(CMAKE_CXX_COMPILER_VERSION VERSION_LESS 12.0)
       MESSAGE(FATAL_ERROR "MSVC version must be at least 12.0!")
    ENDIF()
  ENDIF()
ENDIF()

# to find our dependencies
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

IF(WIN32)
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
    SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -static-libstdc++")
ENDIF()

OPTION(STATIC_LINK "STATIC_LINK" OFF)

# check fpr windows
IF(MSVC OR MSYS OR MINGW)
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
    IF(${STATIC_LINK})
        SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -static-libstdc++")
    ELSE()		
        SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
    endif()	
# check for apple
ELSEIF(APPLE)
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -O3 -fpermissive")
    IF(${STATIC_LINK})
	SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -static-libstdc++ -fpermissive")
    ELSE()
	SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -fpermissive")
    ENDIF()	
# check for linux
ELSE()
    SET(CMAKE_CXX_FLAGS_DEBUG "-g -O3")
    IF(${STATIC_LINK})
	SET(CMAKE_CXX_FLAGS_RELEASE "-O2 -static-libstdc++")
    ELSE()
	SET(CMAKE_CXX_FLAGS_RELEASE "-O2")
    ENDIF()
ENDIF()

OPTION(BUILD_ASL "BUILD_ASL" OFF)
OPTION(BUILD_HSL "BUILD_HSL" OFF)

IF(${BUILD_ASL})
    add_subdirectory(asl_interface)
    add_subdirectory(tests)
ENDIF()

IF(${BUILD_HSL})
    add_subdirectory(hsl_interface)
ENDIF()


