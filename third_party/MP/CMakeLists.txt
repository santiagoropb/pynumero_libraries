include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)
ExternalProject_Add(mp
		    PREFIX ASL_MP
		    GIT_TAG "3eb164b926c6d3e1e4cd600edb337039a7c824e0"
        	    GIT_REPOSITORY https://github.com/ampl/mp
        	    CMAKE_ARGS -DBUILD=asl -DCMAKE_INSTALL_PREFIX=${CMAKE_BINARY_DIR}/third_party/
        	    GIT_SHALLOW 1
		    GIT_PROGRESS 1)
