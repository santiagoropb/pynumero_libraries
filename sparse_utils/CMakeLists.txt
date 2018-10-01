ADD_LIBRARY( pynumero_SPARSE SHARED src/SparseUtils.cpp)
INSTALL(TARGETS pynumero_SPARSE
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
       )
