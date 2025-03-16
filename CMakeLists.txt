##################################################
# DMFTED PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0)
PROJECT(edipack2 Fortran)
SET(EDI ${PROJECT_NAME})
SET(EDI2INEQ edipack2ineq)
SET(EDI2PY edipack2py)
SET(EDI2INEQ2PY edipack2ineq2py)
SET(VERSION 2.0.8)


MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)


# CMakeGraphVizOptions.cmake
set(GRAPHVIZ_CUSTOM_TARGETS TRUE)


# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Set some options the user may choose
OPTION(USE_MPI "Set MPI pre-compilation flag" ON)
OPTION(VERBOSE "Build a verbose Makefile" OFF)
OPTION(WITH_BLAS_LAPACK "Use internal Blas/Lapack" OFF)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)



############################################################
# Look for MPI if needed
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MpiConfig.cmake)



############################################################
# Look for required libraries: scifor
# (includes Blas/Lapack/Scalapack/Quadpack/Minpack/FFTpack/Arpack/...) 
############################################################
# Find Scifor
FIND_PACKAGE(PkgConfig)
PKG_SEARCH_MODULE(SCIFOR REQUIRED scifor)
IF(${SCIFOR_FOUND})
  MESSAGE(STATUS "SCIFOR found at: ${SCIFOR_PREFIX}")
ELSE()
  MESSAGE(FATAL_ERROR "
Can not find SciFortran in this system. Check SciFor installation. Try adding the provided PKG-CONFIG file to the environment:
$ export PKG_CONFIG_PATH=<path_to_scifor_dir>/etc:$PKG_CONFIG_PATH")
ENDIF()
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})


############################################################
# Configure Fortran Compiler Main options
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/FortranCompilerConfig.cmake)




############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/GitConfig.cmake)




############################################################
# Define the actual files and folders that make up the build
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MainConfig.cmake)



############################################################
############################################################
############################################################
# DEFINE HERE THE DETAILS OF THE LIBRARY CONFIGURATION 

# Define some directories
SET(LIB_SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB_ETC ${CMAKE_SOURCE_DIR}/etc)
SET(LIB_TEST ${CMAKE_SOURCE_DIR}/test)
SET(LIB_PY ${CMAKE_SOURCE_DIR}/python)
SET(LIB_INEQ ${CMAKE_SOURCE_DIR}/ineq)
SET(LIB_INEQPY ${CMAKE_SOURCE_DIR}/python/ineq)
SET(LIB_ENV ${CMAKE_SOURCE_DIR}/etc/environment_modules)
#
SET(LIB_TMP_INC ${CMAKE_BINARY_DIR}/include)
SET(LIB_TMP_ETC ${CMAKE_BINARY_DIR}/etc)
SET(LIB_TMP_VER ${CMAKE_BINARY_DIR}/version)
#
SET(LIB_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(LIB_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(LIB_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(LIB_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(LIB_TARGET_DIR ${CMAKE_INSTALL_PREFIX})

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})






####################################################
#       EDIpack2 library
####################################################
ADD_SUBDIRECTORY(${LIB_SRC}/ED_BATH)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_NORMAL)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_SUPERC)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_NONSU2)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_IO)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_FIT)

ADD_SUBDIRECTORY(${LIB_SRC})


ADD_DEPENDENCIES(ED_BATH_LIB
  ED_INIT_LIB
)


ADD_DEPENDENCIES(ED_NORMAL_LIB
  ED_INIT_LIB
  ED_BATH_LIB
)

ADD_DEPENDENCIES(ED_SUPERC_LIB
  ED_INIT_LIB
  ED_BATH_LIB
)


ADD_DEPENDENCIES(ED_NONSU2_LIB
  ED_INIT_LIB
  ED_BATH_LIB
)



ADD_DEPENDENCIES(ED_FIT_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  ED_SUPERC_LIB
  ED_NONSU2_LIB
)



ADD_DEPENDENCIES(ED_GLOB_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  ED_SUPERC_LIB
  ED_NONSU2_LIB
  ED_FIT_LIB
)

ADD_DEPENDENCIES(ED_IO_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  ED_SUPERC_LIB
  ED_NONSU2_LIB
  ED_GLOB_LIB
)


ADD_DEPENDENCIES(ED_LIB_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  ED_SUPERC_LIB
  ED_NONSU2_LIB
  ED_FIT_LIB
  ED_IO_LIB
)

ADD_LIBRARY(${EDI} STATIC
  $<TARGET_OBJECTS:ED_INIT_LIB>
  $<TARGET_OBJECTS:ED_BATH_LIB>
  $<TARGET_OBJECTS:ED_NORMAL_LIB>
  $<TARGET_OBJECTS:ED_SUPERC_LIB>
  $<TARGET_OBJECTS:ED_NONSU2_LIB>
  $<TARGET_OBJECTS:ED_FIT_LIB>
  $<TARGET_OBJECTS:ED_GLOB_LIB>
  $<TARGET_OBJECTS:ED_IO_LIB>
  $<TARGET_OBJECTS:ED_LIB_LIB>
)




####################################################
#   EDI2PY: EDIpack2 Python interface: 
####################################################
ADD_SUBDIRECTORY(${LIB_PY}/c_bindings EXCLUDE_FROM_ALL)
ADD_DEPENDENCIES(ED2PY_LIB ${EDI})
ADD_LIBRARY(${EDI2PY} SHARED EXCLUDE_FROM_ALL $<TARGET_OBJECTS:ED2PY_LIB>)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  TARGET_LINK_OPTIONS(${EDI2PY} PRIVATE -Wl,-z,noexecstack)
endif()
TARGET_LINK_LIBRARIES(${EDI2PY} PRIVATE ${EDI})
TARGET_LINK_LIBRARIES(${EDI2PY} PRIVATE ${SCIFOR_LINK_LIBRARIES})
SET(EDI2PY_NAME ${CMAKE_SHARED_LIBRARY_PREFX}${EDI2PY}${CMAKE_SHARED_LIBRARY_SUFFIX})




####################################################
#   EDIpack2ineq: inequivalent sites extension lib
####################################################
ADD_SUBDIRECTORY(${LIB_INEQ}/E2I_BATH EXCLUDE_FROM_ALL)  
ADD_SUBDIRECTORY(${LIB_INEQ}/E2I_IO   EXCLUDE_FROM_ALL)   
ADD_SUBDIRECTORY(${LIB_INEQ}/E2I_FIT  EXCLUDE_FROM_ALL)  
ADD_SUBDIRECTORY(${LIB_INEQ}          EXCLUDE_FROM_ALL)  
ADD_DEPENDENCIES(E2I_BATH_LIB E2I_INIT_LIB )  
ADD_DEPENDENCIES(E2I_FIT_LIB  E2I_INIT_LIB E2I_BATH_LIB)
ADD_DEPENDENCIES(E2I_IO_LIB   E2I_INIT_LIB E2I_BATH_LIB)
ADD_DEPENDENCIES(E2I_LIB      E2I_INIT_LIB E2I_BATH_LIB E2I_FIT_LIB E2I_IO_LIB)
ADD_LIBRARY(${EDI2INEQ} STATIC EXCLUDE_FROM_ALL 
  $<TARGET_OBJECTS:E2I_INIT_LIB>
  $<TARGET_OBJECTS:E2I_BATH_LIB>
  $<TARGET_OBJECTS:E2I_FIT_LIB>
  $<TARGET_OBJECTS:E2I_IO_LIB>
  $<TARGET_OBJECTS:E2I_LIB>
)
SET(${EDI2INEQ}_NAME ${CMAKE_SHARED_LIBRARY_PREFX}${EDI2INEQ}${CMAKE_SHARED_LIBRARY_SUFFIX})


####################################################
# EDI2INEQ2PY: EDIpack2 + EDIpack2ineq  Python interface
####################################################
ADD_SUBDIRECTORY(${LIB_INEQPY}/c_bindings EXCLUDE_FROM_ALL)
ADD_DEPENDENCIES(EDI2PY_LIB ${EDI2INEQ})
ADD_LIBRARY(${EDI2INEQ2PY} SHARED EXCLUDE_FROM_ALL  $<TARGET_OBJECTS:EDI2PY_LIB>)
ADD_DEPENDENCIES(${EDI2INEQ2PY} ${EDI2INEQ})
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  TARGET_LINK_OPTIONS(${EDI2INEQ2PY} PRIVATE -Wl,-z,noexecstack)
endif()
TARGET_LINK_LIBRARIES(${EDI2INEQ2PY} PRIVATE ${EDI})
TARGET_LINK_LIBRARIES(${EDI2INEQ2PY} PRIVATE ${EDI2INEQ})
TARGET_LINK_LIBRARIES(${EDI2INEQ2PY} PRIVATE ${SCIFOR_LINK_LIBRARIES})
SET(${EDI2INEQ2PY}_NAME ${CMAKE_SHARED_LIBRARY_PREFX}${EDI2INEQ2PY}${CMAKE_SHARED_LIBRARY_SUFFIX})









############################################################
############################################################
############################################################

SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)





#Final messages + additional custom targets:
# - distclean: cleans the generated cmake files
# - uninstall: provides uninstall of (almost) all the installed files
# - doc: creates documentation using sphinx (if found)
# - test: compiles and launches the test in ${LIB_TEST} directory 
INCLUDE(${CMAKE_MODULE_PATH}/PostConfig.cmake)

