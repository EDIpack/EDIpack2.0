# EDIpack2py library
FILE(GLOB EDI_C_SRC
  EDIPACK2_C_COMMON.f90
  EDIPACK2_C.f90  
)
ADD_LIBRARY(EDI_C_LIB OBJECT  ${EDI_C_SRC})

if(WITH_INEQ)
  set(WITH_INEQ_HEADER "int *Nineq")
else()
  set(WITH_INEQ_HEADER "void")
endif()

set(HEADER_OUTPUT_DIR ${CMAKE_BINARY_DIR}/include)
configure_file(${LIB_CBIND}/edipack2_cbinding.h.in
               ${HEADER_OUTPUT_DIR}/edipack2_cbinding.h)
