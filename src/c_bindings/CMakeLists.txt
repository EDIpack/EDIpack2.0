# EDIpack2py library
FILE(GLOB EDI_C_SRC
  EDIPACK2_C_COMMON.f90
  EDIPACK2_C.f90  
)
ADD_LIBRARY(EDI_C_LIB OBJECT  ${EDI_C_SRC})

set(HEADER_OUTPUT_DIR ${CMAKE_BINARY_DIR}/include)
configure_file(${LIB_CBIND}/edipack2_cbinding.h
               ${HEADER_OUTPUT_DIR}/edipack2_cbinding.h)
