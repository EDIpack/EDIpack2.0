FILE(GLOB E2I_BATH_SRC
  E2I_BATH_REPLICA.f90
  E2I_BATH_USER.f90
  E2I_BATH.f90
  )

ADD_LIBRARY(E2I_BATH_LIB OBJECT  ${E2I_BATH_SRC})
