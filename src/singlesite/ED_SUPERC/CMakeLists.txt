# EDIpack library
FILE(GLOB ED_SUPERC_LIB_SRC
  ED_HAMILTONIAN_SUPERC_COMMON.f90
  ED_HAMILTONIAN_SUPERC_STORED_HxV.f90
  ED_HAMILTONIAN_SUPERC_DIRECT_HxV.f90
  ED_HAMILTONIAN_SUPERC.f90
  ED_DIAG_SUPERC.f90
  ED_GF_SUPERC.f90
  ED_OBSERVABLES_SUPERC.f90
  ED_RDM_SUPERC.f90
  )
ADD_LIBRARY(ED_SUPERC_LIB OBJECT  ${ED_SUPERC_LIB_SRC})



