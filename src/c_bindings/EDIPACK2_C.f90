MODULE EDIPACK2_C
  !
  ! A Fortran-C interface module for EDIpack2+EDIpack2ineq built around `ISO_C_BINDING`
  !
  USE EDIPACK2_C_COMMON
  USE ISO_C_BINDING
#ifdef _WINEQ
  USE EDIPACK2INEQ
#endif
  implicit none
  !
  !
contains
  !
  !Most of the Fortran-C interface is additive:
  include "edipack2/edipack_c_binding_read_input.f90"
  include "edipack2/edipack_c_binding_aux_funx.f90"
  include "edipack2/edipack_c_binding_bath.f90"
  include "edipack2/edipack_c_binding_bath_fit.f90"
#ifdef _WINEQ  
  include "edipack2ineq/edipack2ineq_c_binding_aux_funx.f90"
  include "edipack2ineq/edipack2ineq_c_binding_bath.f90"
  include "edipack2ineq/edipack2ineq_c_binding_bath_fit.f90"
#endif
  !
  !But these two which are mutually exclusive
#ifdef _WINEQ  
  include "edipack2ineq/edipack2ineq_c_binding_io.f90"
  include "edipack2ineq/edipack2ineq_c_binding_main.f90"
#else
  include "edipack2/edipack_c_binding_io.f90"
  include "edipack2/edipack_c_binding_main.f90"
#endif
  !
END MODULE EDIPACK2_C
