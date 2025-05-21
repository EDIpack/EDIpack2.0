MODULE EDIPACK_C
  !:synopsis: Fortran-C bindings: main module
  !A Fortran-C interface module for EDIpack+EDIpack2ineq built around `ISO_C_BINDING`
  !
  USE EDIPACK_C_COMMON
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
  include "edipack/edipack_c_bindings_read_input.f90"
  include "edipack/edipack_c_bindings_aux_funx.f90"
  include "edipack/edipack_c_bindings_bath.f90"
  include "edipack/edipack_c_bindings_bath_fit.f90"
  include "edipack/edipack_c_bindings_parse_umatrix.f90"
#ifdef _WINEQ  
  include "edipack2ineq/edipack2ineq_c_bindings_aux_funx.f90"
  include "edipack2ineq/edipack2ineq_c_bindings_bath.f90"
  include "edipack2ineq/edipack2ineq_c_bindings_bath_fit.f90"
#endif
  !
  !But these two which are mutually exclusive
#ifdef _WINEQ  
  include "edipack2ineq/edipack2ineq_c_bindings_io.f90"
  include "edipack2ineq/edipack2ineq_c_bindings_main.f90"
#else
  include "edipack/edipack_c_bindings_io.f90"
  include "edipack/edipack_c_bindings_main.f90"
#endif
  !
END MODULE EDIPACK_C
