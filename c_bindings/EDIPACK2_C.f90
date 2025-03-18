MODULE EDIPACK2_C
  !
  ! A Fortran-C interface module for EDIpack2+EDIpack2ineq built around `ISO_C_BINDING`
  !
  USE COMMON
  USE ISO_C_BINDING
#ifdef _WINEQ
  USE EDIPACK2INEQ
#endif
  implicit none
  !
#ifdef _WINEQ
  integer(c_int),bind(c, name="has_ineq") :: has_ineq=1 
#else
  integer(c_int),bind(c, name="has_ineq") :: has_ineq=0
#endif
  !
contains
  !
  !Most of the Fortran-C interface is additive:
  include "edipack2/edipack2py_read_input.f90"
  include "edipack2/edipack2py_aux_funx.f90"
  include "edipack2/edipack2py_bath.f90"
  include "edipack2/edipack2py_bath_fit.f90"
#ifdef _WINEQ  
  include "edipack2ineq/edipack2ineq2py_aux_funx.f90"
  include "edipack2ineq/edipack2ineq2py_bath.f90"
  include "edipack2ineq/edipack2ineq2py_bath_fit.f90"
#endif
  !
  !But these two which are mutually exclusive
#ifdef _WINEQ  
  include "edipack2ineq/edipack2ineq2py_io.f90"
  include "edipack2ineq/edipack2ineq2py_main.f90"
#else
  include "edipack2/edipack2py_io.f90"
  include "edipack2/edipack2py_main.f90"
#endif
  !
END MODULE EDIPACK2_C
