MODULE EDIPACK2_C
  !
  ! A Fortran-C interface module for EDIpack2 built around `ISO_C_BINDING`
  !
  USE COMMON
  implicit none


contains

  !include library functions
  include "edipack2/edipack2py_read_input.f90"
  include "edipack2/edipack2py_aux_funx.f90"
  !
  include "edipack2/edipack2py_bath.f90"
  !
  include "edipack2/edipack2py_bath_fit.f90"
  !
  include "edipack2/edipack2py_io.f90"
  include "edipack2/edipack2py_main.f90"
END MODULE EDIPACK2_C
