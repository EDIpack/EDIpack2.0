MODULE EDIPACK2PY
  !
  ! A Fortran module including C-binding interfaces to the EDIpack2 procedures
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
END MODULE EDIPACK2PY
