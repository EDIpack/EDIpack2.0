MODULE EDIPACK2INEQ2PY
  USE COMMON
  USE EDIPACK2INEQ
  implicit none


contains

  !include library functions
  include "edipack2/edipack2py_read_input.f90"
  include "edipack2/edipack2py_aux_funx.f90"
  include "edipack2ineq/edipack2ineq2py_aux_funx.f90"
  !
  include "edipack2/edipack2py_bath.f90"
  include "edipack2ineq/edipack2ineq2py_bath.f90"
  !
  include "edipack2/edipack2py_bath_fit.f90"
  include "edipack2ineq/edipack2ineq2py_bath_fit.f90"
  !
  include "edipack2ineq/edipack2ineq2py_io.f90"
  include "edipack2ineq/edipack2ineq2py_main.f90"
END MODULE EDIPACK2INEQ2PY
