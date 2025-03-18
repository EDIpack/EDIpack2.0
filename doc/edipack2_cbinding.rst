.. _edipack2_cbinding:

EDIpack2 C-binding layers
=================================================================

`EDIpack2ineq2py` is made of a single Fortran module
:f:mod:`EDIPACK2INEQ2PY`  providing Fortran-C 
interfaces of the  `EDIpack2` **and** `EDIpack2ineq` procedures. The
module leverages on the Fortran `iso_C_bindings` features. The generated
dynamic library can be used to setup library API of different type, as
for instance Python API in EDIpy2_ . 

This interface layer includes a number of additional procedures, on
top of those already available from :f:mod:`EDIpack2py`, to
cover the case of inequivalent impurities or sites (i.e. Real-Space-DMFT).




   
.. _EDIpy2: https://github.com/edipack/EDIpy2.0
