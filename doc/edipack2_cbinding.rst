.. _edipack2_cbinding:

EDIpack2 Fortran-C interface
=================================================================

The C-binding layer is made of a single Fortran module
:f:mod:`EDIPACK2_C`  which provides a Fortran-C 
interface. This is achieved using the interoperability features
implemented in :f:mod:`ISO_C_BINDINGS` Fortran module.

The :f:mod:`EDIPACK2_C` module essentially consists of a series of
C-type procedures which wraps one-by-one the relevant |edipack2|
functions as required to perform a full calculation.  

The build process generates a dynamic library which can then be used
to setup specific API for different languages. The Python interface is
implemented in EDIpy2_ . 

Because the interfaces included in this module ultimately call
|edipack2| functions using dedicated C-types, we can avoid entering in
the details of the implementation.
The interested reader is addressed to the source files available here
EDIPACK2C_ .



.. _EDIPACK2C: https://github.com/EDIpack/EDIpack2.0/tree/detach_rdmft/c_bindings   
.. _EDIpy2: https://github.com/edipack/EDIpy2.0
