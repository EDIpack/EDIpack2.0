.. _edipack2_cbinding:

EDIpack2 Fortran-C interface
=================================================================

EDIpack2 offers Fortran-C interoperability making use of the
:f:mod:`ISO_C_BINDINGS` Fortran module.

The :f:mod:`EDIPACK2_C` module consists of a series of
C-type procedures which wrap one-by-one the relevant
functions as required to perform a full calculation.  

At compile time, if :code:`make` or :code:`make edipack2_cbinding` is issued,
a dynamic library called :code:`libedipack2_cbinding.so` will be created and
put in the same folder as the static Fortrna library.
An API to another language will need to interface to this library. 

A Python API called EDIpy2.0_ is provided and mantained.
An experimental Julia API called EDijl2.0_ is in development.

A :code:`c++` header file called :code:`edipack2_cbinding.h` is installed into the
include directory at build time. The list of variables and functions therein provided
is documented in the following. 
Caution needs to be applied for functions requiring array input parameters: since arrays 
have to be passed as raw pointers, additional
arrays containing the dimensions of the former need to be passed.

An example :code:`c++` is provided  in the examples_ folder of the GitHub repository.
A more nuanced usecase is available here_.

.. _EDIpy2.0: https://github.com/edipack/EDIpy2.0
.. _EDIjl2.0: https://github.com/edipack/EDIjl2.0
.. _examples: https://github.com/EDIpack/EDIpack2.0/tree/master/examples/cpp
.. _here: https://github.com/lcrippa/prematurata_la_dmft


.. doxygenfile:: edipack2_cbinding.h
   :project: edipack2
