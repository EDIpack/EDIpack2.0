A thin Fortran-Python interface layer
====================================

`EDI2py` is a single Fortran module interfacing the `EDIpack2.0`
library to Python language. 

The goal of this module is to wrap all the library procedures
leveraging  on the `iso-C-bindings` features of the modern Fortran
language. The module is automatically compiled into a shared
library during the installation step and then hosted  alongside the official EDIpack2.0 static
library 

The `EDI2py` shared library is loaded into the Python API of the
library `EDIpy2` (see relative documentation here). 






