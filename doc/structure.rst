Library structure
#################################################################################

The  |edipack2|  library is made of 3 distinct parts. 
Each part builds into a library (static or dynamic) with its own
internal dependencies.  Alongside the library, specific module files
are generated and stored in the same location during installation.  

The organization of the library structure is as follows: 

**EDIpack2**: :ref:`edipack2`
*****************************************************************  
The core library, implementing the parallel Lanczos-based exact diagonalization of a quantum
impurity problem and dynamical correlation functions
calculations.
It builds into the static libary `libedipack2.a` and it is
accessible via the Fortran module :f:mod:`EDIPACK2`


**EDIpack2ineq**: :ref:`edipack2ineq`
******************************************************************
This is a sub-library implementing the inequivalent impurities
extension for |edipack2|, tackling the solution of systems represented by several
independent quantum impurity problems. Either serial or parallel execution
schemes with respect to the inequivalent impurities are provided, as
well as automatic deployment of the inequivalent files. 
This |edipack2ineq| builds into the same static library
`libedipack2.a`, yet it is accessible only through the module
:f:mod:`EDIPACK2INEQ`.
The naming convention of the |edipack2| procedures is preserved. 


**Fortran-C interface** :ref:`edipack2_cbinding`
********************************************************************************
A module implementing the Fortran-C interface for |edipack2|,
including |edipack2ineq| extension where included.
The  language interoperability leverages over the implicit Fortran
`Iso_C_Binding` features and is achieved using C-types to wrap the
relevant `EDIpack2` procedures.
The C-binding module serves as a key tool to setup additional API for 
|edipack2|, as for instance the Python API (see EDIpy2_ ).
This module builds into dynamic libraries
`libedipack2_cbinding.so/.dylib`. 
Dependency: `edipack2` and, if required, `edipack2ineq`.  



.. _EDIpy2: https://github.com/edipack/EDIpy2.0
