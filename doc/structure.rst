Libary structure
#################################################################################

The  **EDIpack2**  library is constructed out of 4 distinct parts. The
modular structure ensures a better development scalability and offers
a better way to provide access to additional features. 

Each part builds into a static or dynamic library with its own
internal dependencies and it handles a specific task.   
All the libraries and module files generated at build time are stored in the same
location for ease of use. 

The organization of the library structure is as follows: 

**EDIpack2**: :ref:`edipack2`
*****************************************************************  
  The core library,
  implementing the parallel Lanczos-based exact diagonalization of a quantum
  impurity problem and dynamical correlation functions
  calculations.
  It builds into the static libary `libedipack2.a`.
 

**EDIpack2ineq**: :ref:`edipack2ineq`
******************************************************************
A Fortran sub-library
  implementing the inequivalent lattice extension for the EDIpack2_
  library. This part enables to solve systems represented by several
  independent quantum impurity problems, ideally one per inequivalent
  site, using different parallelization schemes and automatic
  deployment of the files while keeping the naming convention of
  `edipack2`.
  It Builds into the static library `libedipack2ineq.a`. Dependency: `edipack2`


**Fortran-C interface** :ref:`edipack2_cbinding`
********************************************************************************
  two modules implementing Fortran-C interface. The  interoperability is built by
  wrapping with C-types the relevant `EDIpack2` and, respectively,
  `EDIpack2ineq` procedures leveraging over the implicit Fortran
  `Iso_C_Binding` features. These two modules serve as a key tool to setup additional API for the
`EDIpack2` library, as for instance the Python API (see EDIpy2_ ). The
C-bindings modules build into dynamic libraries
`libedipack2_cbinding.so/.dylib` or
`libedipack2ineq_cbinding.so/.dylib` Dependency: `edipack2` and `edipack2ineq`.  



.. _EDIpy2: https://github.com/edipack/EDIpy2.0
