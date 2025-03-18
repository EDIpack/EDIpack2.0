.. _edipack2py:

EDIpack2 C-binding layer
====================================

`EDIpack2py` is made of a single Fortran module
:f:mod:`EDIPACK2PY`  providing Fortran-C 
interfaces of the  **EDIpack2** procedures leveraging  on the Fortran `iso_C_bindings` features. The generated
dynamic library can be used to setup library API of different type, as
for instance Python API in EDIpy2_ . 


Auxiliary procedures
###########################

This set of procedures interface several auxiliary functions of
`EDIpack2` library, including input reading and other setup functions.  

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2py/01_edipack2py_auxiliary




Bath procedures
###########################

This set of procedures interface the effective bath related functions of
`EDIpack2` library. This includes the replica/general bath setup up
functions and symmetry operations on the user side. 

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2py/02_edipack2py_bath




Input/Output procedures
##################################

This set of procedures interface the input/output functions of 
`EDIpack2`, including those returning dynamical functions.  

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2py/03_edipack2py_io




:math:`\chi^2` Fit procedures
########################################

This set of procedures interface the `EDIpack2` functions implementing
the :math:`\chi^2` fit of the supplied effective bath function.  

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2py/04_edipack2py_fit




Main procedures
#################################

This set of procedures interface the main `EDIpack2` module functions:
the impurity solver initialization, execution, finalization.   

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2py/05_edipack2py_main


   
.. _EDIpy2: https://github.com/edipack/EDIpy2.0
