.. _edipack2ineq2py:

EDIpack2ineq C-binding layer
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


Auxiliary procedures
###########################

Interface input reading and local Hamiltonian setup for both single
and multiple impurities.

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2ineq2py/01_edipack2ineq2py_auxiliary




Bath procedures
###########################

Includes procedures to set the :code:`replica/general` matrix basis
decomposition also in the inequivalent impurities case. Provides
inequivalent baths symmetry operations. 

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2ineq2py/02_edipack2ineq2py_bath




Input/Output procedures
##################################

Extend the input/output procedure to the case of inequivalent
impurities. 

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2ineq2py/03_edipack2ineq2py_io




:math:`\chi^2` Fit procedures
########################################

Extend the bath optimization function to the case of multiple
independent impurities.

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2ineq2py/04_edipack2ineq2py_fit




Main procedures
#################################

This set of procedures interface the main `EDIpack2` and
`EDIpack2ineq` functions: the impurities solver initialization, execution, finalization.   

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2ineq2py/05_edipack2ineq2py_main


   
.. _EDIpy2: https://github.com/edipack/EDIpy2.0
