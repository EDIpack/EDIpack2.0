.. _edipack2ineq:

EDIpack2 inequivalent sites extension
=======================================================

`EDIpack2ineq` is a reduced sub-library of `EDIpack2` which extends
its functionalities of to handle systems  with several 
independent impurity problems, i.e. originating from inequivalent
sites. These are usually emerging from RealSpace-DMFT treatment of
unit cells with different inequivalent atoms or large super-cells with
somehow broken translational symmetries (heterostructures, disordered
systems, etc.).

As for its scope the structure of this sub-library closely mimics that
of `EDIpack2` adding, where needed, new procedures dealing with the
inequivalent sites. A standard fortran interface is then used to group
all the procedures (from `EDIpack2` and `EDIpack2ineq`) under the same
name.      






Library Frontend
###########################

The :f:mod:`EDIPACK2INEQ` module represents the main user interface
(or Fortran API). This modules gives access to all the available
procedures and variables as needed to solve quantum impurity problems.


.. toctree::
   :maxdepth: 1
   :glob:

   edipack2ineq/edipack2ineqmodule
      

   
Core Solver Routines
###########################

The module :f:mod:`E2I_MAIN` extends wrapping of the main algorithms into
three functions, one for each step: initialization,
solution and finalization. 

.. toctree::
   :maxdepth: 1
   :glob:

   edipack2ineq/main


General Environment
###########################

This part of the library includes a set of global variables and
procedures  for `EDIpack2ineq`.  This includes:

:f:mod:`E2I_VARS_GLOBAL`: contains global shared
variables storing site dependent quantities in the memory

Finally, :f:mod:`E2I_AUX_FUNX`  defines procedures which are used throughout the code.

.. toctree::
   :maxdepth: 2
   :glob:

   edipack2ineq/general/*


   

Bath 
###########################

The module :f:mod:`E2I_BATH` contains some extensions to procedures
enabling user side control of the bath setup and the implementation of
conventional symmetry operations. 

      
.. toctree::
   :maxdepth: 2
   :glob:

   edipack2ineq/bath/*
   



   

Input/Output
###########################

The module :f:mod:`E2I_IO` provides to extended input and output
procedures returning observables or dynamical functions per
inequivalent site. 

.. toctree::
   :maxdepth: 1

   edipack2ineq/io




:math:`\chi^2` Fit
###########################

The module :f:mod:`E2I_BATH_FIT` contains the extensions of the
generic function :f:func:`ed_chi2_fitgf` to perform 
multiple independent minimizations for all inequivalent sites.  

.. toctree::
   :maxdepth: 1

   edipack2ineq/fit

 






