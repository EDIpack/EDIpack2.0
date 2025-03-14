EDIpack2.0 Fortran Library
===========================

Here we give an overview of the structure of the EDIpack2.0 library,
with a detailed description of the relevant modules and procedures.


Library Frontend
###########################

The :f:mod:`EDIPACK2` module represents the main user interface (or Fortran API) 
of the **EDIpack2** library. This modules gives access to all the
available procedures and variables as needed to solve quantum impurity problems.  

.. toctree::
   :maxdepth: 1
   :glob:

   structure/edipackmodule
      

Core Solver Routines
###########################

The module :f:mod:`ED_MAIN` wraps the **EDIpack2** algorithm into
three functions, one for each step: initialization,
solution and finalization. 

.. toctree::
   :maxdepth: 1
   :glob:

   structure/main


General Environment
###########################

This part of the library includes a set of modules which contain
variables or procedures serving the rest of the code. 
The set includes:

:f:mod:`ED_INPUT_VARS` which defines global input variables and input
read method.

:f:mod:`ED_VARS_GLOBAL` containing global shared
variables and classes, such as the key class  :f:type:`effective_bath`
which gathers the different bath components according to value of :f:var:`bath_type` and :f:var:`ed_mode`.

:f:mod:`ED_SETUP`  which handles global memory (de)allocation and the evaluation
of the dimensions of the symmetry sectors.

Finally, :f:mod:`ED_AUX_FUNX`  defines procedures which are used throughout the code.

..
   :f:mod:`ed_vars_global`  contains the definition of simple data
   structures, such as  :f:type:`gfmatrix`  storing all weights and poles of the Green's functions    

.. toctree::
   :maxdepth: 2
   :glob:

   structure/general/*


Sparse Matrix
###########################

The module  :f:mod:`ED_SPARSE_MATRIX` contains the implementation of a
suitable Compact Row Stored (CRS) sparse matrix, used to store the sparse Hamiltonians in each symmetry
sector. 


.. toctree::
   :maxdepth: 2

   structure/classes/01_ed_sparse_matrix


   

EigenSpace
###########################

The module :f:mod:`ED_EIGENSPACE` implements an ordered
single linked list to store the lower part of the energy
spectrum.

.. toctree::
   :maxdepth: 2

   structure/classes/02_ed_eigenspace


GFmatrix
###########################

The module :f:mod:`ED_GFMATRIX` contains definition of a specific
class to store the results of the dynamical Lanczos calculation of
dynamical correlation functions in terms of weights and poles of the
Kallen-Lehmann representation.  


.. toctree::
   :maxdepth: 2

   structure/classes/03_ed_gfmatrix


Sectors
###########################

The :f:mod:`ED_SECTOR` module implements the construction of the symmetry
sectors for the three sets of quantum numbers :math:`\vec{Q}` considered in **EDIpack2**, that we
recall are:

* :math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` for which
  either the number of total or orbital spin up and down electrons is
  conserved: **NORMAL** 

* :math:`\vec{Q}=S_z`  with conserved total magnetization:  **SUPERConducting**  

* :math:`\vec{Q}=N_{\rm tot}`  where spin degrees freedom is not fully
  conserved:  **NON-SU(2)**
  

.. toctree::
   :maxdepth: 1
   :glob:

   structure/ed_sector



Bath 
###########################

In **EDIpack2** the bath is handled using a Reverse Communication
Strategy. All the procedures designed to define or handle the
discretized bath as well as those to evaluate suitable functions of
the bath are grouped in set of modules.   


      
.. toctree::
   :maxdepth: 1
   :glob:

   structure/bath
   

Hamiltonian
###########################

This part of the `EDIpack2.0` code implements the setup of the
sector Hamiltonian in each  operational modes,  corresponding to the
choice of one of the  symmetries implemented in the code selected by the variable :f:var:`ed_mode` =  :code:`normal, superc, nosu2`. See
:f:mod:`ed_sector` for more info about the symmetries implemented in
the code.

Any of three different mode is implemented in a distinct class of
modules performing all the  operations required for the construction
of the Hamiltonian and the definition of the corresponding matrix
vector product, required by the Arpack/Lanczos. 

.. toctree::
   :maxdepth: 1
   :glob:

   structure/index_hamiltonian


   

Exact Diagonalization
###########################


This part of the **EDIpack2** code implements the exact
diagonalization of the general, single-site, multi-orbital quantum
impurity problem in each of the specific symmetry implemented in the
code. The operational modes are selected by the variable
:f:var:`ed_mode` =  :code:`normal, superc, nosu2` (see
:f:mod:`ED_SECTOR` for more info about the symmetries implemented
in the code).

As above, we implemented the three different channels in
distinct class of modules. 

.. toctree::
   :maxdepth: 1
   :glob:

   structure/index_diag
   



Green's Functions 
###########################
This part of the **EDIpack2** code implements the calculation of the
impurity interacting Green's functions, self-energy functions and
impurity susceptibilities. Calculations are performed in  each operational mode,  corresponding to the
choice of the specific symmetry implemented in the code, i.e. which
quantum numbers are to be conserved. The operational modes are selected
by the variable :f:var:`ed_mode` =  :code:`normal, superc, nosu2` (see
:f:mod:`ED_SECTOR` for more info about the symmetries implemented in
the code).


.. toctree::
   :maxdepth: 1
   :glob:

   structure/index_greensfunctions


Observables
###########################

This part of the `EDIpack2.0` code implements the calculation of the
impurity observables and static correlations, such as density,
internal energy or double occupation. Calculations are performed in
each operational mode,  corresponding to the choice of the specific
symmetry implemented in the code, i.e. which quantum numbers are to be
conserved. The operational modes are selected by the variable
:f:var:`ed_mode` =  :code:`normal, superc, nosu2`. The observables
are printed in plain text files and are accessible to the user
through the routines lited in :f:mod:`ed_io`  (see
:f:mod:`ED_SECTOR` for more info about the symmetries implemented in
the code).


.. toctree::
   :maxdepth: 1
   :glob:

   structure/index_observables


   
Reduced Density Matrix
####################################

This part of the `EDIpack2.0` code implements the calculation of the
impurity Reduced Density Matrix (RDM). Calculations are performed in
each operational mode,  corresponding to the choice of the specific
symmetry implemented in the code, i.e. which quantum numbers are to be
conserved. The operational modes are selected by the variable
:f:var:`ed_mode` =  :code:`normal, superc, nosu2`.
The evaluation of the RDM exploit an efficient sparse algorithm to
avoid tracing over the exponentially many bath configuration.
The RDM is saved to file and made available to the user 
through the routines listed in :f:mod:`ED_IO`. 


.. toctree::
   :maxdepth: 1
   :glob:

   structure/index_rdm


   

Input/Output
###########################

This module provides access to the results of the exact
diagonalization. All quantities such as dynamical response functions,
self-energy components  or impurity observables  can be retrieved  
using specific functions. Additionally we provide procedure to perform
on-the-fly re-calculation of the impurity Green's functions and
self-energy on a given arbitrary set of points in the complex
frequency domain.    

.. toctree::
   :maxdepth: 1

   structure/io/ed_io




:math:`\chi^2` Fit
###########################

In  this module we provide to the user a generic function
:f:func:`ed_chi2_fitgf` performing the minimization of a user provided
Weiss field against the corresponding model of non-interacting
Anderson Green's function with the aim of updating the user bath parameters.


.. toctree::
   :maxdepth: 1

   structure/fit




