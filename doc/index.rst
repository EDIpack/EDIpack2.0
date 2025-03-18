EDIpack2
################

A massively parallel Exact Diagonalization solver for Quantum Impurity problems.
***************************************************************************************************************

.. warning::
   The documentation is under construction


|edipack2| is a Lanczos based Exact Diagonalization method 
for the solution of generic Quantum Impurity problems,  exploiting MPI
distributed memory parallelization.

This version 2 extends the former EDIpack_ library by enabling the solution of
multi-orbital quantum impurity models with different conserved
quantum numbers :math:`\vec{Q}`:

* :math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` for which
  either the number of total or orbital spin up and down electrons is
  conserved: **NORMAL** 

* :math:`\vec{Q}=S_z`  with conserved total magnetization:  **SUPERConducting**  

* :math:`\vec{Q}=N_{\rm tot}`  where spin degrees freedom is not fully conserved:  **NON-SU(2)**

in |edipack2| these can be selected using the input
variable `ed_mode=normal,superc,nonsu2`

.. note::
   The `superc` mode deals with local *s*-wave pairing although in 
   diagonal and off-diagonal orbital channels. The actual
   implementation does not support long-range magnetic ordering.
   
.. note::
   The `nonsu2` operational mode deals with any situation in which
   spin symmetry group is not fully conserved, for instance in
   presence of local Spin-Orbit Coupling :math:`\vec{L} \cdot \vec{S}`,
   in-plane magnetization :math:`\langle S_x\rangle\gt0`  or in-plane
   triplet excitonic condensation, see `PhysRevB.107.115117`_. 

.. _PhysRevB.107.115117: https://journals.aps.org/prb/abstract/10.1103/PhysRevB.107.115117


In the actual development stage both the `normal` and `superc` modes
include **electron-phonon** coupling to Holstein phonons. 

|edipack2| is designed to obtain the lowest part of the
spectrum of the quantum impurity problem, thus it naturally works at **zero temperature**
yet it also supports determination of **low temperature** properties.  
 
The |edipack2| diagonalization algorithm is based on a massively
parallel execution of matrix-vector products, required in the context
of Lanczos-Arnoldi linear procedures.
However, substantial modifications have been introduced in this
version to address the *Superconducting* and *non-SU(2)* channels.  

.. _EDIPACK: https://github.com/edipack/EDIpack
.. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261


Authors
=================

The `EDIpack` libraries have been developed as a
collective effort by different authors, each contributing to diverse
aspects of the library.

* `Adriano Amaricci`_ (leading author)
  
* `Lorenzo Crippa`_
  
* `Samuele Giuli`_

* `Gabriele Bellomia`_

* `Giacomo Mazza`_

* Alberto Scazzola

* Luca de Medici
  
* Massimo Capone

.. _Adriano Amaricci: https://github.com/aamaricci
.. _Lorenzo Crippa: https://github.com/lcrippa    
.. _Samuele Giuli: https://github.com/SamueleGiuli
.. _Gabriele Bellomia: https://github.com/beddalumia
.. _Giacomo Mazza: https://github.com/GiacMazza


Installation
=================

:doc:`dependencies`
     Software requirements to install the |edipack2| library.
     
:doc:`installation`
     Build, install and configure the library in the OS.
     


Usage
=================

:doc:`quickstart`
     A quick start guide with two simple examples.

:doc:`examples`
     Further examples showcasing some potentialities of the software. 


Structure
======================
:doc:`structure`
     Global view of the |edipack2| library structure.
     
     

EDIpack2
======================
:doc:`edipack2`
     A detailed overview  of the whole library with a thorough 
     description of the relevant modules, data types and procedures.

     
EDIpack2ineq
======================
:doc:`edipack2ineq`
     The inequivalent impurities extension of |edipack2|


   
C-bindings
======================
:doc:`edipack2_cbinding`
     The Fortran-C interface for |edipack2| and |edipack2ineq|
     


Browse Source Code
============================

:doc:`browsecode`
     Browse the software source



.. Hidden TOCs

.. toctree::
   :caption: Installation
   :maxdepth: 2
   :hidden:

   dependencies
   installation

.. toctree::
   :caption: Usage
   :maxdepth: 2
   :hidden:

   quickstart
   examples
   
.. toctree::
   :caption: Structure
   :maxdepth: 1
   :hidden:

   structure

.. toctree::
   :caption: EDIpack2
   :maxdepth: 2
   :hidden:
      
   edipack2


.. toctree::
   :caption: EDIpack2ineq
   :maxdepth: 2
   :hidden:

   edipack2ineq

   
.. toctree::
   :caption: C-Bindings
   :maxdepth: 2
   :hidden:

   edipack2_cbinding
   

.. toctree::
   :caption: Browse code
   :maxdepth: 2
   :hidden:

   browsecode

.. toctree::
   :caption: External Links
   
   EDIpy2.0 on GitHub <https://github.com/edipack/EDIpy2.0>
   SciFortran on GitHub <https://github.com/SciFortran/SciFortran>




   

