EDIpack2.0
################

A massively parallel Exact Diagonalization solver for Quantum Impurity problems.
***************************************************************************************************************

.. warning::
   The documentation is under construction


EDIpack2_ is a Lanczos based Exact Diagonalization method 
for the solution of generic Quantum Impurity problems,  exploiting MPI
distributed memory parallelization.

The 2.0 version extends the former EDIpack_ library enabling the solution of
single-site, multi-orbital models with different conserved
quantum numbers :math:`\vec{Q}` corresponding to separate operational
modes which, in `EDIpack2` software, are selected by the input
variable `ed_mode=normal,superc,nonsu2` as follow: 

* :math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` for which
  either the number of total or orbital spin up and down electrons is
  conserved: **NORMAL** 

* :math:`\vec{Q}=S_z`  with conserved total magnetization:  **SUPERConducting**  

* :math:`\vec{Q}=N_{\rm tot}`  where spin degrees freedom is not fully conserved:  **NON-SU(2)**


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


In the current state both the `normal` and `superc` operational modes
include **electron-phonon** coupling (local or Holstein phonons).

EDIpack2  is designed to obtain the lowest part of the
spectrum of the problem, thus it naturally works at **zero temperature**
but can also be used to explore **low temperature** properties.  
 
The EDIpack2 diagonalization algorithm is based on a massively
parallel execution of matrix-vector products, required in the context
of Lanczos-Arnoldi linear procedures.
However, substantial modifications have been introduced in this
version to address the *Superconducting* and *non-SU(2)* channels.  
An updated manuscript will be released soon. 

.. _EDIPACK2: **EDIPACK2**
.. _EDIPACK: https://github.com/edipack/EDIpack
.. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261


Authors
=================

The `EDIpack` libraries (1.0 and 2.0) have been developed as a
collective effort by different authors, each contributing to diverse
aspects of the library. The following list does not follow any
particular order:  

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
     Software requirements to install `EDIpack2`
     
:doc:`installation`
     Build, install and configure `EDIpack2`
     



     
Usage
=================

:doc:`quickstart`
     A quick start guide to  `EDIpack2` usage

:doc:`examples`
     Some examples illustrating the use of `EDIpack2` for simple test problems


Structure
======================
:doc:`structure`
     Overview of the `EDIpack2` library structure
     
     

EDIpack2
======================

:doc:`edipack2`
     A detailed overview  of the whole library with a thorough 
     description of the relevant modules, data types and procedures.


EDI2py
======================
:doc:`edipack2py`
     An overview of the Fortran-C interface for EDIpack2_


EDIpack2ineq
======================
:doc:`edipack2ineq`
     The inequivalent sites extension of  EDIpack2_ 
     

EDIpack2ineq2py
======================
:doc:`edipack2ineq2py`
     An overview of the Fortran-C interface for `EDIpack2ineq`
     



Browse Source Code
============================

:doc:`browsecode`
     Browse the `EDIpack2` structure



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
   :caption: EDIpack2py
   :maxdepth: 1
   :hidden:

   edipack2py


.. toctree::
   :caption: EDIpack2ineq
   :maxdepth: 2
   :hidden:

   edipack2ineq


.. toctree::
   :caption: EDIpack2ineq2py
   :maxdepth: 2
   :hidden:

   edipack2ineq2py
   

.. toctree::
   :caption: Browse code
   :maxdepth: 2
   :hidden:

   browsecode

.. toctree::
   :caption: External Links
   
   EDIpy2.0 on GitHub <https://github.com/edipack/EDIpy2.0>
   SciFortran on GitHub <https://github.com/SciFortran/SciFortran>




   

