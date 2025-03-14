Impurity Model Bath
########################

The construction and the handling of the bath is a crucial part of the
description of the generic single impurity Anderson problem.
The bath is described by two set of parameters: the local
hamiltonian :math:`\hat{h}^p` and the hybridization :math:`\hat{V}^p`,
for :math:`p=1,\dots,N_{bath}`. The first describes the local
properties of each bath element (be that a single electronic level or
a more complex structure made of few levels), the second describes the
coupling with the impurity.


The user accesses the bath as a double precision rank-1 array
containing in a given order all the parameters. This array is passed
as input to the library procedures and dumped into an internal data
structure :f:var:`effective_bath`, implemented in :f:mod:`ED_BATH_DMFT`.

We implemented different bath topologies, which can be
selected using the variable :f:var:`bath_type` = 
:code:`normal, hybrid, replica, general` .


For :f:var:`bath_type` = :code:`normal` (:code:`hybrid`) a number :f:var:`nbath` of electronic
levels are coupled to each orbital level (to any orbital level) of the
impurity site. The bath local Hamiltonian is diagonal  
:math:`\hat{h}^p\equiv\epsilon^p_a\delta_{ab}` while the
hybridizations are:  :math:`\hat{V}^p=V^p_{a}\delta_{ab}`
(:math:`\hat{V}^p=V^p_{ab}`).  
If  :f:var:`ed_mode` = :code:`superc` the bath includes a set of parameters
:math:`\Delta_p` describing the superconductive amplitude on each bath
level. 


For :f:var:`bath_type` = :code:`replica` (:code:`general`) a number :f:var:`nbath` of copies of the
impurity structure are  coupled to the impurity itself.  Each bath
element is made of a number :math:`N_{orb}` of electronic levels,
i.e. the number of orbitals in the impurity site.
The hybridization to the impurity site is
:math:`\hat{V}^p=V^p_{a}\delta_{ab}` (:math:`\hat{V}^p=V^p_{ab}`).  
The local bath  Hamiltonian is :math:`\hat{h}^p = \sum_{m=1}^{M} \lambda^p_m O_m`.
The set  :math:`\{O\}_m` is a user defined matrix basis for the impurity
Hamiltonian or, equally, for the local Hamiltonian of the lattice
problem. The numbers  :math:`\lambda^p_m\in{\mathbb R}` are
variational parameters.  

.. note::
   The enumeration of the total bath electronic levels is different
   among the different cases. This number is automatically evaluated
   upon calling :f:func:`get_bath_dimension`, see
   :f:mod:`ED_BATH_DIM`.


.. note::
   The :code:`replica, general` bath topologies are available also for
   the superconductive case :f:var:`ed_mode` = :code:`superc`. In this case the
   structure of the matrix basis should be set to the proper
   multi-orbital Nambu basis, so that off-diagonal blocks corresponds
   to anomalous components.


   
.. f:automodule::   ed_bath
   :hide-output: True

Bath Auxiliary
++++++++++++++++

In this set of modules we implement a number of auxiliary procedures
which are required to enumerate the bath levels, performs all the
relevant checks on the user input bath, build up of the
:code:`replica/general` matrix basis or apply given symmetry operation
on the user side. 

.. toctree::
   :maxdepth: 1

   bath_aux/00_ed_vars_global
   bath_aux/01_ed_bath_aux
   bath_aux/02_ed_bath_dim
   bath_aux/03_ed_bath_user
   bath_aux/04_ed_bath_replica
	      
Bath DMFT
+++++++++++++

In :f:mod:`ED_BATH_DMFT` we implement operations on the  :f:type:`effective_bath` data
structure: a suitable representation of the effective bath used
internally in the code. We refer to the generic shared instance of this bath
as :f:var:`dmft_bath`.  Depending on the value :f:var:`bath_type` this quantity collect
different bath parameters which can be directly accessed in the
construction of symmetry sectors Hamiltonian.  

.. toctree::
   :maxdepth: 1

   bath_dmft/01_ed_bath_dmft


Bath Functions
++++++++++++++++

In :f:mod:`ED_BATH_FUNCTIONS` we implement on-the-fly construction of
the hybridization functions :math:`\Delta(z) = \sum_p
\hat{V}^p\left[z-\hat{h}^p \right]^{-1}\hat{V}^p`, as well as of the
non-interacting Anderson Green's functions
:math:`G_0(z) = \left[z +\mu - H_{loc} - \Delta(z) \right]^{-1}` for
all different cases selected by :f:var:`ed_mode` and :f:var:`bath_type`.

.. toctree::
   :maxdepth: 1

   bath_functions/ed_bath_functions

   
