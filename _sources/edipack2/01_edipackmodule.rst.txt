.. _edipackmodule:

EDIpack2 Fortran Module
==========================


:f:mod:`EDIPACK2` is the top module of the |EDIpack2| library. It provides access to the all the relevant procedures of the library, realizing the Fortran API. The user needs to invoke use of this module to get access to |edipack2| as:

   .. code-block:: fortran

      program test
          USE EDIPACK2
	  ...

   		   
The module also contains a subset of the global and input variables that can be accessed in the userspace. 

.. preferred-crossrefs::
   :ed_bath/set_hreplica: f/ed_bath_replica/set_hreplica
   :ed_bath/set_hgeneral: f/ed_bath_replica/set_hgeneral
   :ed_bath/spin_symmetrize_bath: f/ed_bath_user/spin_symmetrize_bath
   :ed_bath/orb_symmetrize_bath: f/ed_bath_user/orb_symmetrize_bath
   :ed_bath/orb_equality_bath: f/ed_bath_user/orb_equality_bath
   :ed_bath/ph_symmetrize_bath: f/ed_bath_user/ph_symmetrize_bath
   :ed_bath/ph_trans_bath: f/ed_bath_user/ph_trans_bath
   :ed_bath/break_symmetry_bath: f/ed_bath_user/break_symmetry_bath
   :ed_bath/enforce_normal_bath: f/ed_bath_user/enforce_normal_bath
   :ed_bath/save_array_as_bath: f/ed_bath_user/save_array_as_bath


.. f:automodule::   edipack2
