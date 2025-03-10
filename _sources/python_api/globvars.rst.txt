Global variables
=================

These are global variables of the :code:`edipy2` module. They form a subset of the global variables of the EDIpack2 library. 
Along with all the other global variables, they can be set in the input file, 
and are read when calling the :func:`edipy2.global_env.read_input` function.

The exposed global variables can be accessed as properties of the :class:`edipy2.global_env` class.

.. code-block:: python

    import numpy as np
    from edipy2 import global_env as ed
   
    ed.Nspin = 1            # set a global variable
    mylocalvar = ed.Nspin   # assing to a local variable (the value of mylocalvar will not change if ed.Nspin changes)
    print(ed.Nspin)         # all functions can have global variables as arguments
    np.arange(ed.Nspin)


.. data:: edipy2.global_env.beta

   Value of the inverse temperature, at T=0 is used as a IR cut-off
   
   Accesses :f:var:`beta`
   
   :type: float
   :default: 1000.0

.. data:: edipy2.global_env.Jh

   Value of the Hund's coupling
   
   Accesses :f:var:`jh`
   
   :type: float
   :default: 0.0
   
.. data:: edipy2.global_env.dmft_error

   Error threshold for DMFT convergence
   
   Accesses :f:var:`dmft_error`
   
   :type: float
   :default: 1e-05
   
.. data:: edipy2.global_env.ed_total_ud

   Flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw.
   
   Accesses :f:var:`ed_total_ud`
   
   :type: bool
   :default: True
   
.. data:: edipy2.global_env.ed_twin

   Flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry
   
   Accesses :f:var:`ed_twin`
   
   :type: bool
   :default: False
   
.. data:: edipy2.global_env.eps

   Broadening on the real-axis
   
   Accesses :f:var:`eps`
   
   :type: float
   :default: 1e-02

.. data:: edipy2.global_env.Jx

   Value of the spin exchange coupling
   
   Accesses :f:var:`jx`
   
   :type: float
   :default: 0.0

.. data:: edipy2.global_env.Jp

   Value of the pair hopping coupling
   
   Accesses :f:var:`jp`
   
   :type: float
   :default: 0.0
   
.. data:: edipy2.global_env.Lmats

   Number of frequencies, Matsubara axis
   
   Accesses :f:var:`lmats`
   
   :type: int
   :default: 4096
  
.. data:: edipy2.global_env.LOGfile

   Log unit
   
   Accesses :f:var:`logfile`
   
   :type: int
   :default: 6
   
.. data:: edipy2.global_env.Lpos

   Number of points for the lattice PDF
   
   Accesses :f:var:`lpos`
   
   :type: int
   :default: 100
  
.. data:: edipy2.global_env.Lreal

   Number of frequencies, real axis
   
   Accesses :f:var:`lreal`
   
   :type: int
   :default: 5000

.. data:: edipy2.global_env.Ltau

   Number of imaginary time points
   
   Accesses :f:var:`ltau`
   
   :type: int
   :default: 1024

.. data:: edipy2.global_env.Nbath

   Number of bath levels. See the specifics of the bath geometries
   
   Accesses :f:var:`nbath`
   
   :type: int
   :default: 6
   
.. data:: edipy2.global_env.Nloop

   Maximum number of DMFT loops
   
   Accesses :f:var:`nloop`
   
   :type: int
   :default: 100

.. data:: edipy2.global_env.Norb

   Number of correlated orbitals. Maximum 5 orbitals are supported
   
   Accesses :f:var:`norb`
   
   :type: int
   :default: 1

.. data:: edipy2.global_env.Nph

   Max number of phonons allowed (cut off)
   
   Accesses :f:var:`nph`
   
   :type: int
   :default: 0
   
.. data:: edipy2.global_env.nread

   Value of the target density for fixed density calculations. If valued 0, it is discarded.
   
   Accesses :f:var:`nread`
   
   :type: float
   :default: 0.0

.. data:: edipy2.global_env.Nspin

   Number of explicitly defined spin degrees of freedom. If Nspin=1, the two spin block 
   of the Hamiltonian, Green's function, self-energy and so on are assumed equal.
   If Nspin=2 they may differ (e.g. for non-SU(2) or magnetic systems).
   The superconductive variant of the code requires Nspin=1
   
   Accesses :f:var:`nspin`
   
   :type: int
   :default: 1
   
.. data:: edipy2.global_env.Nsuccess

   Number of successive iterations below threshold for convergence
   
   Accesses :f:var:`nsuccess`
   
   :type: int
   :default: 1
   
.. data:: edipy2.global_env.sb_field

   Value of a symmetry breaking field for magnetic solutions
   
   Accesses :f:var:`sb_field`
   
   :type: float
   :default: 0.1


.. data:: edipy2.global_env.Uloc

   Values of the local interaction per orbital (max 5). 
   If less values are provided, the array is filled in increasing order
   
   Accesses :f:var:`uloc`
   
   :type: float
   :default: [2.0, 0.0, 0.0, 0.0, 0.0]
   
.. data:: edipy2.global_env.Ust

   Value of the inter-orbital interaction term
   
   Accesses :f:var:`ust`
   
   :type: float
   :default: 0.0
   
.. data:: edipy2.global_env.wini

   Value of the smallest real-axis frequency
   
   Accesses :f:var:`wini`
   
   :type: float
   :default: -5.0
   
.. data:: edipy2.global_env.wfin

   Value of the largest real-axis frequency
   
   Accesses :f:var:`wfin`
   
   :type: float
   :default: -5.0
   
.. data:: edipy2.global_env.xmin

   Value for the smallest position for the lattice PDF
   
   Accesses :f:var:`xmin`
   
   :type: float
   :default: -3.0

.. data:: edipy2.global_env.xmax

   Value for the largest position for the lattice PDF
   
   Accesses :f:var:`xmax`
   
   :type: float
   :default: 3.0

   
.. data:: edipy2.global_env.xmu

   Value of the chemical potential. If :f:var:`hfmode` = :code:`T`, :f:var:`xmu` = :code:`0` satisfied the half-filling condition
   
   Accesses :f:var:`xmu`
   
   :type: float
   :default: 0.0

