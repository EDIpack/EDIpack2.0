MODULE ED_INPUT_VARS
  !
  !Contains all global input variables which can be set by the user through the input file. A specific preocedure :f:func:`ed_read_input` should be called to read the input file using :f:func:`parse_input_variable` procedure from SciFortran. All variables are automatically set to a default, looked for and updated by reading into the file and, sequentially looked for and updated from command line (std.input) using the notation `variable_name=variable_value(s)` (case independent).  
  !
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str,free_unit
  USE ED_VERSION
  use iso_c_binding
  implicit none


  !input variables
  !=========================================================
  integer(c_int), bind(c, name="Nbath")                              :: Nbath               !Number of bath sites (per orbital or not depending on bath_type)
  integer(c_int), bind(c, name="Norb")                               :: Norb                !Number of impurity orbitals
  integer(c_int), bind(c, name="Nspin")                              :: Nspin               !Number spin degeneracy (max 2)

  integer(c_int), bind(c, name="Nloop")                              :: Nloop               !max dmft loop variables
  integer(c_int), bind(c, name="Nph")                                :: Nph                 !max number of phonons allowed (cut off)
  real(c_double),dimension(5),bind(c, name="Uloc")                   :: Uloc                !local interactions
  real(c_double),bind(c, name="Ust")                                 :: Ust                 !intra-orbitals interactions
  real(c_double),bind(c, name="Jh")                                  :: Jh                  !J_Hund: Hunds' coupling constant 
  real(c_double),bind(c, name="Jx")                                  :: Jx                  !J_X: coupling constant for the spin-eXchange interaction term
  real(c_double),bind(c, name="Jp")                                  :: Jp                  !J_P: coupling constant for the Pair-hopping interaction term 
  real(c_double),bind(c, name="xmu")                                 :: xmu                 !chemical potential
  real(c_double),bind(c, name="beta")                                :: beta                !inverse temperature
  !
  !
  integer(c_int), bind(c, name="Nsuccess")                           :: Nsuccess            !Number of repeated success to fall below convergence threshold  
  real(c_double),bind(c, name="dmft_error")                          :: dmft_error          !dmft convergence threshold
  real(c_double),bind(c, name="eps")                                 :: eps                 !broadening
  real(c_double),bind(c, name="wini")                                :: wini                !frequency range min
  real(c_double),bind(c, name="wfin")                                :: wfin                !frequency range max
  real(c_double),bind(c, name="xmin")                                :: xmin                !x-range for the local lattice probability distribution function (phonons)
  real(c_double),bind(c, name="xmax")                                :: xmax                !x-range for the local lattice probability distribution function (phonons)
  real(c_double),bind(c, name="sb_field")                            :: sb_field            !symmetry breaking field
  real(c_double),bind(c, name="nread")                               :: nread               !fixed density. if 0.d0 fixed chemical potential calculation.
  !
  !these variable need an equivalent "internal" one because parse_input_variable is not capable of reading logical(c_bool)
  logical(c_bool),bind(c, name="ed_total_ud")                        :: ed_total_ud         !flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw
  logical                                                            :: ed_total_ud_
  logical(c_bool),bind(c, name="ed_twin")                            :: ed_twin             !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
  logical                                                            :: ed_twin_
  !
  logical              :: HFmode              !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)              :: cutoff              !cutoff for spectral summation
  real(8)              :: gs_threshold        !Energy threshold for ground state degeneracy loop up
  real(8)              :: deltasc             !breaking symmetry field
  !
  integer              :: ph_type             !shape of the e part of the e-ph interaction: 1=orbital occupation, 2=orbital hybridization
  real(8)              :: A_ph                !phonon field coupled to displacement operator (constant)
  complex(8),allocatable  :: g_ph(:,:)        !electron-phonon coupling constant all
  real(8)              :: w0_ph               !phonon frequency (constant)
  real(8),allocatable  :: g_ph_diag(:)        !electron-phonon coupling constant diagonal (density)
  !
  real(8),allocatable  :: spin_field_x(:)        !magnetic field per orbital coupling to X-spin component
  real(8),allocatable  :: spin_field_y(:)        !magnetic field per orbital coupling to Y-spin component
  real(8),allocatable  :: spin_field_z(:)        !magnetic field per orbital coupling to Z-spin component
  real(8),allocatable  :: pair_field(:)          !pair field per orbital coupling to s-wave order parameter component
  real(8),dimension(4) :: exc_field           !external field coupling to exciton order parameter
  !
  logical              :: rdm_flag              !evaluate impurity RDM
  logical              :: chispin_flag        !evaluate spin susceptibility
  logical              :: chidens_flag        !evaluate dens susceptibility
  logical              :: chipair_flag        !evaluate pair susceptibility
  logical              :: chiexct_flag        !evaluate excitonic susceptibility
  !
  character(len=7)     :: ed_mode             !flag to set ed symmetry type: normal=normal (default), superc=superconductive, nonsu2=broken SU(2)
  logical              :: ed_finite_temp      !flag to select finite temperature method. note that if T then lanc_nstates_total must be > 1 
  logical              :: ed_sparse_H         !flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--
  logical              :: ed_solve_offdiag_gf !flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal 
  logical              :: ed_print_Sigma      !flag to print impurity Self-energies
  logical              :: ed_print_G          !flag to print impurity Green`s functions
  logical              :: ed_print_G0         !flag to print impurity non-interacting Green`s functions
  logical              :: ed_print_chispin    !flag to print impurity spin susceptibility
  logical              :: ed_print_chidens    !flag to print impurity dens susceptibility
  logical              :: ed_print_chipair    !flag to print impurity pair susceptibility
  logical              :: ed_print_chiexct    !flag to print impurity exct susceptibility
  logical              :: ed_all_G            !flag to evaluate all the components of the impurity Green`s functions irrespective of the symmetries
  logical              :: ed_sectors          !flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift
  integer              :: ed_sectors_shift    !shift to the ed_sectors scan
  integer              :: ed_verbose          !verbosity level: 0=almost nothing --> 5:all. Really: all
  real(8)              :: ed_offset_bath      !half-bandwidth for the bath initialization: flat in -hwband:hwband
  real(8)              :: ed_hw_bath          !half-bandwidth for the bath initialization: flat in -hwband:hwband
  logical              :: ed_obs_all          !flag to print observables for every loop
  !
  character(len=12)    :: lanc_method         !select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only) 
  real(8)              :: lanc_tolerance      !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos. 
  integer              :: lanc_niter          !Max number of Lanczos iterations
  integer              :: lanc_ngfiter        !Max number of iteration in resolvant tri-diagonalization
  integer              :: lanc_ncv_factor     !Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_ncv_add        !Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_nstates_sector !Max number of required eigenvalues per sector
  integer              :: lanc_nstates_total  !Max number of states hold in the finite T calculation
  integer              :: lanc_nstates_step   !Number of states added at each step to determine the optimal spectrum size at finite T
  integer              :: lanc_dim_threshold  !Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.
  !
  character(len=5)     :: cg_Scheme           !fit scheme: delta (default), weiss for G0
  integer              :: cg_method           !fit routine type:0=CGnr (default), 1=minimize (old f77)
  integer              :: cg_grad             !gradient evaluation: 0=analytic, 1=numeric
  integer              :: cg_Niter            !Max number of iteration in the fit
  real(8)              :: cg_Ftol             !Tolerance in the cg fit
  integer              :: cg_stop             !fit stop condition:0-3, 0=C1.AND.C2, 1=C1, 2=C2 with C1= :math:`\vert F_{n-1} -F_{n} \vert < tol*(1+F_{n})`, C2= :math:`\vert\vert x_{n-1} -x_{n} \vert\vert <tol*(1+ \vert\vert x_{n} \vert\vert`).
  integer              :: cg_Weight           !CGfit mode 0=1, 1=1/n , 2=1/w_n weight
  integer              :: cg_pow              !fit power to generalize the distance as  :math:`\vert G0 - G0_{and} \vert ^{cg\_pow}`
  character(len=12)    :: cg_norm             !frobenius/elemental (for now only in general bath)
  logical              :: cg_minimize_ver     !flag to pick old (Krauth) or new (Lichtenstein) version of the minimize CG routine
  real(8)              :: cg_minimize_hh      !unknown parameter used in the CG minimize procedure.  
  !
  logical              :: finiteT             !flag for finite temperature calculation
  character(len=7)     :: bath_type           !flag to set bath type: normal (1bath/imp), hybrid(1bath), replica(1replica/imp), general(replica++)
  !
  real(8)              :: nerr                !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)              :: ndelta              !initial chemical potential step
  real(8)              :: ncoeff              !multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff)
  integer              :: niter               !
  logical              :: Jz_basis            !"Flag to enable the Jz basis"
  logical              :: Jz_max              !"Flag to enable a maximum value for Jz"
  real(8)              :: Jz_max_value        !"Maximum value for Jz"

  !Some parameters for function dimension:
  integer(c_int),bind(c, name="Lmats")             :: Lmats !Number of Matsubara frequencies
  integer(c_int),bind(c, name="Lreal")             :: Lreal !Number of real-axis frequencies
  integer(c_int),bind(c, name="Lfit")              :: Lfit  !Number of frequencies for bath fitting

  integer(c_int),bind(c, name="Ltau")              :: Ltau  !Number of imaginary time points
  integer(c_int),bind(c, name="Lpos")              :: Lpos  !Number of points in PDF lattice

  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=100)   :: Hfile  !File where to retrieve/store the bath parameters.
  character(len=100)   :: HLOCfile !File read the input local H
  character(len=100)   :: SectorFile !File where to retrieve/store the sectors contributing to the spectrum
  character(len=100)   :: GPHfile !File of Phonon couplings. Put NONE to use only density couplings.
  integer(c_int),bind(c, name="LOGfile"),save             :: LOGfile  !Logfile unit

  !THIS IS JUST A RELOCATED GLOBAL VARIABLE
  character(len=200)                                 :: ed_input_file="" !Name of input file


contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine ed_read_input(INPUTunit)
    !This functions reads the input file provided by :code:`INPUTunit` and sets the global variables accordingly
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    logical          :: master=.true.,bool
    integer          :: i,iorb,rank=0
    integer          :: unit_xmu, unit_gph
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A,A)")"DEBUG ed_read_input: read input from",trim(INPUTunit)
#endif
#ifdef _MPI    
    if(check_MPI())then
       master=get_Master_MPI()
       rank  =get_Rank_MPI()
    endif
#endif
    !
    !Store the name of the input file:
    ed_input_file=str(INPUTunit)
    !
    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals (max 5).")
    call parse_input_variable(Nbath,"NBATH",INPUTunit,default=6,comment="Number of bath sites:(normal=>Nbath per orb)(hybrid=>Nbath total)(replica/general=>Nbath=Nreplica/Ngeneral)")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(Nph,"NPH",INPUTunit,default=0,comment="Max number of phonons allowed (cut off)")
    call parse_input_variable(bath_type,"BATH_TYPE",INPUTunit,default='normal',comment="flag to set bath type: normal (1bath/imp), hybrid(1bath), replica(1replica/imp), general(replica++)")
    !
    !allocate(Uloc(Norb)) #TODO: put me back!
    !call parse_input_variable(uloc,"ULOC",INPUTunit,default=(/( 2d0,i=1,size(Uloc) )/),comment="Values of the local interaction per orbital")
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2d0,0d0,0d0,0d0,0d0],comment="Values of the local interaction per orbital (max 5)")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    !
    call parse_input_variable(nloop,"NLOOP",INPUTunit,default=100,comment="Max number of DMFT iterations.")
    call parse_input_variable(nsuccess,"NSUCCESS",INPUTunit,default=1,comment="Number of successive iterations below threshold for convergence")
    call parse_input_variable(dmft_error,"DMFT_ERROR",INPUTunit,default=0.00001d0,comment="Error threshold for DMFT convergence")
    call parse_input_variable(sb_field,"SB_FIELD",INPUTunit,default=0.1d0,comment="Value of a symmetry breaking field for magnetic solutions.")
    call parse_input_variable(deltasc,"DELTASC",INPUTunit,default=0.02d0,comment="Value of the SC symmetry breaking term.")

    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")

    if(allocated(g_ph))deallocate(g_ph)
    if(allocated(g_ph_diag))deallocate(g_ph_diag)
    allocate(g_ph(Norb,Norb)) ! THIS SHOULD BE A MATRIX Norb*Norb
    allocate(g_ph_diag(Norb)) ! THIS SHOULD BE A MATRIX Norb*Norb
    call parse_input_variable(g_ph_diag,"G_PH",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="Electron-phonon coupling density constant")
    call parse_input_variable(w0_ph,"W0_PH",INPUTunit,default=0.d0,comment="Phonon frequency")
    call parse_input_variable(A_ph,"A_PH",INPUTunit,default=0.d0,comment="Forcing field coupled to phonons displacement operator")
    call parse_input_variable(GPHfile,"GPHfile",INPUTunit,default="NONE",comment="File of Phonon couplings. Put NONE to use only density couplings.")

    if(allocated(spin_field_x))deallocate(spin_field_x)
    if(allocated(spin_field_y))deallocate(spin_field_y)
    if(allocated(spin_field_z))deallocate(spin_field_z)
    if(allocated(pair_field))deallocate(pair_field)
    allocate(spin_field_x(Norb))
    allocate(spin_field_y(Norb))
    allocate(spin_field_z(Norb))
    allocate(pair_field(Norb))
    call parse_input_variable(spin_field_x,"SPIN_FIELD_X",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="magnetic field per orbital coupling to X-spin component")
    call parse_input_variable(spin_field_y,"SPIN_FIELD_Y",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="magnetic field per orbital coupling to Y-spin component")
    call parse_input_variable(spin_field_z,"SPIN_FIELD_Z",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="magnetic field per orbital coupling to Z-spin component")
    call parse_input_variable(exc_field,"EXC_FIELD",INPUTunit,default=[0d0,0d0,0d0,0d0],comment="external field coupling to exciton order parameters")
    call parse_input_variable(pair_field,"PAIR_FIELD",INPUTunit,default=(/( 0d0,i=1,Norb )/),comment="pair field per orbital coupling to s-wave order parameter component")
    !
    call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(chidens_flag,"CHIDENS_FLAG",INPUTunit,default=.false.,comment="Flag to activate density susceptibility calculation.")
    call parse_input_variable(chipair_flag,"CHIPAIR_FLAG",INPUTunit,default=.false.,comment="Flag to activate pair susceptibility calculation.")
    call parse_input_variable(chiexct_flag,"CHIEXCT_FLAG",INPUTunit,default=.false.,comment="Flag to activate excitonis susceptibility calculation.")
    !
    call parse_input_variable(ed_mode,"ED_MODE",INPUTunit,default='normal',comment="Flag to set ED type: normal=normal, superc=superconductive, nonsu2=broken SU(2)")
    call parse_input_variable(ed_finite_temp,"ED_FINITE_TEMP",INPUTunit,default=.false.,comment="flag to select finite temperature method. note that if T then lanc_nstates_total must be > 1")
    call parse_input_variable(ed_sectors,"ED_SECTORS",INPUTunit,default=.false.,comment="flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift.")
    call parse_input_variable(ed_sectors_shift,"ED_SECTORS_SHIFT",INPUTunit,1,comment="shift to ed_sectors")
    call parse_input_variable(ed_sparse_H,"ED_SPARSE_H",INPUTunit,default=.true.,comment="flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE ")
    !
    call parse_input_variable(ed_total_ud_,"ED_TOTAL_UD",INPUTunit,default=.true.,comment="flag to select which type of quantum numbers have to be considered: T (default) total Nup-Ndw, F orbital based Nup-Ndw")
    ed_total_ud = ed_total_ud_
    call parse_input_variable(ed_twin_,"ED_TWIN",INPUTunit,default=.false.,comment="flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.")
    ed_twin = ed_twin_
    call parse_input_variable(ed_obs_all,"ED_OBS_ALL",INPUTunit,default=.true.,comment="flag to print observables for every loop.")
    !
    call parse_input_variable(ed_solve_offdiag_gf,"ED_SOLVE_OFFDIAG_GF",INPUTunit,default=.false.,comment="flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal") 
    call parse_input_variable(ed_print_Sigma,"ED_PRINT_SIGMA",INPUTunit,default=.true.,comment="flag to print impurity Self-energies")
    call parse_input_variable(ed_print_G,"ED_PRINT_G",INPUTunit,default=.true.,comment="flag to print impurity Greens function")
    call parse_input_variable(ed_print_G0,"ED_PRINT_G0",INPUTunit,default=.true.,comment="flag to print non-interacting impurity Greens function")
    call parse_input_variable(ed_print_chispin,"ED_PRINT_CHISPIN",INPUTunit,default=.true.,comment="flag to print impurity spin susceptibility")
    call parse_input_variable(ed_print_chidens,"ED_PRINT_CHIDENS",INPUTunit,default=.true.,comment="flag to print impurity dens susceptibility")
    call parse_input_variable(ed_print_chipair,"ED_PRINT_CHIPAIR",INPUTunit,default=.true.,comment="flag to print impurity pair susceptibility")
    call parse_input_variable(ed_print_chiexct,"ED_PRINT_CHIEXCT",INPUTunit,default=.true.,comment="flag to print impurity exct susceptibility")
    call parse_input_variable(ed_all_G,"ED_ALL_G",INPUTunit,default=.true.,comment="flag to evaluate all the components of the impurity Green`s functions irrespective of the symmetries")
    call parse_input_variable(ed_verbose,"ED_VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 5:all. Really: all")
    call parse_input_variable(ed_hw_bath,"ed_hw_bath",INPUTunit,default=2d0,comment="half-bandwidth for the bath initialization: flat in -ed_hw_bath:ed_hw_bath")
    call parse_input_variable(ed_offset_bath,"ed_offset_bath",INPUTunit,default=1d-1,comment="offset for the initialization of diagonal terms in replica/general bath: -offset:offset")

    !
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=4096,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1024,comment="Number of imaginary time points.")
    call parse_input_variable(Lfit,"LFIT",INPUTunit,default=1000,comment="Number of Matsubara frequencies used in the \Chi2 fit.")
    call parse_input_variable(Lpos,"LPOS",INPUTunit,default=100,comment="Number of points for the lattice PDF.")
    !
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0,comment="Objective density for fixed density calculations.")
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4,comment="Error threshold for fixed density calculations.")
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0,comment="Initial step for fixed density calculations.")
    call parse_input_variable(ncoeff,"NCOEFF",INPUTunit,default=1d0,comment="multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff).")
    !
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    call parse_input_variable(xmin,"XMIN",INPUTunit,default=-3.d0,comment="Smallest position for the lattice PDF")
    call parse_input_variable(xmax,"XMAX",INPUTunit,default=3.d0,comment="Largest position for the lattice PDF")
    call parse_input_variable(rdm_flag,"RDM_FLAG",INPUTunit,default=.true.,comment="Flag to activate RDM calculation.")
    call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(chidens_flag,"CHIDENS_FLAG",INPUTunit,default=.false.,comment="Flag to activate density susceptibility calculation.")
    call parse_input_variable(chipair_flag,"CHIPAIR_FLAG",INPUTunit,default=.false.,comment="Flag to activate pair susceptibility calculation.")
    call parse_input_variable(chiexct_flag,"CHIEXCT_FLAG",INPUTunit,default=.false.,comment="Flag to activate excitonis susceptibility calculation.")
    !
    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    !    
    call parse_input_variable(lanc_method,"LANC_METHOD",INPUTunit,default="arpack",comment="select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only), DVDSON (no MPI)")
    call parse_input_variable(lanc_nstates_sector,"LANC_NSTATES_SECTOR",INPUTunit,default=2,comment="Initial number of states per sector to be determined.")
    call parse_input_variable(lanc_nstates_total,"LANC_NSTATES_TOTAL",INPUTunit,default=1,comment="Initial number of total states to be determined.")
    call parse_input_variable(lanc_nstates_step,"LANC_NSTATES_STEP",INPUTunit,default=2,comment="Number of states added to the spectrum at each step.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=10,comment="Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=0,comment="Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,comment="Number of Lanczos iteration in spectrum determination.")
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=1d-18,comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,default=1024,comment="Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.")
    !
    call parse_input_variable(cg_method,"CG_METHOD",INPUTunit,default=0,comment="Conjugate-Gradient method: 0=NumericalRecipes, 1=minimize.")
    call parse_input_variable(cg_grad,"CG_GRAD",INPUTunit,default=0,comment="Gradient evaluation method: 0=analytic (default), 1=numeric.")
    call parse_input_variable(cg_ftol,"CG_FTOL",INPUTunit,default=0.00001d0,comment="Conjugate-Gradient tolerance.")
    call parse_input_variable(cg_stop,"CG_STOP",INPUTunit,default=0,comment="Conjugate-Gradient stopping condition: 0-2, 0=C1.AND.C2, 1=C1, 2=C2 with C1=|F_n-1 -F_n|<tol*(1+F_n), C2=||x_n-1 -x_n||<tol*(1+||x_n||).")
    call parse_input_variable(cg_niter,"CG_NITER",INPUTunit,default=500,comment="Max. number of Conjugate-Gradient iterations.")
    call parse_input_variable(cg_weight,"CG_WEIGHT",INPUTunit,default=1,comment="Conjugate-Gradient weight form: 1=1.0, 2=1/n , 3=1/w_n.")
    call parse_input_variable(cg_scheme,"CG_SCHEME",INPUTunit,default='weiss',comment="Conjugate-Gradient fit scheme: delta or weiss.")
    call parse_input_variable(cg_norm,"CG_NORM",INPUTunit,default='elemental',comment="Conjugate-Gradient norm definition: elemental (default) or frobenius.")
    call parse_input_variable(cg_pow,"CG_POW",INPUTunit,default=2,comment="Fit power for the calculation of the generalized distance as |G0 - G0and|**cg_pow")
    call parse_input_variable(cg_minimize_ver,"CG_MINIMIZE_VER",INPUTunit,default=.false.,comment="Flag to pick old/.false. (Krauth) or new/.true. (Lichtenstein) version of the minimize CG routine")
    call parse_input_variable(cg_minimize_hh,"CG_MINIMIZE_HH",INPUTunit,default=1d-4,comment="Unknown parameter used in the CG minimize procedure.")
    !
    call parse_input_variable(Jz_basis,"JZ_BASIS",INPUTunit,default=.false.,comment="Flag to enable the Jz basis")
    call parse_input_variable(Jz_max,"JZ_MAX",INPUTunit,default=.false.,comment="Whether to cutoff Jz")
    call parse_input_variable(Jz_max_value,"JZ_MAX_VALUE",INPUTunit,default=1000.d0,comment="Maximum Jz")
    !
    call parse_input_variable(SectorFile,"SectorFile",INPUTunit,default="sectors",comment="File where to retrieve/store the sectors contributing to the spectrum.")
    call parse_input_variable(Hfile,"Hfile",INPUTunit,default="hamiltonian",comment="File where to retrieve/store the bath parameters.")
    call parse_input_variable(HLOCfile,"HLOCfile",INPUTunit,default="inputHLOC.in",comment="File read the input local H.")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")

    if(nph>0)then
       !Here the non-diagonal (non-density) phononic coupling are read
       g_ph=0.d0
       if(trim(GPHfile).eq."NONE")then
          do iorb=1,Norb
             g_ph(iorb,iorb)=g_ph_diag(iorb)
          enddo
       else
          inquire(file=trim(GPHfile),EXIST=bool)
          if(bool)then
             open(free_unit(unit_gph),file=GPHfile)
             do iorb=1,Norb
                read(unit_gph,*) g_ph(iorb,:)
             enddo
             close(unit_gph)
             !maybe an assert_hermitian would be globally useful
             if(any(g_ph /= transpose(conjg(g_ph))))then
                stop "ERROR: non hermitian phonon coupling matrix (g_ph) in input"
             end if
          else
             stop "GPHfile/=NONE but there is no GPHfile with the provided name"
          endif
       endif
       !TO BE PUT SOMEWHERE ELSE
       open(free_unit(unit_gph),file="GPHinput.used")
       do iorb=1,Norb
          write(unit_gph,*) g_ph(iorb,:)
       enddo
       close(unit_gph)
    end if

#ifdef _MPI
    if(check_MPI())then
       if(.not.master)then
          LOGfile=1000-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".ed")
          do i=1,get_Size_MPI()
             if(i==rank)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    !
    Ltau=max(int(beta),Ltau)
    if(master)then
       call print_input()
       call save_input(INPUTunit)
       call scifor_version()
       call code_version(version)
    endif
    !
    if(nread .ne. 0d0) then
       inquire(file="xmu.restart",EXIST=bool)
       if(bool)then
          open(free_unit(unit_xmu),file="xmu.restart")
          read(unit_xmu,*)xmu,ndelta
          ndelta=abs(ndelta)*ncoeff
          close(unit_xmu)
          write(*,"(A,F9.7,A)")"Adjusting XMU to ",xmu," as per provided xmu.restart "
       endif
    endif
    !Act on the input variable only after printing.
    !In the new parser variables are hard-linked into the list:
    !any change to the variable is immediately copied into the list... (if you delete .ed it won't be printed out)
    call substring_delete(Hfile,".restart")
    call substring_delete(Hfile,".ed")
  end subroutine ed_read_input

  subroutine ed_update_input(name,vals)
    !This functions updates some variables in the input file, namely
    !:f:var:`exc_field`, :f:var:`pair_field`, :f:var:`exc_field`,
    !:f:var:`spin_field_x`, :f:var:`spin_field_y`, and :f:var:`spin_field_z`.
    character(len=*)      :: name !the name of the variable to update
    real(8),dimension(:)  :: vals !the new value of the variable
    select case (name)
    case default
       stop "WRONG NAME ON ED_UPDATE_INPUT"
    case ("EXC_FIELD")
       if(size(vals)/=4)stop "WRONG SIZE IN ED_UPDATE_EXC_FIELD"
       exc_field=vals
    case ("PAIR_FIELD")
       if(size(vals)/=Norb)stop "WRONG SIZE IN ED_UPDATE_PAIR_FIELD"
       pair_field=vals
    case ("SPIN_FIELD_X")
       if(size(vals)/=Norb)stop "WRONG SIZE IN ED_UPDATE_SPIN_FIELD_X"
       spin_field_x=vals
    case ("SPIN_FIELD_Y")
       if(size(vals)/=Norb)stop "WRONG SIZE IN ED_UPDATE_SPIN_FIELD_Y"
       spin_field_y=vals
    case ("SPIN_FIELD_Z")
       if(size(vals)/=Norb)stop "WRONG SIZE IN ED_UPDATE_SPIN_FIELD_Z"
       spin_field_z=vals
    end select
    
  end subroutine ed_update_input




  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE ED_INPUT_VARS
