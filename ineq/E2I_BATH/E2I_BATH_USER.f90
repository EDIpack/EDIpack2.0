MODULE E2I_BATH_USER
  USE EDIPACK2
  !Implements functions the user can use to enforce specific symmetry operations on the bath array.
  !
  USE E2I_VARS_GLOBAL
  USE E2I_AUX_FUNX
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  implicit none

  private

  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  !explicit symmetries:
  interface break_symmetry_bath
     !
     ! Function to impose a specific symmetry breaking pattern into the energy levels of the bath. A common case is to find magnetic solution by breaking spin degeneracy of the levels.
     !
     ! * :f:var:`bath` is  rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     !
     module procedure break_symmetry_bath_lattice
  end interface break_symmetry_bath

  interface spin_symmetrize_bath
     !
     ! Function to impose a spin symmetry to the parameters of the bath. Enforces a non-magnetic solution
     !
     ! * :f:var:`bath` is  rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure spin_symmetrize_bath_lattice
  end interface spin_symmetrize_bath

  interface orb_symmetrize_bath
     !
     ! Function to impose a orbital symmetry to the parameters of the bath. Enforces an orbital non-polarized solution. If two orbital indices :f:var:`orb1` and :f:var:`orb2` are passed symmetry is imposed only among such two orbitals
     !
     !.. warning::
     !
     !   This operation requires the orbital to be degenerate.
     !
     ! * :f:var:`bath` is  rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure orb_symmetrize_bath_lattice
     module procedure orb_symmetrize_bath_lattice_o1o2
  end interface orb_symmetrize_bath

  interface orb_equality_bath
     !
     ! Function to impose a orbital equality on the parameters of the bath. 
     !
     ! * :f:var:`bath` is  rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure orb_equality_bath_lattice
  end interface orb_equality_bath

  interface ph_symmetrize_bath
     !
     ! Function to impose particle-hole symmetry to the parameters of the bath. 
     !
     ! * :f:var:`bath` is  rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure ph_symmetrize_bath_lattice
  end interface ph_symmetrize_bath

  interface ph_trans_bath
     !
     ! Function to perform particle-hole transformation to the parameters of the bath. 
     !
     ! * :f:var:`bath` is  rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure ph_trans_bath_lattice
  end interface ph_trans_bath

  interface enforce_normal_bath
     !
     ! Function to impose normal solution to the parameters of the bath, i.e. suppressed superconductivity if any.
     !
     ! * :f:var:`bath` is  rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure enforce_normal_bath_lattice
  end interface enforce_normal_bath

  interface save_array_as_bath
     !
     ! Write the bath parameters to a file following the convention of the internal data structure :f:var:`effective_bath`. 
     !
     ! * :f:var:`bath` is  rank-2 [:, :f:var:`Nlat` ] array type: double precision with assumed size
     !
     module procedure save_array_as_bath_lattice
  end interface save_array_as_bath




  !##################################################################
  !
  !     USER BATH ROUTINES:
  !
  !##################################################################
  public :: break_symmetry_bath
  public :: spin_symmetrize_bath
  public :: orb_symmetrize_bath
  public :: orb_equality_bath
  public :: ph_symmetrize_bath
  public :: ph_trans_bath
  public :: enforce_normal_bath
  public :: save_array_as_bath





contains


  !##################################################################
  !
  !     USER BATH  SYMMETRIES: PREDEFINED AND USER CONTROLLED
  !
  !##################################################################

  !+-------------------------------------------------------------------+
  !PURPOSE  : given a bath array apply a specific transformation or
  ! impose a given symmetry:
  ! - break spin symmetry by applying a symmetry breaking field
  ! - given a bath array set both spin components to have
  !    the same bath, i.e. impose non-magnetic solution
  ! - given a bath array enforces the particle-hole symmetry
  !    by setting the positive energies in modulo identical to the negative
  !    ones.
  ! - given a bath enforce normal (i.e. non superconducting) solution
  ! - given a dmft bath pull/push the components W^{ss'}_\a(l) of the Hybridization
  !    matrix
  ! - given a dmft bath pull/push the nonsu2 components
  !+-------------------------------------------------------------------+
  subroutine break_symmetry_bath_lattice(bath_,field,sign,save)
    real(8),dimension(:,:) :: bath_
    real(8)                :: field
    real(8)                :: sign
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_break_symmetry_bath(bath_(ilat,:),field,sign,save_)
    enddo
    ed_file_suffix=""
  end subroutine break_symmetry_bath_lattice


  !---------------------------------------------------------!

  !
  subroutine spin_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_spin_symmetrize_bath(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine spin_symmetrize_bath_lattice


  !---------------------------------------------------------!

  subroutine orb_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "orb_symmetry_bath ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_orb_symmetrize_bath(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_symmetrize_bath_lattice

  subroutine orb_symmetrize_bath_lattice_o1o2(bath_,orb1,orb2,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat,orb1,orb2
    if(bath_type=="replica")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_orb_symmetrize_bath(bath_(ilat,:),orb1,orb2,save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_symmetrize_bath_lattice_o1o2

  !---------------------------------------------------------!

  subroutine orb_equality_bath_lattice(bath_,indx,save)
    real(8),dimension(:,:) :: bath_
    integer,optional       :: indx
    logical,optional       :: save
    integer                :: indx_
    logical                :: save_
    integer                :: iorb
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "orb_equality_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "orb_equality_bath_site ERROR: can not be used with bath_type=general"
    indx_=1     ;if(present(indx))indx_=indx
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_orb_equality_bath(bath_(ilat,:),indx_,save_)
    enddo
    ed_file_suffix=""
  end subroutine orb_equality_bath_lattice



  !---------------------------------------------------------!

  subroutine ph_symmetrize_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_symmetry_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_ph_symmetrize_bath(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine ph_symmetrize_bath_lattice

  !---------------------------------------------------------!

  subroutine ph_trans_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    if(bath_type=="replica")stop "ph_trans_bath_site ERROR: can not be used with bath_type=replica"
    if(bath_type=="general")stop "ph_trans_bath_site ERROR: can not be used with bath_type=general"
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_ph_trans_bath(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine ph_trans_bath_lattice

  !---------------------------------------------------------!

  subroutine enforce_normal_bath_lattice(bath_,save)
    real(8),dimension(:,:) :: bath_
    logical,optional       :: save
    logical                :: save_
    integer                :: Nsites,ilat
    save_=.true.;if(present(save))save_=save
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_enforce_normal_bath(bath_(ilat,:),save_)
    enddo
    ed_file_suffix=""
  end subroutine enforce_normal_bath_lattice





  !+------------------------------------------------------------------+
  !PURPOSE  : given a array, save it as a bath. Do nothing else.
  !+------------------------------------------------------------------+
  subroutine save_array_as_bath_lattice(bath_)
    real(8),dimension(:,:) :: bath_
    integer                :: Nsites,ilat
    Nsites=size(bath_,1)
    do ilat=1,Nsites
       ed_file_suffix=reg(ineq_site_suffix)//reg(str(ilat,site_indx_padding))
       call ed_save_array_as_bath(bath_(ilat,:))
    enddo
    ed_file_suffix=""
  end subroutine save_array_as_bath_lattice




END MODULE E2I_BATH_USER













