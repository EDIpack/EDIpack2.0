MODULE E2I_BATH_REPLICA
  !
  !:synopsis: Replica bath construction and manipulation routines, inequivalent sites version
  !This module implements the functions to set the matrix basis :math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` and the initial variational parameters :math:`\vec{\lambda}` used to decompose each local bath hamiltonian for the  :f:var:`replica` and :f:var:`general` bath types.  
  !The functions extends the functionalities to the case of multiple inequivalent impurities.
  !
  USE EDIPACK2
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





  interface set_Hreplica
     ! This function sets the matrix basis :math:`\{ \hat{O}_i \}_{i=1,\dots,N_{sym}}` used to decompose the single bath hamiltonian :math:`h^p`. It also sets the initial values of the variational parameters :math:`\vec{\lambda}` in the :f:var:`replica` ath type.
     !
     ! Input: :f:var:`Hvec` 
     !  * rank-5, dimensions: [ |Nns| , |Nns|, |Norb| , |Norb| , |Nsym| ]
     !  * rank-3, dimensions: [ |Nnso| , |Nnso|, |Nsym| ]
     !
     ! Input: :f:var:`lambdavec`
     !  * rank-2, dimensions: [ |Nbath| , |Nsym| ]
     !  * rank-3, dimensions: [ |Nlat| , |Nbath| , |Nsym| ]
     !
     module procedure init_Hreplica_symmetries_lattice_d5
     module procedure init_Hreplica_symmetries_lattice_d3
  end interface set_Hreplica

  interface set_Hgeneral
     !
     !A clone of  :f:var:`set_Hreplica`
     !
     module procedure init_Hreplica_symmetries_lattice_d5
     module procedure init_Hreplica_symmetries_lattice_d3
  end interface set_Hgeneral


  real(8),dimension(:,:,:),allocatable :: Hlambda_ineq  !rank-3 array, dimensions: [ |Nlat| , |Nbath| , |Nsym| ]. Stores the initial values of parameters the :math:`\lambda_i` for the baths matrix decomposition

  !Internal use
  public :: Hlambda_ineq
  public :: Hreplica_site
  !
  !To be published
  public :: set_Hreplica
  public :: set_Hgeneral


  integer :: ibath

contains



  subroutine init_Hreplica_symmetries_lattice_d5(Hvec,lambdavec)
    complex(8),dimension(:,:,:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:,:,:)        :: lambdavec ![Nlat,Nbath,Nsym]
    integer                         :: ilat,Nlat
    integer                         :: isym,Nsym
    logical                         :: bool
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
#endif
    if(ed_mode=="superc")Nnambu=2
    Nlat=size(lambdavec,1)
    Nsym=size(lambdavec,3)
    !
    call assert_shape(Hvec,[Nnambu*Nspin,Nnambu*Nspin,Norb,Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default  ;bool = check_herm(nn2so_reshape(Hvec(:,:,:,:,isym),Nspin,Norb),Nspin*Norb)
       case("superc");bool = check_nambu(nn2so_reshape(Hvec(:,:,:,:,isym),Nnambu*Nspin,Norb),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    call ed_allocate_Hreplica(Nsym) !< allocate the memory for EDIpack.ED_BATH.Hb
    !
    if(allocated(Hlambda_ineq))deallocate(Hlambda_ineq)
    allocate(Hlambda_ineq(Nlat,Nbath,Nsym))
    do isym=1,Nsym
       call ed_set_Hsym_Hreplica(isym,Hvec(:,:,:,:,isym))
       Hlambda_ineq(:,:,isym) = lambdavec(:,:,isym)
    enddo
    !
    if(ed_verbose>2)then
       do ilat=1,Nlat
          write(LOGfile,"(A)")"Inequivalent #"//str(ilat)//":"
          do ibath=1,Nbath
             write(LOGfile,"(A)")"> Hreplica #"//str(ibath)//":"
             call ed_print_Hreplica( ed_build_Hreplica( Hlambda_ineq(ilat,ibath,:) ) )
          enddo
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_lattice_d5

  subroutine init_Hreplica_symmetries_lattice_d3(Hvec,lambdavec)
    complex(8),dimension(:,:,:) :: Hvec      ![size(H),Nsym]
    real(8),dimension(:,:,:)    :: lambdavec ![Nlat,Nbath,Nsym]
    integer                     :: ilat,Nlat
    integer                     :: isym,Nsym
    logical                     :: bool
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG init_Hreplica_symmetries_lattice: from ({[Hs,Lam]}_b)_site"
#endif
    if(ed_mode=="superc")Nnambu=2
    Nlat=size(lambdavec,1)
    Nsym=size(lambdavec,3)
    !
    call assert_shape(Hvec,[Nnambu*Nspin*Norb,Nnambu*Nspin*Norb,Nsym],"init_Hreplica_symmetries","Hvec")
    !
    !CHECK NAMBU and HERMITICTY of each Hvec
    do isym=1,Nsym
       select case(ed_mode)
       case default
          bool = check_herm(Hvec(:,:,isym),Nspin*Norb)
       case("superc")
          bool = check_nambu(Hvec(:,:,isym),Nspin*Norb)
       end select
       if(.not.bool)then
          write(LOGfile,"(A)")"init_Hreplica_symmetries_site ERROR: not Hermitian/Nambu of replica basis O_"//str(isym)
          stop
       endif
    enddo
    !
    if(allocated(Hlambda_ineq))deallocate(Hlambda_ineq)
    allocate(Hlambda_ineq(Nlat,Nbath,Nsym))
    call ed_allocate_Hreplica(Nsym)
    !
    do isym=1,Nsym
       Hlambda_ineq(:,:,isym) = lambdavec(:,:,isym)
       call ed_set_Hsym_Hreplica(isym, so2nn_reshape(Hvec(:,:,isym),Nnambu*Nspin,Norb) )
    enddo
    !
    if(ed_verbose>2)then
       do ilat=1,Nlat
          write(LOGfile,"(A)")"Inequivalent #"//str(ilat)//":"
          do ibath=1,Nbath
             write(LOGfile,"(A)")"> Hreplica #"//str(ibath)//":"
             call ed_print_Hreplica( ed_build_Hreplica( Hlambda_ineq(ilat,ibath,:) ) )
          enddo
       enddo
    endif
    !
  end subroutine init_Hreplica_symmetries_lattice_d3





  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################



  subroutine Hreplica_site(site)
    !
    !Set the current initial value of the replica bath parameters to the current inequivalent impurity :math:`\lambda_0=\lambda_{site}` 
    !
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    integer :: site             !the index of the inequivalent impurity
    if(site<1.OR.site>size(Hlambda_ineq,1))stop "ERROR Hreplica_site: site not in [1,Nlat]"
    if(.not.allocated(Hlambda_ineq))stop "ERROR Hreplica_site: Hreplica_lambda_ineq not allocated"
    call ed_set_linit_Hreplica(Hlambda_ineq(site,:,:))
  end subroutine Hreplica_site




  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################
  !##################################################################



  function check_herm(A,N,error) result(bool)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    integer,intent(in)                   :: N
    complex(8),dimension(N,N),intent(in) :: A
    logical                              :: bool
    real(8),optional                     :: error
    real(8)                              :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    bool   = all(abs(A - conjg(transpose(A)))<error_)
  end function check_herm


  function check_nambu(A,N,error) result(bool)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    integer,intent(in)                       :: N
    complex(8),dimension(2*N,2*N),intent(in) :: A
    complex(8),dimension(N,N)                :: h11,h22
    logical                                  :: bool
    real(8),optional                         :: error
    real(8)                                  :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    h11    = A(1:N    ,1:N)
    h22    = A(N+1:2*N,N+1:2*N)
    bool   = check_herm(A,2*N,error_) !this checks also for F = A_12, s.t. A_21=herm(A_12)
    bool   = bool.AND.( all(abs(h22 + conjg(h11))<error_) )
  end function check_nambu







END MODULE E2I_BATH_REPLICA
