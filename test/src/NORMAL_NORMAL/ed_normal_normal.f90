program ed_normal_normal
  USE EDIPACK2
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  implicit none
  integer                                     :: i,iw,jo,js,Nso,Nmomenta
  integer                                     :: unit,unit_
  real(8)                                     :: w, Im, Re
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,jspin,inso,print_mode
  real(8),allocatable                         :: Bath(:),Wlist(:)
  !GFs and Sigma:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8),dimension(:),allocatable            :: H0     ![Nso]
  !variables for the model:
  real(8)                                     :: Delta
  character(len=16)                           :: finput
  !NORMAL variables:
  real(8),allocatable                         :: dens(:),docc(:),energy(:),imp(:), Smats11mom(:)
  !CHECK variables:
  real(8),allocatable                         :: dens_(:),docc_(:),energy_(:),imp_(:),Smats11mom_(:)
  !
  !MPI Vars:
  integer                                     :: irank,comm,rank,size2,ierr
  logical                                     :: master
  !
  ! MPI initialization
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  size2 = get_Size_MPI(comm)
  master = get_Master_MPI(comm)
  !
  !Parse additional variables && read Input
  call parse_cmd_variable(finput,"FINPUT",default="inputED.in")
  call parse_input_variable(delta,"DELTA",finput,default=0.d0)
  !
  !
  call ed_read_input(trim(finput))
  !
  if(bath_type/="normal")stop "Wrong setup from input file: non normal bath"
  if(ed_mode/="normal")stop "Wrong setup from input file: non normal ed_mode"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=1 )stop "Wrong setup from input file: Nspin/=1"
  Nso=Nspin*Norb
  Nmomenta=4
  !
  !Allocate Weiss Field:
  allocate(Weiss(1,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(dens(Norb),dens_(Norb))
  allocate(docc(Norb),docc_(Norb))
  allocate(energy(8),energy_(8))
  allocate(imp(2),imp_(2))
  allocate(Wlist(size(Smats,5)))  
  allocate(Smats11mom(Nmomenta),Smats11mom_(Nmomenta))
  Wlist = pi/beta*(2*arange(1,Lmats)-1)
  !
  !
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(H0(Nso))
  Hloc = zero
  H0   = zero
  do js=1,Nspin
     Hloc(js,js,:,:)= Delta*pauli_sigma_z
     do jo=1,Norb
        H0(jo+2*(js-1)) =Hloc(js,js,jo,jo)
     end do
  end do
  !
  print_mode=3
  !
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  !
  !
  !SPARSE + HK
  call run_test(.true.,.false.)
  
  ! !SPARSE + UMATRIX
  ! call run_test(.true.,.true.)

  ! !DIRECT + HK
  ! call run_test(.false.,.false.)

  ! !DIRECT + UMATRIX
  ! call run_test(.false.,.true.)
  
  
  call finalize_MPI()



contains


  subroutine run_test(sparse,umatrix)
    logical :: sparse,umatrix
    ! ED_SPARSE_H    =sparse
    ! ED_READ_UMATRIX=umatrix
    ! LOGFILE=100
    ! if(sparse)LOGFILE=LOGFILE+10
    ! if(umatrix)LOGFILE=LOGFILE+1
    call ed_init_solver(bath)
    call ed_set_Hloc(hloc)
    write(*,*) ""
    write(*,*) "ED_MODE = NORMAL   |   BATH_TYPE = NORMAL"
    write(*,*) "SPARSE_H= "//str(sparse)//"        |   Umatrix = "//str(umatrix)
    write(*,*) "Checking..."
    call ed_solve(bath)
    call test_results()
    call ed_finalize_solver()
  end subroutine run_test



  

  subroutine test_results()
    !
    !get ED
    call ed_get_sigma(Smats,axis="m",type="n")
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    call ed_get_eimp(energy(1:4))
    call ed_get_doubles(energy(5:8))
    call ed_get_imp_info(imp)   
    do i=1,Nmomenta
       call compute_momentum(Wlist,Smats(1,1,1,1,:),i,Smats11mom(i))
    enddo
    !
    call read_array("dens_last.check",dens_)
    call read_array("docc_last.check",docc_)
    call read_array("energy_last.check",energy_)
    call read_array("imp_last.check",imp_)
    print*,imp
    print*,imp_
    call read_array("impSigma_l11_s1_iw.momenta.check",Smats11mom_)
    !
    !
    call assert(dens,dens_,"dens(:)")
    call assert(docc,docc_,"docc(:)")
    call assert(energy,energy_,"energy(:)")
    call assert(imp,imp_,"imp(:)")
    call assert(Smats11mom/Smats11mom_,dble(ones(Nmomenta)),"Sigma_matsubara_l11(:)",tol=1.0d-8)
    !
  end subroutine test_results







  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index


  function so2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so

  ! Subroutine to compute momenta
  ! 
  ! ( sum_w abs(F(w))*w**n ) / ( sum_w abs(F(w)) )
  subroutine compute_momentum(x,Fx,n,momentum)
    real(8)   ,dimension(:),intent(in)       :: x
    complex(8),dimension(:),intent(in)       :: Fx
    integer   ,intent(in)                    :: n
    real(8)   ,intent(out)                   :: momentum
    !
    integer                                  :: iw
    real(8)                                  :: num,den
    num=0.0; den=0.0
    do iw=1,size(x,1)
       num = num + abs(Fx(iw))*x(iw)**n
       den = den + abs(Fx(iw))
    enddo
    momentum=num/den
  end subroutine compute_momentum


end program



