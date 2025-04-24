program ed_hybrid_nonsu2
  USE EDIPACK2
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  implicit none
  integer                :: i,iw,jo,js,Nso,Nsymm,Nmomenta
  !Bath:
  integer                :: Nb,iorb,jorb,ispin,jspin,inso,print_mode
  real(8),allocatable    :: Bath(:),Wlist(:)
  !GFs and Sigma:
  complex(8),allocatable :: Smats(:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable :: Hloc(:,:,:,:)
  !variables for the model:
  real(8)                :: mh,lambda
  character(len=16)      :: finput
  real(8),allocatable    :: evals(:),evalsR(:)
  !Replica variables:
  real(8),allocatable    :: dens(:),docc(:),magX(:),energy(:),doubles(:),imp(:)
  real(8),allocatable    :: densR(:),doccR(:),magXR(:),energyR(:),doublesR(:),impR(:)
  real(8),allocatable    :: S11mom(:,:),S12mom(:,:)
  real(8),allocatable    :: S11momR(:,:),S12momR(:,:)
  !
  complex(8),dimension(4,4)                   :: Gamma1,Gamma2,Gamma5,GammaN,GammaS
  complex(8),dimension(4,4)                   :: GammaE0,GammaEx,GammaEy,GammaEz
  real(8),dimension(:),allocatable            :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable :: Hsym_basis
  !MPI Vars:
  integer :: irank,comm,rank,size2,ierr
  logical :: master
  logical :: dsave
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
  call parse_input_variable(mh,"MH",finput,default=1.d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_cmd_variable(dsave,"dsave",default=.false.)
  !
  !
  call ed_read_input(trim(finput))
  !
  if(bath_type/="hybrid")stop "Wrong setup from input file: non hybrid bath"
  if(ed_mode/='nonsu2')stop "Wrong setup from input file: ed_mode != nonsu2"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=2 )stop "Wrong setup from input file: Nspin/=2"
  Nso=Nspin*Norb
  Nmomenta=4
  !
  ! Matrices for replica hamiltonian
  gamma1=kron( pauli_sigma_z, pauli_tau_x)
  gamma2=kron( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron( pauli_sigma_0, pauli_tau_z)
  gammaN=kron( pauli_sigma_0, pauli_tau_0)
  !
  gammaE0=kron( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron( pauli_sigma_z, pauli_tau_x )
  !
  !Allocate Weiss Field:
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(dens(Norb))
  allocate(docc(Norb))
  allocate(magX(Norb))
  allocate(energy(4))
  allocate(doubles(4))
  allocate(imp(2))
  allocate(S11mom(Norb,Nmomenta))
  allocate(S12mom(Norb,Nmomenta))
  !
  allocate(Wlist(Lmats))  
  Wlist = pi/beta*(2*arange(1,Lmats)-1)
  !
  !
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc = j2so(Mh*Gamma5)
  !
  !
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  !
  !
  call run_test(sparse=.true.,umatrix=.false.)  
  call run_test(sparse=.false.,umatrix=.false.)
  call run_test(sparse=.true.,umatrix=.true.)  
  call run_test(sparse=.false.,umatrix=.true.)
  !
  call finalize_MPI()



contains

  subroutine run_test(sparse,umatrix)
    logical :: sparse,umatrix
    ED_SPARSE_H    =sparse
    ED_READ_UMATRIX=umatrix
    if(umatrix)ED_USE_KANAMORI=.false.
    call ed_init_solver(bath)
    call ed_set_Hloc(hloc)
    write(*,*) ""
    write(*,*) "ED_MODE = SUPERC   |   BATH_TYPE = HYBRID"
    write(*,*) "SPARSE_H= "//str(sparse)
    write(*,*) "U_MATRIX= "//str(umatrix)
    call ed_solve(bath)
    call ed_get_sigma(Smats,axis="m",type="n")
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    call ed_get_mag(magX,'x')
    call ed_get_eimp(energy)
    call ed_get_doubles(doubles)
    call ed_get_imp_info(imp)
    call ed_get_evals(evals)
    do i=1,Nmomenta
       do iorb=1,Norb
          call compute_momentum(Wlist,Smats(1,1,iorb,iorb,:),i,S11mom(iorb,i))
          call compute_momentum(Wlist,Smats(1,2,iorb,iorb,:),i,S12mom(iorb,i))
       enddo
    enddo
    call ed_finalize_solver()
    if(dsave)then
       write(*,*)"Saving results to .check files and exit"
       call save_results()
       stop
    endif
    call test_results(sparse,umatrix)
  end subroutine run_test



  subroutine read_results()
    integer :: L
    L = file_length("evals.check")
    if(allocated(evals))deallocate(evals)
    allocate(evals(L))
    call read_array("evals.check",evals)
    call read_array("dens.check",dens)
    call read_array("docc.check",docc)
    call read_array("magX.check",magX)
    call read_array("energy.check",energy)
    call read_array("doubles.check",doubles)
    call read_array("imp.check",imp)
    call read_array("Sigma11_momenta.check",S11mom)
    call read_array("Sigma12_momenta.check",S12mom)
    if(allocated(densR))deallocate(densR)
    if(allocated(doccR))deallocate(doccR)
    if(allocated(magXR))deallocate(magXR)
    if(allocated(energyR))deallocate(energyR)
    if(allocated(doublesR))deallocate(doublesR)
    if(allocated(impR))deallocate(impR)
    if(allocated(evalsR))deallocate(evalsR)
    if(allocated(S11momR))deallocate(S11momR)
    if(allocated(S12momR))deallocate(S12momR)
    allocate(densR, source=dens)
    allocate(doccR, source=docc)
    allocate(magXR, source=magX)
    allocate(energyR, source=energy)
    allocate(doublesR, source=doubles)
    allocate(impR, source=imp)
    allocate(evalsR,source=evals)
    allocate(S11momR, source=S11mom)
    allocate(S12momR, source=S12mom)
  end subroutine read_results


  subroutine test_results(sparse,umatrix)
    logical  :: sparse,umatrix
    write(*,*)""
    write(*,*) "Check RESULTS sparse_H, read_umatrix="//str(sparse)//","//str(umatrix)
    call read_results()
    call assert(dens,densR,"dens")
    call assert(docc,doccR,"docc")
    call assert(magX,magXR,"magX")
    call assert(energy,energyR,"energy")
    call assert(doubles,doublesR,"doubles")
    call assert(imp,impR,"imp")
    call assert(evals,evalsR,"evals")
    call assert(S11mom/S11momR,dble(ones(Norb,Nmomenta)),"Sigma11_momenta 1:4",tol=1.0d-8)
    call assert(S12mom/S12momR,dble(ones(Norb,Nmomenta)),"Sigma12_momenta 1:4",tol=1.0d-8)
    write(*,*)""
    write(*,*)""
    write(*,*)""
    write(*,*)""
    call wait(500)
  end subroutine test_results



  subroutine save_results()
    call save_array("evals.check",evals)
    call save_array("dens.check",dens)
    call save_array("docc.check",docc)
    call save_array("magX.check",magX)
    call save_array("energy.check",energy)
    call save_array("doubles.check",doubles)
    call save_array("imp.check",imp)
    call save_array("Sigma11_momenta.check",S11mom)
    call save_array("Sigma12_momenta.check",S12mom)
  end subroutine save_results











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


end program ed_hybrid_nonsu2



