program ed_replica_normal
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
  real(8)                :: Delta
  character(len=16)      :: finput
  real(8),allocatable    :: evals(:),evalsR(:)
  real(8),allocatable    :: dens(:),docc(:),exciton(:),energy(:),doubles(:),imp(:),Smom(:,:)
  real(8),allocatable    :: densR(:),doccR(:),excitonR(:),energyR(:),doublesR(:),impR(:),SmomR(:,:)
  !
  complex(8),dimension(4,4)                   :: GammaN,GammaE0
  real(8),dimension(:,:),allocatable          :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable :: Hsym_basis
  !MPI Vars:
  integer                                     :: irank,comm,rank,size2,ierr
  logical                                     :: master
  logical                                     :: dsave
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
  call parse_cmd_variable(dsave,"dsave",default=.false.)
  !
  call ed_read_input(trim(finput))
  !
  if(bath_type/="replica")stop "Wrong setup from input file: non replica bath"
  if(ed_mode/="normal")stop "Wrong setup from input file: not normal ed_mode"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=1 )stop "Wrong setup from input file: Nspin/=1"
  Nso=Nspin*Norb
  Nmomenta=4
  !
  !Allocate Weiss Field:
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))

  ! Matrices for replica hamiltonian
  gammaN =kron( pauli_sigma_0, pauli_tau_0)
  gammaE0=kron( pauli_sigma_0, pauli_tau_x )

  allocate(dens(Norb))
  allocate(docc(Norb))
  allocate(exciton(2))
  allocate(energy(4))
  allocate(doubles(4))
  allocate(imp(2))
  allocate(Smom(Norb,Nmomenta))


  allocate(Wlist(Lmats))  
  Wlist = pi/beta*(2*arange(1,Lmats)-1)
  !
  !

  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc = zero
  if(Norb==2)then
     do js=1,Nspin
        Hloc(js,js,:,:)= Delta*pauli_sigma_z
     end do
  endif
  !
  !
  ! Set up replica hamiltonian
  Nsymm=2
  allocate(lambdasym_vector(Nbath,Nsymm))
  allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,Nsymm))
  !
  ! N
  Hsym_basis(:,:,:,:,1)=j2so(GammaN(:Nso,:Nso))
  do i=1,Nbath
     lambdasym_vector(i,1) = -1.0 + 2.0*dble(i-1)/dble(Nbath-1)
  end do
  !
  ! E0
  Hsym_basis(:,:,:,:,2)=j2so(GammaE0(:Nso,:Nso))
  lambdasym_vector(:,2)=0.1d0
  !
  call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
  Nb=ed_get_bath_dimension(Nsymm)
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
    call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
    call ed_init_solver(bath)
    call ed_set_Hloc(hloc)
    write(*,*) ""
    write(*,*) "ED_MODE = NORMAL   |   BATH_TYPE = REPLICA"
    write(*,*) "SPARSE_H= "//str(sparse)
    write(*,*) "U_MATRIX= "//str(umatrix)
    call ed_solve(bath)
    call ed_get_sigma(Smats,axis="m",type="n")
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    call ed_get_exct(exciton(1),'0')
    call ed_get_exct(exciton(2),'z')
    call ed_get_eimp(energy)
    call ed_get_doubles(doubles)
    call ed_get_imp_info(imp)
    call ed_get_evals(evals)
    do i=1,Nmomenta
       do iorb=1,Norb
          call compute_momentum(Wlist,Smats(1,1,iorb,iorb,:),i,Smom(iorb,i))
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
    call read_array("exciton.check",exciton)
    call read_array("energy.check",energy)
    call read_array("doubles.check",doubles)
    call read_array("imp.check",imp)
    call read_array("Sigma_momenta.check",Smom)
    if(allocated(densR))deallocate(densR)
    if(allocated(doccR))deallocate(doccR)
    if(allocated(excitonR))deallocate(excitonR)
    if(allocated(energyR))deallocate(energyR)
    if(allocated(doublesR))deallocate(doublesR)
    if(allocated(impR))deallocate(impR)
    if(allocated(evalsR))deallocate(evalsR)
    if(allocated(SmomR))deallocate(SmomR)
    allocate(densR, source=dens)
    allocate(doccR, source=docc)
    allocate(excitonR, source=exciton)
    allocate(energyR, source=energy)
    allocate(doublesR, source=doubles)
    allocate(impR, source=imp)
    allocate(evalsR,source=evals)
    allocate(SmomR, source=Smom)
  end subroutine read_results




  subroutine test_results(sparse,umatrix)
    logical  :: sparse,umatrix
    write(*,*)""
    write(*,*) "Check RESULTS sparse_H, read_umatrix="//str(sparse)//","//str(umatrix)
    call read_results()
    call assert(dens,densR,"dens")
    call assert(docc,doccR,"docc")
    call assert(exciton,excitonR,"docc")
    call assert(energy,energyR,"energy")
    call assert(doubles,doublesR,"doubles")
    call assert(imp,impR,"imp")
    call assert(evals,evalsR,"evals")
    call assert(Smom/SmomR,dble(ones(Norb,Nmomenta)),"Sigma_momenta 1:4",tol=1.0d-8)
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
    call save_array("exciton.check",exciton)
    call save_array("energy.check",energy)
    call save_array("doubles.check",doubles)
    call save_array("imp.check",imp)
    call save_array("Sigma_momenta.check",Smom)
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


end program ed_replica_normal



