program ed_replica_nonsu2
  USE EDIPACK
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  USE COMMON
  implicit none
  integer                                     :: i,js,Nso,Nsymm,Nmomenta
  integer                                     :: Nb,iorb,jorb,ispin,jspin
  complex(8),allocatable                      :: Smats(:,:,:,:,:)
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8)                                     :: mh,lambda
  character(len=16)                           :: finput
  logical                                     :: dsave
  complex(8),dimension(4,4)                   :: Gamma1,Gamma2,Gamma5,GammaN,GammaS
  complex(8),dimension(4,4)                   :: GammaE0,GammaEx,GammaEy,GammaEz
  real(8),dimension(:,:),allocatable          :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable :: Hsym_basis

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
  if(ed_mode/='nonsu2')stop "Wrong setup from input file: ed_mode != nonsu2"
  if(bath_type/="replica")stop "Wrong setup from input file: non replica bath"
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

  !Allocate Weiss Field:
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(dens(Norb))
  allocate(docc(Norb))
  allocate(exciton(4))
  allocate(energy(4))
  allocate(doubles(4))
  allocate(imp(2))
  allocate(SmomNN(Nspin,Nspin,Norb,Norb,Nmomenta))
  !
  allocate(Wlist(Lmats))  
  Wlist = pi/beta*(2*arange(1,Lmats)-1)

  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc = j2so(Mh*Gamma5)
  !
  ! Set up replica hamiltonian
  Nsymm=4
  allocate(lambdasym_vector(Nbath,Nsymm))
  allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,Nsymm))
  Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(:,1)=Mh
  Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(:,2)=sb_field
  Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(:,3)=sb_field
  Hsym_basis(:,:,:,:,4)=j2so(GammaEx) ;lambdasym_vector(:,4)=-sb_field
  call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
  !
  Nb=ed_get_bath_dimension(Nsymm)
  allocate(Bath(Nb))
  !
  call run_test(sparse=.true.,umatrix=.false.,hk=.true.)  
  call run_test(sparse=.false.,umatrix=.false.,hk=.true.)
  call run_test(sparse=.true.,umatrix=.true.,hk=.false.)  
  call run_test(sparse=.false.,umatrix=.true.,hk=.false.)
  call run_test(sparse=.true.,umatrix=.false.,hk=.false.)  
  call run_test(sparse=.false.,umatrix=.false.,hk=.false.)
  !
  call finalize_MPI()


contains

  subroutine run_test(sparse,umatrix,hk)
    logical :: sparse,umatrix,hk
    ED_SPARSE_H    =sparse
    ED_READ_UMATRIX=umatrix
    ED_USE_KANAMORI=hk
    if(.not.umatrix.AND..not.hk)call set_twobody_hk()
    call print_status()
    call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
    call ed_init_solver(bath)
    call ed_set_Hloc(hloc)
    call ed_solve(bath)
    call ed_get_sigma(Smats,axis="m",type="n")
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    call ed_get_exct(exciton)
    call ed_get_eimp(energy)
    call ed_get_doubles(doubles)
    call ed_get_imp_info(imp)
    call ed_get_evals(evals)
    do i=1,Nmomenta
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   call compute_momentum(Wlist,Smats(ispin,jspin,iorb,jorb,:),i,SmomNN(ispin,jspin,iorb,jorb,i))
                enddo
             enddo
          enddo
       enddo
    enddo
    call ed_finalize_solver()
    if(dsave)then
       write(*,*)"Saving results to .check files and exit"
       call save_results()
       stop
    endif
    call test_results()
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
    call read_array("Sigma_momenta.check",SmomNN)
    if(allocated(densR))deallocate(densR)
    if(allocated(doccR))deallocate(doccR)
    if(allocated(excitonR))deallocate(excitonR)
    if(allocated(energyR))deallocate(energyR)
    if(allocated(doublesR))deallocate(doublesR)
    if(allocated(impR))deallocate(impR)
    if(allocated(evalsR))deallocate(evalsR)
    if(allocated(SmomNNR))deallocate(SmomNNR)
    allocate(densR, source=dens)
    allocate(doccR, source=docc)
    allocate(excitonR, source=exciton)
    allocate(energyR, source=energy)
    allocate(doublesR, source=doubles)
    allocate(impR, source=imp)
    allocate(evalsR,source=evals)
    allocate(SmomNNR, source=SmomNN)
  end subroutine read_results



  subroutine test_results()
    call read_results()
    write(*,*)
    write(*,"(A50)") "Summary RESULTS:"
    call assert(dens,densR,"dens")
    call assert(docc,doccR,"docc")
    call assert(exciton,excitonR,"exciton")
    call assert(energy,energyR,"energy")
    call assert(doubles,doublesR,"doubles")
    call assert(imp,impR,"imp")
    call assert(evals,evalsR,"evals")
    call assert(SmomNN/SmomNNR,dble(ones(Nspin,Nspin,Norb,Norb,Nmomenta)),"Sigma_momenta 1:4",tol=1.0d-8)
    call print_status()
    write(*,*)""
    write(*,*)""
    write(*,*)""
    call wait(1000)
  end subroutine test_results


  subroutine save_results()
    call save_array("evals.check",evals)
    call save_array("dens.check",dens)
    call save_array("docc.check",docc)
    call save_array("exciton.check",exciton)
    call save_array("energy.check",energy)
    call save_array("doubles.check",doubles)
    call save_array("imp.check",imp)
    call save_array("Sigma_momenta.check",SmomNN)
  end subroutine save_results




  subroutine set_twobody_hk()
    call ed_add_twobody_operator(1,"u",1,"d",1,"u",1,"d",1.000000d0)
    call ed_add_twobody_operator(1,"d",1,"u",1,"d",1,"u",1.000000d0)
    call ed_add_twobody_operator(2,"u",2,"d",2,"u",2,"d",1.000000d0)
    call ed_add_twobody_operator(2,"d",2,"u",2,"d",2,"u",1.000000d0)
    call ed_add_twobody_operator(1,"d",2,"u",1,"d",2,"u",1.000000d0)
    call ed_add_twobody_operator(1,"u",2,"d",1,"u",2,"d",1.000000d0)
    call ed_add_twobody_operator(2,"d",1,"u",2,"d",1,"u",1.000000d0)
    call ed_add_twobody_operator(2,"u",1,"d",2,"u",1,"d",1.000000d0)
    call ed_add_twobody_operator(1,"u",2,"u",1,"u",2,"u",1.000000d0)
    call ed_add_twobody_operator(1,"d",2,"d",1,"d",2,"d",1.000000d0)
    call ed_add_twobody_operator(2,"u",1,"u",2,"u",1,"u",1.000000d0)
    call ed_add_twobody_operator(2,"d",1,"d",2,"d",1,"d",1.000000d0)
    call ed_add_twobody_operator(1,"u",2,"u",2,"u",1,"u",0.010000d0)
    call ed_add_twobody_operator(1,"d",2,"d",2,"d",1,"d",0.010000d0)
    call ed_add_twobody_operator(2,"u",1,"u",1,"u",2,"u",0.010000d0)
    call ed_add_twobody_operator(2,"d",1,"d",1,"d",2,"d",0.010000d0)
    call ed_add_twobody_operator(1,"d",2,"u",2,"d",1,"u",0.010000d0)
    call ed_add_twobody_operator(1,"u",2,"d",2,"u",1,"d",0.010000d0)
    call ed_add_twobody_operator(2,"d",1,"u",1,"d",2,"u",0.010000d0)
    call ed_add_twobody_operator(2,"u",1,"d",1,"u",2,"d",0.010000d0)
    call ed_add_twobody_operator(1,"d",1,"u",2,"d",2,"u",0.010000d0)
    call ed_add_twobody_operator(1,"u",1,"d",2,"u",2,"d",0.010000d0)
    call ed_add_twobody_operator(2,"d",2,"u",1,"d",1,"u",0.010000d0)
    call ed_add_twobody_operator(2,"u",2,"d",1,"u",1,"d",0.010000d0)
  end subroutine set_twobody_hk


end program ed_replica_nonsu2
