program ed_replica_superc
  USE EDIPACK
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  USE COMMON
  implicit none
  integer                                     :: i,js,Nso,Nson,Nsymm,Nmomenta
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,jspin
  complex(8),allocatable                      :: Smats(:,:,:,:,:,:)
  complex(8),allocatable                      :: Hloc(:,:,:,:)
  real(8)                                     :: Delta
  character(len=16)                           :: finput
  logical                                     :: dsave  
  complex(8),dimension(4,4)                   :: GammaN,GammaPhiAA,GammaPhiAB
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
  call parse_input_variable(delta,"DELTA",finput,default=0.d0)
  call parse_cmd_variable(dsave,"dsave",default=.false.)
  !
  !
  call ed_read_input(trim(finput))
  !
  if(bath_type/="replica")stop "Wrong setup from input file: non replica bath"
  if(ed_mode/="superc")stop "Wrong setup from input file: non superc ed_mode"
  if(Norb/=2)stop "Wrong setup from input file: Norb!=2"
  if(Nspin/=1 )stop "Wrong setup from input file: Nspin/=1"
  Mnambu=2
  Nso=Nspin*Norb
  Nson=Nso*Mnambu
  Nmomenta=4
  !
  !Allocate Weiss Field:
  allocate(Smats(2,Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(dens(Norb))
  allocate(docc(Norb))
  allocate(energy(4))
  allocate(doubles(4))
  allocate(phisc(Norb,Norb))
  allocate(imp(2))
  allocate(Smom(Norb,Nmomenta))
  allocate(ASmomAB(Norb,Norb,Nmomenta))
  !
  allocate(Wlist(Lmats))  
  Wlist = pi/beta*(2*arange(1,Lmats)-1)
  !
  ! Matrices for replica hamiltonian in Nambu representation
  gammaN    =kron( pauli_sigma_z, pauli_tau_0)
  gammaPhiAA=kron( pauli_sigma_x, pauli_tau_0 )
  gammaPhiAB=kron( pauli_sigma_x, pauli_tau_x )
  !
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  Hloc = zero
  if(Norb==2)then
     do js=1,Nspin
        Hloc(js,js,:,:)= Delta*pauli_sigma_z
     end do
  endif
  !
  ! Set up replica hamiltonian
  Nsymm=3
  allocate(lambdasym_vector(Nbath,Nsymm))
  allocate(Hsym_basis(Nspin*Mnambu,Nspin*Mnambu,Norb,Norb,Nsymm))
  !
  ! N
  Hsym_basis(:,:,:,:,1)=j2mso(GammaN(:Nson,:Nson))
  do i=1,Nbath
     lambdasym_vector(i,1) = -1.0 + 2.0*dble(i-1)/dble(Nbath-1)
  end do
  !
  ! PhiAA
  Hsym_basis(:,:,:,:,2)=j2mso(GammaPhiAA(:Nson,:Nson))
  lambdasym_vector(:,2)=0.1d0
  !
  ! PhiAB
  Hsym_basis(:,:,:,:,3)=j2mso(GammaPhiAB(:Nson,:Nson))
  lambdasym_vector(:,3)=0.2d0
  !
  !
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
    call ed_get_sigma(Smats(1,:,:,:,:,:),axis="m",type="n")
    call ed_get_sigma(Smats(2,:,:,:,:,:),axis="m",type="a")
    call ed_get_sigma(Smats(1,:,:,:,:,:),axis="m",type="n")
    call ed_get_sigma(Smats(2,:,:,:,:,:),axis="m",type="a")
    call ed_get_dens(dens)
    call ed_get_docc(docc)
    call ed_get_phi(phisc)
    call ed_get_eimp(energy)
    call ed_get_doubles(doubles)
    call ed_get_imp_info(imp)
    call ed_get_evals(evals)
    do i=1,Nmomenta
       do iorb=1,Norb
          call compute_momentum(Wlist,Smats(1,1,1,iorb,iorb,:),i,Smom(iorb,i))
          do jorb=1,Norb
             call compute_momentum(Wlist,Smats(2,1,1,iorb,jorb,:),i,ASmomAB(iorb,jorb,i))
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
    call read_array("phisc.check",phisc)
    call read_array("energy.check",energy)
    call read_array("doubles.check",doubles)
    call read_array("imp.check",imp)
    call read_array("Sigma_momenta.check",Smom)
    call read_array("Self_momenta.check",ASmomAB)
    if(allocated(densR))deallocate(densR)
    if(allocated(doccR))deallocate(doccR)
    if(allocated(phiscR))deallocate(phiscR)
    if(allocated(energyR))deallocate(energyR)
    if(allocated(doublesR))deallocate(doublesR)
    if(allocated(impR))deallocate(impR)
    if(allocated(evalsR))deallocate(evalsR)
    if(allocated(SmomR))deallocate(SmomR)
    if(allocated(ASmomABR))deallocate(ASmomABR)
    allocate(densR, source=dens)
    allocate(doccR, source=docc)
    allocate(phiscR, source=phisc)
    allocate(energyR, source=energy)
    allocate(doublesR, source=doubles)
    allocate(impR, source=imp)
    allocate(evalsR,source=evals)
    allocate(SmomR, source=Smom)
    allocate(ASmomABR, source=ASmomAB)
  end subroutine read_results



  subroutine test_results()
    call read_results()
    write(*,*)
    write(*,"(A50)") "Summary RESULTS:"
    call assert(dens,densR,"dens")
    call assert(docc,doccR,"docc")
    call assert(phisc,phiscR,"phisc")
    call assert(energy,energyR,"energy")
    call assert(doubles,doublesR,"doubles")
    call assert(imp,impR,"imp")
    call assert(evals,evalsR,"evals")
    call assert(Smom/SmomR,dble(ones(Norb,Nmomenta)),"Sigma_momenta 1:4",tol=1.0d-8)
    call assert(ASmomAB/ASmomABR,dble(ones(Norb,Norb,Nmomenta)),"Self_momenta 1:4",tol=1.0d-8)
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
    call save_array("phisc.check",phisc)
    call save_array("energy.check",energy)
    call save_array("doubles.check",doubles)
    call save_array("imp.check",imp)
    call save_array("Sigma_momenta.check",Smom)
    call save_array("Self_momenta.check",ASmomAB)
  end subroutine save_results




  subroutine set_twobody_hk()
    call ed_add_twobody_operator(1,"u",1,"d",1,"u",1,"d",-2.000000d0)
    call ed_add_twobody_operator(1,"d",1,"u",1,"d",1,"u",-2.000000d0)
    call ed_add_twobody_operator(2,"u",2,"d",2,"u",2,"d",-2.000000d0)
    call ed_add_twobody_operator(2,"d",2,"u",2,"d",2,"u",-2.000000d0)
    call ed_add_twobody_operator(1,"d",2,"u",1,"d",2,"u",-1.500000d0)
    call ed_add_twobody_operator(1,"u",2,"d",1,"u",2,"d",-1.500000d0)
    call ed_add_twobody_operator(2,"d",1,"u",2,"d",1,"u",-1.500000d0)
    call ed_add_twobody_operator(2,"u",1,"d",2,"u",1,"d",-1.500000d0)
    call ed_add_twobody_operator(1,"u",2,"u",1,"u",2,"u",-1.500000d0)
    call ed_add_twobody_operator(1,"d",2,"d",1,"d",2,"d",-1.500000d0)
    call ed_add_twobody_operator(2,"u",1,"u",2,"u",1,"u",-1.500000d0)
    call ed_add_twobody_operator(2,"d",1,"d",2,"d",1,"d",-1.500000d0)
    call ed_add_twobody_operator(1,"u",2,"u",2,"u",1,"u",0.250000d0)
    call ed_add_twobody_operator(1,"d",2,"d",2,"d",1,"d",0.250000d0)
    call ed_add_twobody_operator(2,"u",1,"u",1,"u",2,"u",0.250000d0)
    call ed_add_twobody_operator(2,"d",1,"d",1,"d",2,"d",0.250000d0)
    call ed_add_twobody_operator(1,"d",2,"u",2,"d",1,"u",0.250000d0)
    call ed_add_twobody_operator(1,"u",2,"d",2,"u",1,"d",0.250000d0)
    call ed_add_twobody_operator(2,"d",1,"u",1,"d",2,"u",0.250000d0)
    call ed_add_twobody_operator(2,"u",1,"d",1,"u",2,"d",0.250000d0)
    call ed_add_twobody_operator(1,"d",1,"u",2,"d",2,"u",0.250000d0)
    call ed_add_twobody_operator(1,"u",1,"d",2,"u",2,"d",0.250000d0)
    call ed_add_twobody_operator(2,"d",2,"u",1,"d",1,"u",0.250000d0)
    call ed_add_twobody_operator(2,"u",2,"d",1,"u",1,"d",0.250000d0)
  end subroutine set_twobody_hk




end program ed_replica_superc
