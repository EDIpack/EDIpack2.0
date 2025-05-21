program ed_normal_normal_ineq
  USE EDIPACK
  USE EDIPACK2INEQ
  USE SCIFOR
  USE MPI
  USE SF_MPI
  USE ASSERTING
  implicit none
  integer                            :: i,iw,jo,js,Nso,Nmomenta,Nlat,Nlso
  !Bath:
  integer                            :: Nb,iorb,jorb,ispin,jspin,inso,ilat
  real(8),allocatable                :: Bath(:,:),Wlist(:)
  !GFs and Sigma:
  complex(8),allocatable             :: Smats(:,:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable             :: Hloc(:,:,:,:,:)
  !variables for the model:
  real(8)                            :: Delta
  character(len=16)                  :: finput
  real(8),allocatable                :: evals(:),evalsR(:)
  real(8),allocatable,dimension(:,:) :: dens,docc,energy,doubles
  real(8),allocatable,dimension(:,:) :: densR,doccR,energyR,doublesR
  real(8),allocatable,dimension(:,:) :: Smom
  real(8),allocatable,dimension(:,:) :: SmomR
  !
  !MPI Vars:
  integer                            :: irank,comm,rank,size2,ierr
  logical                            :: master
  logical                            :: dsave,mpi_lanc
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
  if(bath_type/="normal")stop "Wrong setup from input file: non normal bath"
  if(ed_mode/="normal")stop "Wrong setup from input file: non normal ed_mode"
  if(Norb/=1)stop "Wrong setup from input file: Norb!=1"
  if(Nspin/=2 )stop "Wrong setup from input file: Nspin/=2"
  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nso
  Nmomenta=4
  !
  !Allocate Weiss Field:
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  !
  allocate(dens(Nlat,Norb))
  allocate(docc(Nlat,Norb))
  allocate(energy(Nlat,4))
  allocate(doubles(Nlat,4))
  allocate(Smom(Nlat,Nmomenta))
  !
  allocate(Wlist(Lmats))  
  Wlist = pi/beta*(2*arange(1,Lmats)-1)
  !
  !
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  Hloc = lso2nnn_reshape(kron(pauli_sigma_z,pauli_tau_0),Nlat,Nspin,Norb)
  !
  !
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nlat,Nb))
  !
  mpi_lanc=.false.
  call run_test(sparse=.true.,umatrix=.false.)  
  call run_test(sparse=.false.,umatrix=.false.)
  call run_test(sparse=.true.,umatrix=.true.)  
  call run_test(sparse=.false.,umatrix=.true.)
  !
  mpi_lanc=.true.
  call run_test(sparse=.true.,umatrix=.false.)  
  call run_test(sparse=.false.,umatrix=.false.)
  call run_test(sparse=.true.,umatrix=.true.)  
  call run_test(sparse=.false.,umatrix=.true.)
  !
  call finalize_MPI()

contains


  subroutine run_test(sparse,umatrix)
    logical :: sparse,umatrix
    integer :: ip
    ED_SPARSE_H    =sparse
    ED_READ_UMATRIX=umatrix
    if(umatrix)ED_USE_KANAMORI=.false.
    call ed_init_solver(bath)
    call ed_break_symmetry_bath(bath,sb_field, (/( (-1d0)**(ip+1), ip=1,Nlat)/) )
    call ed_set_Hloc(hloc,Nlat)
    if(master)then
       write(*,*) ""
       write(*,*) "ED_MODE = NORMAL   |   BATH_TYPE = NORMAL"
       write(*,*) "SPARSE_H= "//str(sparse)
       write(*,*) "U_MATRIX= "//str(umatrix)
       write(*,*) "MPI_LANC= "//str(mpi_lanc)
    endif
    call ed_solve(bath,mpi_lanc)
    call ed_get_sigma(Smats,Nlat,axis="m",type="n")
    call ed_get_dens(dens,Nlat)
    call ed_get_docc(docc,Nlat)
    call ed_get_eimp(energy,Nlat)
    call ed_get_doubles(doubles,Nlat)
    ! call ed_get_evals(evals)
    do i=1,Nmomenta
       do ilat=1,Nlat
          call compute_momentum(Wlist,Smats(ilat,1,1,1,1,:),i,Smom(ilat,i))
       enddo
    enddo
    call ed_finalize_solver()
    if(dsave)then
       write(*,*)"Saving results to .check files and exit"
       call save_results()
       stop
    endif
    if(master)then
       write(*,*) "TEST RESULTS"
       write(*,*) "ED_MODE = NORMAL   |   BATH_TYPE = NORMAL"
       write(*,*) "SPARSE_H= "//str(sparse)
       write(*,*) "U_MATRIX= "//str(umatrix)
       write(*,*) "MPI_LANC= "//str(mpi_lanc)
    endif
    call test_results(sparse,umatrix)
  end subroutine run_test


  subroutine read_results()
    call read_array("dens.check",dens)
    call read_array("docc.check",docc)
    call read_array("energy.check",energy)
    call read_array("doubles.check",doubles)
    call read_array("Sigma_momenta.check",Smom)
    if(allocated(densR))deallocate(densR)
    if(allocated(doccR))deallocate(doccR)
    if(allocated(energyR))deallocate(energyR)
    if(allocated(doublesR))deallocate(doublesR)
    if(allocated(SmomR))deallocate(SmomR)
    allocate(densR, source=dens)
    allocate(doccR, source=docc)
    allocate(energyR, source=energy)
    allocate(doublesR, source=doubles)
    allocate(SmomR, source=Smom)
  end subroutine read_results


  subroutine test_results(sparse,umatrix)
    logical  :: sparse,umatrix
    if(master)then
       write(*,*)""
       write(*,*) "Check RESULTS sparse_H, read_umatrix="//str(sparse)//","//str(umatrix)
    endif
    call read_results()
    call assert(dens,densR,"dens")
    call assert(docc,doccR,"docc")
    call assert(energy,energyR,"energy")
    call assert(doubles,doublesR,"doubles")
    call assert(Smom/SmomR,dble(ones(Nlat,Nmomenta)),"Sigma_momenta 1:4",tol=1.0d-8)
    if(master)then
       write(*,*)""
       write(*,*)""
       write(*,*)""
       write(*,*)""
    endif
    call wait(500)
  end subroutine test_results



  subroutine save_results()
    call save_array("dens.check",dens)
    call save_array("docc.check",docc)
    call save_array("energy.check",energy)
    call save_array("doubles.check",doubles)
    call save_array("Sigma_momenta.check",Smom)
  end subroutine save_results





  !--------------------------------------------------------------------!
  !Reshaping functions:                                                !
  !--------------------------------------------------------------------!

  function nnn2lso_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fin
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
       Fout(is,js) = Fin(ilat,ispin,jspin,iorb,jorb)
    enddo
  end function nnn2lso_reshape

  function lso2nnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fin
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do concurrent(ilat=1:Nlat,ispin=1:Nspin,jspin=1:Nspin,iorb=1:Norb,jorb=1:Norb)
       is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
       js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
       Fout(ilat,ispin,jspin,iorb,jorb) = Fin(is,js)
    enddo
  end function lso2nnn_reshape



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


end program ed_normal_normal_ineq



