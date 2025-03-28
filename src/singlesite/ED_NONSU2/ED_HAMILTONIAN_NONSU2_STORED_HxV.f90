! > BUILD STORED SPARSE HAMILTONIAN of the SECTOR
MODULE ED_HAMILTONIAN_NONSU2_STORED_HxV
!:synopsis: Routines for sparese matrix-vector product, :code:`NONSU2` case
  USE ED_HAMILTONIAN_NONSU2_COMMON
  implicit none
  private

  !>Sparse Matric constructors
  public :: ed_buildH_nonsu2_main


  !>Sparse Mat-Vec product using stored sparse matrix 
  public  :: spMatVec_nonsu2_main
#ifdef _MPI
  public  :: spMatVec_MPI_nonsu2_main
#endif





contains



  !####################################################################
  !             BUILD SPARSE HAMILTONIAN of the SECTOR
  !####################################################################
  subroutine ed_buildH_nonsu2_main(isector,Hmat)
    !
    !
    !
    ! Builds the sector Hamiltonian :math:`H` and save each term in a suitable sparse matrix instance for :f:var:`ed_total_ed` = :code:`True`. If the dimension :f:var:`dim` of the sector are smaller than :f:var:`lanc_dim_threshold` the global matrix is dumped to the optional variable :f:var:`hmat`.
    !
    ! All the different electronic terms  are collected in the same sparse matrix, possibly using rows splitting and local / non-local blocks according to the :f:var:`MPI_Allgatherv` algorithm: 
    !  * :math:`H_{\rm int} \rightarrow` :f:var:`sph0` : interaction part of the electronic Hamiltonian
    !  * :math:`H_{\rm imp} \rightarrow` :f:var:`sph0` : impurity part of the eletronic Hamiltonian 
    !  * :math:`H_{\rm bath} \rightarrow` :f:var:`sph0`: bath levels part of the eletronic Hamiltonian
    !  * :math:`H_{\rm hyb} \rightarrow` :f:var:`sph0` : impurity - bath coupling part of the eletronic Hamiltonian
    !
    integer                                           :: isector
    complex(8),dimension(:,:),optional                :: Hmat !optional dense matrix
    integer,dimension(Nlevels)                        :: ib
    integer,dimension(Ns)                             :: ibup,ibdw
    real(8),dimension(Norb)                           :: nup,ndw
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath) :: Hbath_tmp
    integer                                           :: first_state,last_state
    integer                                           :: first_state_up,last_state_up
    integer                                           :: first_state_dw,last_state_dw
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG ed_buildH_main NONSU2: build H"
#endif
    !
    if(.not.Hsector%status)stop "ed_buildh_main ERROR: Hsector NOT allocated"
    isector=Hsector%index
    !
    Dim = getdim(isector)
    !
    if(present(Hmat))call assert_shape(Hmat,[Dim,Dim],"ed_buildh_main","Hmat")
    !
    !
    !Get diagonal hybridization, bath energy
    if(allocated(diag_hybr))deallocate(diag_hybr)
    if(allocated(bath_diag))deallocate(bath_diag)
    select case (bath_type)
    case default
       Nfoo = size(dmft_bath%e,2)
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Nfoo,Nbath));bath_diag=0d0       
       do ibath=1,Nbath
          do ispin=1,Nspin             
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%v(ispin,iorb,ibath)
             enddo
             do iorb=1,Nfoo
                bath_diag(ispin,iorb,ibath)=dmft_bath%e(ispin,iorb,ibath)
             enddo
          enddo
       enddo
    case ("replica")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=0d0
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hreplica(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%v
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    case ("general")
       allocate(diag_hybr(Nspin,Norb,Nbath));diag_hybr=0d0
       allocate(bath_diag(Nspin,Norb,Nbath));bath_diag=0d0
       do ibath=1,Nbath
          Hbath_tmp(:,:,:,:,ibath) = build_Hgeneral(dmft_bath%item(ibath)%lambda)
          do ispin=1,Nspin
             do iorb=1,Norb
                diag_hybr(ispin,iorb,ibath)=dmft_bath%item(ibath)%vg(iorb+Norb*(ispin-1))
                bath_diag(ispin,iorb,ibath)=Hbath_tmp(ispin,ispin,iorb,iorb,ibath)
             enddo
          enddo
       enddo
    end select
    !
#ifdef _MPI
    if(MpiStatus)then
       call sp_set_mpi_matrix(MpiComm,spH0,mpiIstart,mpiIend,mpiIshift)
       call sp_init_matrix(MpiComm,spH0,Dim)
    else
       call sp_init_matrix(spH0,Dim)
    endif
#else
    call sp_init_matrix(spH0,Dim)
#endif
    !
    !-----------------------------------------------!
    !
    !IMPURITY  HAMILTONIAN
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NONSU2: stored/Himp"
#endif
    include "stored/Himp.f90"
    !
    !LOCAL INTERACTION
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NONSU2: stored/Hint"
#endif
    include "stored/Hint.f90"
    !
    !BATH HAMILTONIAN
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NONSU2: stored/Hbath"
#endif
    include "stored/Hbath.f90"
    !
    !IMPURITY- BATH HYBRIDIZATION
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")"DEBUG ed_buildH_NONSU2: stored/Himp_bath"
#endif
    include "stored/Himp_bath.f90"
    !
    !
    !-----------------------------------------------!
    !
    !
    if(present(Hmat))then
#ifdef _MPI
       if(MpiStatus)then
          call sp_dump_matrix(MpiComm,spH0,Hmat)
       else
          call sp_dump_matrix(spH0,Hmat)
       endif
#else
       call sp_dump_matrix(spH0,Hmat)
#endif          
    endif
    !
  end subroutine ed_buildH_nonsu2_main










  !####################################################################
  !        SPARSE MAT-VEC PRODUCT USING STORED SPARSE MATRIX 
  !####################################################################
  !+------------------------------------------------------------------+
  !PURPOSE: Perform the matrix-vector product H*v used in the
  ! - serial cmplx(H)*cmplx(V)
  ! - MPI cmplx(H)*cmplx(V)
  !+------------------------------------------------------------------+
  subroutine spMatVec_nonsu2_main(Nloc,v,Hv)
    !
    ! Serial version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in Arpack/Lanczos algorithm.
    ! This procedures applies one by one each term of the global Hamiltonian to an input vector using the stored sparse matrices.  
    !
    integer                         :: Nloc !Global dimension of the problem. :code:`size(v)=Nloc=size(Hv)`
    complex(8),dimension(Nloc)      :: v    !input vector (passed by Arpack/Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)      :: Hv   !output vector (required by Arpack/Lanczos) :math:`\vec{w}`
    integer                         :: i,j
    Hv=zero
    do i=1,Nloc
       matmul: do j=1,spH0%row(i)%Size
          Hv(i) = Hv(i) + spH0%row(i)%cvals(j)*v(spH0%row(i)%cols(j))
       end do matmul
    end do
  end subroutine spMatVec_nonsu2_main


#ifdef _MPI
  subroutine spMatVec_mpi_nonsu2_main(Nloc,v,Hv)
    !
    ! MPI parallel version of the matrix-vector product :math:`\vec{w}=H\times\vec{v}` used in P-Arpack/P-Lanczos algorithm.
    ! This procedures applies one by one each term of the global Hamiltonian to a part of the vector own by the thread using the stored sparse matrices.
    !
    integer                             :: Nloc !Local dimension of the vector chunk. :code:`size(v)=Nloc` with :math:`\sum_p` :f:var:`Nloc` = :f:var:`Dim`
    complex(8),dimension(Nloc)          :: v    !input vector part (passed by P-Arpack/P-Lanczos) :math:`\vec{v}`
    complex(8),dimension(Nloc)          :: Hv   !!output vector (required by P-Arpack/P-Lanczos) :math:`\vec{w}`
    integer                             :: i,j,mpiIerr
    integer                             :: N,MpiShift
    complex(8),dimension(:),allocatable :: vin
    integer,allocatable,dimension(:)    :: Counts,Offset
    !
    !
    if(MpiComm==MPI_UNDEFINED)stop "spHtimesV_mpi_cc ERRROR: MpiComm = MPI_UNDEFINED"
    if(.not.MpiStatus)stop "spMatVec_mpi_cc ERROR: MpiStatus = F"
    !
    MpiRank = get_Rank_MPI(MpiComm)
    MpiSize = get_Size_MPI(MpiComm)
    !
    N = 0
    call AllReduce_MPI(MpiComm,Nloc,N)
    !
    !Evaluate the local contribution: Hv_loc = Hloc*v
    MpiShift = spH0%Ishift
    Hv=0d0
    do i=1,Nloc
       local: do j=1,spH0%loc(i)%Size
          Hv(i) = Hv(i) + spH0%loc(i)%cvals(j)*v(spH0%loc(i)%cols(j)-MpiShift)
       end do local
    end do
    !
    allocate(Counts(0:MpiSize-1)) ; Counts(0:)=0
    allocate(Offset(0:MpiSize-1)) ; Offset(0:)=0
    !
    Counts(0:)        = N/MpiSize
    Counts(MpiSize-1) = N/MpiSize+mod(N,MpiSize)
    !
    do i=1,MpiSize-1
       Offset(i) = Counts(i-1) + Offset(i-1)
    enddo
    !
    allocate(vin(N)) ; vin = zero
    call MPI_Allgatherv(&
         v(1:Nloc),Nloc,MPI_Double_Complex,&
         vin      ,Counts,Offset,MPI_Double_Complex,&
         MpiComm,MpiIerr)
    !
    do i=1,Nloc                 !==spH0%Nrow
       matmul: do j=1,spH0%row(i)%Size
          Hv(i) = Hv(i) + spH0%row(i)%cvals(j)*vin(spH0%row(i)%cols(j))
       end do matmul
    end do
    !
  end subroutine spMatVec_mpi_nonsu2_main
#endif











end MODULE ED_HAMILTONIAN_NONSU2_STORED_HXV
