MODULE COMMON
  USE SCIFOR
  USE EDIPACK
  implicit none


  integer :: Mnambu=1

  real(8),allocatable    :: Bath(:),Wlist(:)



  real(8),allocatable    :: evals(:)
  real(8),allocatable    :: dens(:)
  real(8),allocatable    :: docc(:)
  real(8),allocatable    :: phisc(:,:)
  real(8),allocatable    :: magX(:)
  real(8),allocatable    :: exciton(:) 
  real(8),allocatable    :: energy(:)
  real(8),allocatable    :: doubles(:)
  real(8),allocatable    :: imp(:)
  real(8),allocatable    :: Smom(:,:)
  real(8),allocatable    :: ASmom(:,:)
  real(8),allocatable    :: ASmomAB(:,:,:)
  real(8),allocatable    :: S11mom(:,:),S12mom(:,:)
  real(8),allocatable    :: SmomNN(:,:,:,:,:)
  complex(8),allocatable :: rdm(:,:)

  real(8),allocatable    :: evalsR(:)
  real(8),allocatable    :: densR(:)
  real(8),allocatable    :: doccR(:)
  real(8),allocatable    :: phiscR(:,:)
  real(8),allocatable    :: magXR(:)
  real(8),allocatable    :: excitonR(:)
  real(8),allocatable    :: energyR(:)
  real(8),allocatable    :: doublesR(:)
  real(8),allocatable    :: impR(:)
  real(8),allocatable    :: SmomR(:,:)
  real(8),allocatable    :: ASmomR(:,:)
  real(8),allocatable    :: ASmomABR(:,:,:)
  real(8),allocatable    :: S11momR(:,:),S12momR(:,:)
  real(8),allocatable    :: SmomNNR(:,:,:,:,:)
  complex(8),allocatable :: rdmR(:,:)
  
  integer                :: irank,comm,rank,size2,ierr
  logical                :: master








contains





  subroutine print_status()
    logical :: umatrix,hk
    umatrix= ED_READ_UMATRIX
    hk     = ED_USE_KANAMORI
    write(*,"(A50,A)")"ED_MODE   = "//str(ED_MODE)
    write(*,"(A50,A)")"BATH_TYPE = "//str(BATH_TYPE)
    write(*,"(A50,A)")"SPARSE_H  = "//str(ED_SPARSE_H)
    write(*,"(A50,A)")"U_MATRIX  = "//str(UMATRIX)
    write(*,"(A50,A)")"USE HK    = "//str(HK)
    call wait(1000)
  end subroutine print_status









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



  function mso2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin*Mnambu)stop "error so2j_index: ispin>Nspin*Mnambu"
    isporb=(ispin-1)*Nspin*Mnambu + iorb
  end function mso2j_index


  function mso2j(fg) result(g)
    complex(8),dimension(Nspin*Mnambu,Nspin*Mnambu,Norb,Norb) :: fg
    complex(8),dimension(Nspin*Norb*Mnambu,Nspin*Norb*Mnambu) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin*Mnambu
       do jspin=1,Nspin*Mnambu
          do iorb=1,Norb
             do jorb=1,Norb
                i=mso2j_index(ispin,iorb)
                j=mso2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function mso2j

  function j2mso(fg) result(g)
    complex(8),dimension(Nspin*Norb*Mnambu,Nspin*Norb*Mnambu) :: fg
    complex(8),dimension(Nspin*Mnambu,Nspin*Mnambu,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin*Mnambu
       do jspin=1,Nspin*Mnambu
          do iorb=1,Norb
             do jorb=1,Norb
                i=mso2j_index(ispin,iorb)
                j=mso2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2mso







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




END MODULE COMMON
