! STRUCTURE OF THE MODULE
!
! chi2_fit    : normal
! chi2_fit    : superc
!
! delta chi2  : normal
! delta chi2  : superc
! delta Dchi2 : normal
! delta       : normal
! delta       : superc
! Ddelta      : normal
!
! weiss chi2  : normal
! weiss chi2  : superc
! weiss Dchi2 : normal
! weiss       : normal
! weiss       : superc
! Dweiss      : normal
!
MODULE ED_FIT_REPLICA
!:synopsis: Routines for bath fitting, :code:`REPLICA` case
  USE ED_FIT_COMMON
  USE SF_SPIN, only: pauli_sigma_z
  USE SF_LINALG, only: kron, diag, trace

  implicit none
  private


  public :: chi2_fitgf_replica
  public :: chi2_fitgf_replica_superc


  complex(8),dimension(:,:,:,:,:),allocatable       :: FGmatrix
  complex(8),dimension(:,:,:,:,:),allocatable       :: FFmatrix

contains


  !+-------------------------------------------------------------+
  !PURPOSE  : Chi^2 interface for REPLICA BATH
  !+-------------------------------------------------------------+
  subroutine chi2_fitgf_replica(fg,bath_)
    complex(8),dimension(:,:,:,:,:)             :: fg ![Nspin][Nspin][Norb][Norb][Lmats]
    logical,dimension(Nspin,Nspin,Norb,Norb)    :: Hmask
    real(8),dimension(:),intent(inout)          :: bath_
    real(8),dimension(:),allocatable            :: array_bath
    integer                                     :: i,j,iorb,jorb,ispin,jspin,io,jo,ibath
    integer                                     :: iter,stride,counter,Asize
    real(8)                                     :: chi
    logical                                     :: check
    character(len=256)                          :: suffix
    integer                                     :: unit
    complex(8),dimension(:,:,:,:,:),allocatable :: fgand ![Nspin][][Norb][][Ldelta]  
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_replica: Fit"
#endif
    !
    call assert_shape(fg,[Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_replica","fg")
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_replica error: wrong bath dimensions"
    !
    if(cg_pow/=2.AND.cg_norm=="frobenius")then
       print *, "WARNING: CG_POW must be 2 for a meaningful definition of the Frobenius norm."
       print *, "         we'll still let you go ahead with the desired input, but please be "
       print *, "         be aware that CG_POW is not doing what you would expect for a chi^q"
    endif
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    allocate(array_bath(size(bath_)-1))
    Nsym  =bath_(1)
    array_bath=bath_(2:)
    !
    if(ed_all_g)then
       Hmask=.true.
       ! Aren't we sure about hermiticity?
       ! -> Hmask=Hreplica_mask(wdiag=.false.,uplo=.true.)
    else
       Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
    endif
    !
    totNso=count(Hmask)
    !
    allocate(getIspin(totNso),getJspin(totNso))
    allocate(getIorb(totNso) ,getJorb(totNso))
    counter=0
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                if (Hmask(ispin,jspin,iorb,jorb))then
                   counter=counter+1
                   getIspin(counter) = ispin
                   getIorb(counter)  = iorb
                   getJspin(counter) = jspin
                   getJorb(counter)  = jorb
                endif
             enddo
          enddo
       enddo
    enddo
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,5))Ldelta=size(fg,5)
    !
    allocate(FGmatrix(Nspin,Nspin,Norb,Norb,Ldelta))
    allocate(Gdelta(totNso,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    write(LOGfile,*)"Fitted functions",totNso
    do i=1,totNso
       Gdelta(i,1:Ldelta) = fg(getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
    enddo
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG chi2_fitgf_replica: cg_method:"//str(cg_method)//&
         ", cg_grad:"//str(cg_grad)//&
         ", cg_scheme:"//str(cg_scheme)
#endif
    !
    FGmatrix=fg
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default
       if(cg_grad==0)then
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_replica,grad_chi2_weiss_replica,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_replica,grad_chi2_delta_replica,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
          end select
       else
          select case (cg_scheme)
          case ("weiss")
             call fmin_cg(array_bath,chi2_weiss_replica,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case ("delta")
             call fmin_cg(array_bath,chi2_delta_replica,&
                  iter,chi,&
                  itmax=cg_niter,&
                  ftol=cg_Ftol,  &
                  istop=cg_stop, &
                  iverbose=(ed_verbose>3))
          case default
             stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
          end select
       endif
       !
       !
    case (1)
       if(cg_grad==0)then
          write(*,*) "                                                                                "
          write(*,*) "WARNING: analytic gradient not available with cg-method=1 (minimize f77 routine)"
          write(*,*) "         > we will force cg_grad=1 (so let the routine estimate the gradient)   "
          write(*,*) "                                                                                "
       endif
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,chi2_weiss_replica,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_replica,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_replica error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5,A)")"chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  All Orbs, All Spins"
    !
    suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    bath_(2:size(bath_))=array_bath
    call set_dmft_bath(bath_) ! *** array_bath --> dmft_bath *** (per write fit result)
    call write_dmft_bath(LOGfile)
    call save_dmft_bath()
    !
    allocate(fgand(Nspin,Nspin,Norb,Norb,Ldelta))
    if(cg_scheme=='weiss')then
       fgand = g0and_bath_function(xi*Xdelta(:))
    else
       fgand = delta_bath_function(xi*Xdelta(:))
    endif
    call write_fit_result()
    deallocate(fgand)
    !
    call get_dmft_bath(bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath()
    deallocate(FGmatrix,Gdelta,Xdelta,Wdelta)
    deallocate(getIspin,getJspin)
    deallocate(getIorb,getJorb)
    deallocate(array_bath)
    !
  contains
    !
    subroutine write_fit_result()
      integer   :: i,j,s,l,iorb,jorb,ispin,jspin
      !
      do l=1,totNso
         iorb = getIorb(l)
         jorb = getJorb(l)
         ispin = getIspin(l)
         jspin = getJspin(l)
         suffix="_l"//reg(txtfy(iorb))//&
              "_m"//reg(txtfy(jorb))//&
              "_s"//reg(txtfy(ispin))//&
              "_r"//reg(txtfy(jspin))//reg(ed_file_suffix)
         unit=free_unit()
         if(cg_scheme=='weiss')then
            open(unit,file="fit_weiss"//reg(suffix)//".ed")
         else
            open(unit,file="fit_delta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(unit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(ispin,jspin,iorb,jorb,i)),dimag(fgand(ispin,jspin,iorb,jorb,i)),&
                 dreal(fg(ispin,jspin,iorb,jorb,i)),dreal(fgand(ispin,jspin,iorb,jorb,i))
         enddo
         close(unit)
      enddo
    end subroutine write_fit_result
    !
  end subroutine chi2_fitgf_replica






  subroutine chi2_fitgf_replica_superc(fg,bath_)
    complex(8),dimension(:,:,:,:,:,:)                      :: fg ![2][Nspin][Nspin][Norb][Norb][Lmats]
    logical,dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: Hmask
    real(8),dimension(:),intent(inout)                     :: bath_
    real(8),dimension(:),allocatable                       :: array_bath
    integer                                                :: i,j,iorb,jorb,ispin,jspin,io,jo,ibath,inambu
    integer                                                :: iter,stride,counter,Asize
    real(8)                                                :: chi
    logical                                                :: check
    character(len=256)                                     :: suffix
    integer                                                :: unit
    complex(8),dimension(:,:,:,:,:,:),allocatable          :: fgand ![2][Nspin][][Norb][][Ldelta]  
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A)")"DEBUG chi2_fitgf_replica_superc: Fit"
#endif
    !
    call assert_shape(fg,[2,Nspin,Nspin,Norb,Norb,Lmats],"chi2_fitgf_replica_superc","fg")
    !
    check= check_bath_dimension(bath_)
    if(.not.check)stop "chi2_fitgf_replica_superc error: wrong bath dimensions"
    !
    call allocate_dmft_bath()
    call set_dmft_bath(bath_)
    allocate(array_bath(size(bath_)-1))
    Nsym  =bath_(1)
    array_bath=bath_(2:)
    !
    if(ed_all_g)then
       Hmask=.true.
    else
       Hmask=Hreplica_mask(wdiag=.true.,uplo=.false.)
    endif
    !
    ! Taking only upper-diagonal: superc => Nspin=1
    do iorb=1,Norb
       do jorb=1,Norb
          if(iorb>jorb)Hmask(1,1,iorb,jorb)=.false.
       enddo
    enddo
    !
    !Counting only the normal part as we get anomalous independently.
    totNso=count(Hmask(1,1,:,:))
    !
    !
    allocate(getIspin(totNso),getJspin(totNso))
    allocate(getIorb(totNso) ,getJorb(totNso))
    !
    counter=0
    !HERE NSPIN=1
    getIspin = 1
    getJspin = 1
    do iorb=1,Norb
       do jorb=1,Norb
          if (Hmask(1,1,iorb,jorb))then
             counter=counter+1
             getIorb(counter)   = iorb
             getJorb(counter)   = jorb
          endif
       enddo
    enddo
    !
    Ldelta = Lfit ; if(Ldelta>size(fg,6))Ldelta=size(fg,6)
    !
    allocate(FGmatrix(Nspin,Nspin,Norb,Norb,Ldelta))
    allocate(FFmatrix(Nspin,Nspin,Norb,Norb,Ldelta))
    allocate(Gdelta(totNso,Ldelta))
    allocate(Fdelta(totNso,Ldelta))
    allocate(Xdelta(Ldelta))
    allocate(Wdelta(Ldelta))
    !
    Xdelta = pi/beta*(2*arange(1,Ldelta)-1)
    !
    select case(cg_weight)
    case default
       Wdelta=1d0
    case(2)
       Wdelta=1d0*arange(1,Ldelta)
    case(3)
       Wdelta=Xdelta
    end select
    !
    write(LOGfile,*)"Fitted functions",totNso
    do i=1,totNso
       Gdelta(i,1:Ldelta) = fg(1,getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
       Fdelta(i,1:Ldelta) = fg(2,getIspin(i),getJspin(i),getIorb(i),getJorb(i),1:Ldelta)
    enddo
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A)")&
         "DEBUG chi2_fitgf_replica_superc: cg_method:"//str(cg_method)//&
         ", cg_grad:"//str(cg_grad)//&
         ", cg_scheme:"//str(cg_scheme)
#endif
    !
    FGmatrix=fg(1,:,:,:,:,:)
    FFmatrix=fg(2,:,:,:,:,:)
    !
    !
    if(cg_grad==0)then
       cg_grad=1
       write(LOGfile,*)"WARNING: REPLICA Superc does not support minimization with analytical gradient (cg_grad=0)."
       write(LOGfile,*)"         To spare you some time we fall back to numerical gradient: cg_grad=1"
       call sleep(2)
    endif
    !
    select case(cg_method)     !0=NR-CG[default]; 1=CG-MINIMIZE
    case default
       select case (cg_scheme)
       case ("weiss")
          call fmin_cg(array_bath,chi2_weiss_replica_superc,&
               iter,chi,&
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               istop=cg_stop, &
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cg(array_bath,chi2_delta_replica_superc,&
               iter,chi,&
               itmax=cg_niter,&
               ftol=cg_Ftol,  &
               istop=cg_stop, &
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_replica_superc error: cg_scheme != [weiss,delta]"
       end select
       !
       !
    case (1)
       select case (cg_scheme)
       case ("weiss")
          call fmin_cgminimize(array_bath,chi2_weiss_replica_superc,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case ("delta")
          call fmin_cgminimize(array_bath,chi2_delta_replica_superc,&
               iter,chi,itmax=cg_niter,ftol=cg_Ftol,&
               new_version=cg_minimize_ver,hh_par=cg_minimize_hh,&
               iverbose=(ed_verbose>3))
       case default
          stop "chi2_fitgf_replica_superc error: cg_scheme != [weiss,delta]"
       end select
       !
    end select
    !
    write(LOGfile,"(A,ES18.9,A,I5,A)")"chi^2|iter"//reg(ed_file_suffix)//'= ',chi," | ",iter,"  <--  All Orbs, All Spins"
    !
    suffix="_ALLorb_ALLspins"//reg(ed_file_suffix)
    unit=free_unit()
    open(unit,file="chi2fit_results"//reg(suffix)//".ed",position="append")
    write(unit,"(ES18.9,1x,I5)") chi,iter
    close(unit)
    !
    bath_(2:size(bath_))=array_bath
    call set_dmft_bath(bath_) ! *** array_bath --> dmft_bath *** (per write fit result)
    call write_dmft_bath(LOGfile)
    call save_dmft_bath()
    !
    allocate(fgand(2,Nspin,Nspin,Norb,Norb,Ldelta))
    if(cg_scheme=='weiss')then
       fgand(1,:,:,:,:,:) = g0and_bath_function(xi*Xdelta(:))
       fgand(2,:,:,:,:,:) = f0and_bath_function(xi*Xdelta(:))
    else
       fgand(1,:,:,:,:,:) = delta_bath_function(xi*Xdelta(:))
       fgand(2,:,:,:,:,:) =fdelta_bath_function(xi*Xdelta(:))
    endif
    call write_fit_result()
    deallocate(fgand)
    !
    call get_dmft_bath(bath_)                ! ***  dmft_bath --> bath_ ***    (bath in output)
    call deallocate_dmft_bath()
    deallocate(FGmatrix,FFmatrix,Gdelta,Fdelta,Xdelta,Wdelta)
    deallocate(getIspin,getJspin)
    deallocate(getIorb,getJorb)
    deallocate(array_bath)
    !
  contains
    !
    subroutine write_fit_result()
      integer            :: i,j,s,l,iorb,jorb,ispin,jspin,funit,gunit
      !
      do l=1,totNso
         iorb = getIorb(l)
         jorb = getJorb(l)
         ispin = getIspin(l)
         suffix=&
              "_l"//reg(txtfy(iorb))//&
              "_m"//reg(txtfy(jorb))//&
              reg(ed_file_suffix)
         if(cg_scheme=='weiss')then
            open(free_unit(gunit),file="fit_weiss"//reg(suffix)//".ed")
            open(free_unit(funit),file="fit_fweiss"//reg(suffix)//".ed")
         else
            open(free_unit(gunit),file="fit_delta"//reg(suffix)//".ed")
            open(free_unit(funit),file="fit_fdelta"//reg(suffix)//".ed")
         endif
         do i=1,Ldelta
            write(gunit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(1,ispin,ispin,iorb,jorb,i)),dimag(fgand(1,ispin,ispin,iorb,jorb,i)),&
                 dreal(fg(1,ispin,ispin,iorb,jorb,i)),dreal(fgand(1,ispin,ispin,iorb,jorb,i))
            write(funit,"(5F24.15)")Xdelta(i),&
                 dimag(fg(2,ispin,ispin,iorb,jorb,i)),dimag(fgand(2,ispin,ispin,iorb,jorb,i)),&
                 dreal(fg(2,ispin,ispin,iorb,jorb,i)),dreal(fgand(2,ispin,ispin,iorb,jorb,i))
         enddo
         close(gunit)
         close(funit)
      enddo
    end subroutine write_fit_result
    !
    !
  end subroutine chi2_fitgf_replica_superc








  !##################################################################
  !##################################################################
  !
  !            DELTA FIT
  !
  !##################################################################
  !##################################################################




  !######################################################################
  ! DELTA CHI2 FIT : NORMAL/NONSU2
  !######################################################################
  function chi2_delta_replica(a) result(chi2)
    real(8),dimension(:)                                         :: a
    real(8)                                                      :: chi2
    !
    select case(cg_norm)
    case ("elemental");chi2 = chi2_delta_replica_elemental(a)
    case ("frobenius");chi2 = chi2_delta_replica_frobenius(a)
    case default;stop "chi2_fitgf_replica error: cg_norm != [elemental,frobenius]"
    end select
    !
  end function chi2_delta_replica

  function chi2_delta_replica_elemental(a) result(chi2)
    real(8),dimension(:)                               :: a
    real(8)                                            :: chi2
    real(8),dimension(totNso)                          :: chi2_so
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    real(8),dimension(Ldelta)                          :: Ctmp
    integer                                            :: i,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_replica. a:",a
#endif
    !
    Delta = delta_replica(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ctmp =  abs(Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)
    chi2=chi2/Ldelta/totNso
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_replica. Chi**2:",chi2
#endif
    !
  end function chi2_delta_replica_elemental

  !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
  function chi2_delta_replica_frobenius(a) result(chi2)
    real(8),dimension(:)                                         :: a
    real(8)                                                      :: chi2
    real(8),dimension(Ldelta)                                    :: chi2_freq
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)           :: Delta
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                  :: Delta_so
    integer                                                      :: l
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_delta_replica_frobenius. a:",a
#endif
    !
    Delta = delta_replica(a)
    !
    do l=1,Ldelta
       Delta_so    =  nn2so_reshape(delta(:,:,:,:,l) - FGmatrix(:,:,:,:,l),Nspin,Norb)
       chi2_freq(l) =  sqrt(trace(matmul(Delta_so,conjg(transpose(Delta_so)))))
    enddo
    !
    chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
    chi2 = chi2/Ldelta/(Nspin*Norb) !Normalization over {iw} and Nlso
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_replica_frobenius. chi2:",chi2
#endif
    !
  end function chi2_delta_replica_frobenius




  !######################################################################
  ! DELTA CHI2 FIT : SUPERC
  !######################################################################
  function chi2_delta_replica_superc(a) result(chi2)
    !Evaluate the \chi^2 distance of \Delta_Anderson function.
    real(8),dimension(:) :: a
    real(8)              :: chi2
    !
    select case(cg_norm)
    case default;stop "chi2_fitgf_replica_normal error: cg_norm != [elemental,frobenius]"
    case ("elemental");chi2 = chi2_delta_replica_superc_elemental(a)
    case ("frobenius");chi2 = chi2_delta_replica_superc_frobenius(a)
    end select
  end function chi2_delta_replica_superc

  function chi2_delta_replica_superc_elemental(a) result(chi2)
    real(8),dimension(:)                                 :: a
    real(8)                                              :: chi2
    real(8),dimension(totNso)                            :: chi2_so
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    real(8),dimension(Ldelta)                            :: Ctmp,Ftmp
    integer                                              :: i,l,iorb,jorb,ispin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_replica_superc. a:",a
#endif
    !
    Delta = delta_replica_superc(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       !
       Ctmp =  abs(Gdelta(l,:)-Delta(1,ispin,ispin,iorb,jorb,:))
       Ftmp =  abs(Fdelta(l,:)-Delta(2,ispin,ispin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta ) + sum( Ftmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)
    chi2=chi2/Ldelta/totNso
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_delta_replica_superc. Chi**2:",chi2
#endif
    !
  end function chi2_delta_replica_superc_elemental

  !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
  function chi2_delta_replica_superc_frobenius(a) result(chi2)
    real(8),dimension(:)                                 :: a
    real(8)                                              :: chi2
    real(8),dimension(Ldelta)                            :: chi2_freq
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    complex(8),dimension(Nspin*Norb,Nspin*Norb)          :: D,F
    integer                                              :: l
    !
    Delta = delta_replica_superc(a)
    !
    do l=1,Ldelta
       D    =  nn2so_reshape(delta(1,:,:,:,:,l) - FGmatrix(:,:,:,:,l),Nspin,Norb)
       F    =  nn2so_reshape(delta(2,:,:,:,:,l) - FFmatrix(:,:,:,:,l),Nspin,Norb)          
       chi2_freq(l) =  &
            sqrt(trace(matmul(D,conjg(transpose(D)))))  + &
            sqrt(trace(matmul(F,conjg(transpose(F)))))
    enddo
    !
    chi2 = sum(chi2_freq**cg_pow/Wdelta)
    chi2 = chi2/Ldelta/(Nspin*Norb)
    !
  end function chi2_delta_replica_superc_frobenius






  !######################################################################
  ! DELTA GRAD chi2: NORMAL
  !######################################################################
  function grad_chi2_delta_replica(a) result(dchi2)
    real(8),dimension(:)                                         :: a
    real(8),dimension(size(a))                                   :: dchi2
    !
    select case(cg_norm)
    case ("elemental");dchi2 = grad_chi2_delta_replica_elemental(a)
    case ("frobenius");dchi2 = grad_chi2_delta_replica_frobenius(a)
    case default;stop "chi2_fitgf_replica error: cg_norm != [elemental,frobenius]"
    end select
    !
  end function grad_chi2_delta_replica

  !+---------------------------------------------------------------------+
  !PURPOSE: Evaluate the gradient \Grad\chi^2 of \Delta_Anderson function.
  !+---------------------------------------------------------------------+
  function grad_chi2_delta_replica_elemental(a) result(dchi2)
    real(8),dimension(:)                                       :: a
    real(8),dimension(size(a))                                 :: dchi2
    real(8),dimension(totNso,size(a))                          :: df
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: Delta
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    complex(8),dimension(Ldelta)                               :: Ftmp
    real(8),dimension(Ldelta)                                  :: Ctmp
    integer                                                    :: i,j,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_replica. a:",a
#endif
    !
    Delta  = delta_replica(a)
    dDelta = grad_delta_replica(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ftmp = Gdelta(l,:)-Delta(ispin,jspin,iorb,jorb,:)
       Ctmp = abs(Ftmp)**(cg_pow-2)
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Ftmp)*dreal(dDelta(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dDelta(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta )
       enddo
    enddo
    !
    dchi2 = -cg_pow*sum(df,1)/Ldelta/totNso     !sum over all orbital indices
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_replica. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_delta_replica_elemental


  !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
  function grad_chi2_delta_replica_frobenius(a) result(dchi2)
    real(8),dimension(:)                                       :: a
    real(8),dimension(size(a))                                 :: dchi2
    real(8),dimension(Ldelta,size(a))                          :: df
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: Delta
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    complex(8),dimension(Ldelta)                               :: Ftmp
    real(8),dimension(Ldelta,size(a))                          :: dChi_freq
    integer                                                    :: i,j,idelta,iorb,jorb,ispin,jspin
    !
    Delta  = delta_replica(a)
    dDelta = grad_delta_replica(a)
    Ftmp=zero
    df=zero
    !
    do idelta=1,Ldelta
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   !
                   Ftmp(idelta) = Ftmp(idelta) + abs(Delta(ispin,jspin,iorb,jorb,idelta)-FGmatrix(ispin,jspin,iorb,jorb,idelta))**2
                   do j=1,size(a)
                      df(idelta,j) = df(idelta,j) + &
                           dreal(Delta(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                           dreal(dDelta(ispin,jspin,iorb,jorb,idelta,j)) + &
                           dimag(Delta(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                           dimag(dDelta(ispin,jspin,iorb,jorb,idelta,j))
                   enddo
                enddo
             enddo
          enddo
       enddo
       Ftmp(idelta) = cg_pow * (sqrt(Ftmp(idelta))**(cg_pow-2)) / Wdelta(idelta)
       dchi_freq(idelta,:) = Ftmp(idelta) * df(idelta,:)
    enddo
    !
    dchi2 = sum(dchi_freq,1)/Ldelta/(Nspin*Norb)
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_delta_replica_frobenius. dChi2:",dchi2
#endif
    !
  end function grad_chi2_delta_replica_frobenius



  !######################################################################
  ! DELTA GRAD chi2: SUPERC ??
  !######################################################################





  !######################################################################
  ! DELTA + GRAD DELTA: NORMAL
  !######################################################################
  function delta_replica(a) result(Delta)
    real(8),dimension(:)                               :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: Delta
    integer                                            ::  ibath
    integer                                            :: i,stride
    complex(8),dimension(Nspin*Norb,Nspin*Norb)        :: Haux,Htmp
    complex(8),dimension(Nspin,Nspin,Norb,Norb)        :: invH_knn
    real(8)                                            :: V
    real(8),dimension(Nbath)                           :: Vk
    type(nsymm_vector),dimension(Nbath)                :: Lk
    !
    !Get Hs
    stride = 0
    do ibath=1,Nbath
       allocate(Lk(ibath)%element(Nsym))
       !
       !Get Vs
       stride           = stride + 1
       Vk(ibath)        = a(stride)
       !Get Lambdas
       Lk(ibath)%element= a(stride+1:stride+Nsym)
       stride           = stride + Nsym
    enddo
    !
    Delta=zero
    do ibath=1,Nbath
       V    = Vk(ibath)
       Htmp = nn2so_reshape(build_Hreplica(Lk(ibath)%element),Nspin,Norb)
       do i=1,Ldelta
          Haux     = zeye(Nspin*Norb)*xi*Xdelta(i) - Htmp
          call inv(Haux)
          Delta(:,:,:,:,i)=Delta(:,:,:,:,i) + V*so2nn_reshape(Haux,Nspin,Norb)*V
       enddo
       deallocate(Lk(ibath)%element)
    enddo
    !
  end function delta_replica


  function grad_delta_replica(a) result(dDelta)
    real(8),dimension(:)                                       :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    integer                                                    :: ispin,iorb,jorb,ibath
    integer                                                    :: i,k,l,stride
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                :: H_reconstructed,Htmp,Hbasis_so
    complex(8),dimension(Nspin*Norb,Nspin*Norb,Ldelta)         :: Haux
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: invH_knn
    real(8),dimension(Nbath)                                   :: Vk
    type(nsymm_vector),dimension(Nbath)                        :: Lk
    !
    !
    !Get Hs
    stride = 0
    do ibath=1,Nbath
       allocate(Lk(ibath)%element(Nsym))
       !
       !Get Vs
       stride           = stride + 1
       Vk(ibath)        = a(stride)
       !Get Lambdas
       Lk(ibath)%element= a(stride+1:stride+Nsym)
       stride           = stride + Nsym
    enddo
    !
    dDelta=zero
    stride=0
    do ibath=1,Nbath
       H_reconstructed= nn2so_reshape( build_Hreplica(Lk(ibath)%element) ,Nspin,Norb)
       do i=1,Ldelta
          Haux(:,:,i) = zeye(Nspin*Norb)*xi*Xdelta(i) - H_reconstructed
          call inv(Haux(:,:,i))
          invH_knn(:,:,:,:,i) = so2nn_reshape(Haux(:,:,i),Nspin,Norb)
       enddo
       !Derivate_Vp
       stride = stride + 1
       dDelta(:,:,:,:,:,stride)=2d0*Vk(ibath)*invH_knn(:,:,:,:,:)
       !
       !Derivate_lambda_p
       do k=1,Nsym
          stride = stride + 1
          Hbasis_so=nn2so_reshape(Hb%basis(k)%O,Nspin,Norb)
          do l=1,Ldelta
             Htmp = ((Haux(:,:,l) .x. Hbasis_so)) .x. Haux(:,:,l)
             dDelta(:,:,:,:,l,stride)=so2nn_reshape(Vk(ibath)**2*Htmp,Nspin,Norb)
          enddo
       enddo
       !
    enddo
  end function grad_delta_replica





  !######################################################################
  ! DELTA + GRAD DELTA: SUPERC
  !######################################################################
  function delta_replica_superc(a) result(Delta)
    real(8),dimension(:)                                      :: a
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta)      :: Delta
    integer                                                   :: ispin,jspin,iorb,jorb,ibath
    integer                                                   :: i,stride
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: Haux,Htmp
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: invH_knn
    real(8),dimension(Nbath)                                  :: Vk
    type(nsymm_vector),dimension(Nbath)                       :: Lk
    ! complex(8),dimension(Nnambu*Norb,Nnambu*Norb)             :: JJ = sig_z x 1
    !
    !Get Hs
    stride = 0
    do ibath=1,Nbath
       allocate(Lk(ibath)%element(Nsym))
       !
       !Get Vs
       stride           = stride + 1
       Vk(ibath)        = a(stride)
       !Get Lambdas
       Lk(ibath)%element= a(stride+1:stride+Nsym)
       stride           = stride + Nsym
    enddo
    !
    Delta=zero
    do ibath=1,Nbath
       Htmp = nn2so_reshape(build_Hreplica(Lk(ibath)%element),Nnambu*Nspin,Norb)
       do i=1,Ldelta
          Haux  = zeye(Nnambu*Nspin*Norb)*xi*Xdelta(i) - Htmp
          call inv(Haux)
          ! Haux =matmul(matmul(JJ,Haux),JJ) !?? !AA: Based on my calculations this is not needed
          invH_knn = so2nn_reshape(Haux,Nnambu*Nspin,Norb)
          Delta(1,1,1,:,:,i)=Delta(1,1,1,:,:,i) + Vk(ibath)*invH_knn(1,1,:,:)*Vk(ibath)
          Delta(2,1,1,:,:,i)=Delta(2,1,1,:,:,i) + Vk(ibath)*invH_knn(1,2,:,:)*Vk(ibath)
       enddo
       deallocate(Lk(ibath)%element)
    enddo
    !
  end function delta_replica_superc














  !##################################################################
  !##################################################################
  !
  !            WEISS FIT
  !
  !##################################################################
  !##################################################################


  !######################################################################
  ! WEISS CHI2 FIT : NORMAL
  !######################################################################
  function chi2_weiss_replica(a) result(chi2)
    real(8),dimension(:)                                         :: a
    real(8)                                                      :: chi2
    !
    select case(cg_norm)
    case ("elemental");chi2 = chi2_weiss_replica_elemental(a)
    case ("frobenius");chi2 = chi2_weiss_replica_frobenius(a)
    case default;stop "chi2_fitgf_replica error: cg_norm != [elemental,frobenius]"
    end select
    !
  end function chi2_weiss_replica

  function chi2_weiss_replica_elemental(a) result(chi2)
    real(8),dimension(:)                               :: a
    real(8)                                            :: chi2
    real(8),dimension(totNso)                          :: chi2_so
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
    real(8),dimension(Ldelta)                          :: Ctmp
    integer                                            :: i,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_replica. a:",a
#endif
    !
    g0and = g0and_replica(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ctmp = abs(Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)
    chi2=chi2/Ldelta/totNso
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_replica. Chi**2:",chi2
#endif
    !
  end function chi2_weiss_replica_elemental

  !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
  function chi2_weiss_replica_frobenius(a) result(chi2)
    real(8),dimension(:)                               :: a
    real(8)                                            :: chi2
    real(8),dimension(Ldelta)                          :: chi2_freq
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
    complex(8),dimension(Nspin*Norb,Nspin*Norb)        :: Delta_so
    integer                                            :: l
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_replica_frobenius. a:",a
#endif
    !
    g0and = g0and_replica(a)
    !
    do l=1,Ldelta
       Delta_so    =  nn2so_reshape(g0and(:,:,:,:,l) - FGmatrix(:,:,:,:,l),Nspin,Norb)
       chi2_freq(l) =  sqrt(trace(matmul(Delta_so,conjg(transpose(Delta_so)))))
    enddo
    !
    chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
    chi2 = chi2/Ldelta/(Nspin*Norb) !Normalization over {iw} and Nlso
    !
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_replica_frobenius. chi2:",chi2
#endif
    !
  end function chi2_weiss_replica_frobenius





  !######################################################################
  ! WEISS CHI2 FIT : SUPERC
  !######################################################################
  function chi2_weiss_replica_superc(a) result(chi2)
    !Evaluate the \chi^2 distance of \Weiss_Anderson function.
    real(8),dimension(:) :: a
    real(8)              :: chi2
    !
    select case(cg_norm)
    case default;stop "chi2_fitgf_replica_normal error: cg_norm != [elemental,frobenius]"
    case ("elemental");chi2 = chi2_weiss_replica_superc_elemental(a)
    case ("frobenius");chi2 = chi2_weiss_replica_superc_frobenius(a)
    end select
  end function chi2_weiss_replica_superc

  function chi2_weiss_replica_superc_elemental(a) result(chi2)
    real(8),dimension(:)                                 :: a
    real(8)                                              :: chi2
    real(8),dimension(totNso)                            :: chi2_so
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
    real(8),dimension(Ldelta)                            :: Ctmp,Ftmp
    integer                                              :: i,l,iorb,jorb,ispin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG chi2_weiss_replica_superc. a:",a
#endif
    !
    g0and = g0and_replica_superc(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       !
       Ctmp =  abs(Gdelta(l,:)-g0and(1,ispin,ispin,iorb,jorb,:))
       Ftmp =  abs(Fdelta(l,:)-g0and(2,ispin,ispin,iorb,jorb,:))
       chi2_so(l) = sum( Ctmp**cg_pow/Wdelta ) + sum( Ftmp**cg_pow/Wdelta )
    enddo
    !
    chi2=sum(chi2_so)
    chi2=chi2/Ldelta/totNso
#ifdef _DEBUG
    if(ed_verbose>3)write(Logfile,"(A,ES10.2)")"DEBUG chi2_weiss_replica_superc. Chi**2:",chi2
#endif
    !
  end function chi2_weiss_replica_superc_elemental

  function chi2_weiss_replica_superc_frobenius(a) result(chi2)
    real(8),dimension(:)                                 :: a
    real(8)                                              :: chi2
    real(8),dimension(Ldelta)                            :: chi2_freq
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta) :: g0and
    complex(8),dimension(Nspin*Norb,Nspin*Norb)          :: D,F
    integer                                              :: l
    !
    g0and = g0and_replica_superc(a)
    !
    do l=1,Ldelta
       D    =  nn2so_reshape(g0and(1,:,:,:,:,l) - FGmatrix(:,:,:,:,l),Nspin,Norb)
       F    =  nn2so_reshape(g0and(2,:,:,:,:,l) - FFmatrix(:,:,:,:,l),Nspin,Norb)          
       chi2_freq(l) =  &
            sqrt(trace(matmul(D,conjg(transpose(D)))))  + &
            sqrt(trace(matmul(F,conjg(transpose(F)))))
    enddo
    !
    chi2 = sum(chi2_freq**cg_pow/Wdelta) !Weighted sum over matsubara frqs
    chi2 = chi2/Ldelta/(Nspin*Norb)      !Normalization over {iw} and Nlso
    !
  end function chi2_weiss_replica_superc_frobenius






  !######################################################################
  ! WEISS GRAD chi2: NORMAL
  !######################################################################
  function grad_chi2_weiss_replica(a) result(dchi2)
    real(8),dimension(:)                                         :: a
    real(8),dimension(size(a))                                   :: dchi2
    !
    select case(cg_norm)
    case ("elemental");dchi2 = grad_chi2_weiss_replica_elemental(a)
    case ("frobenius");dchi2 = grad_chi2_weiss_replica_frobenius(a)
    case default;stop "chi2_fitgf_replica error: cg_norm != [elemental,frobenius]"
    end select
    !
  end function grad_chi2_weiss_replica

  function grad_chi2_weiss_replica_elemental(a) result(dchi2)
    real(8),dimension(:)                                       :: a
    real(8),dimension(size(a))                                 :: dchi2
    real(8),dimension(totNso,size(a))                          :: df
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: g0and
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dg0and
    complex(8),dimension(Ldelta)                               :: Ftmp
    real(8),dimension(Ldelta)                                  :: Ctmp
    integer                                                    :: i,j,l,iorb,jorb,ispin,jspin
    !
#ifdef _DEBUG
    if(ed_verbose>5)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_replica. a:",a
#endif
    !
    g0and  = g0and_replica(a)
    dg0and = grad_g0and_replica(a)
    !
    do l=1,totNso
       iorb = getIorb(l)
       jorb = getJorb(l)
       ispin = getIspin(l)
       jspin = getJspin(l)
       !
       Ftmp = Gdelta(l,:)-g0and(ispin,jspin,iorb,jorb,:)
       Ctmp = abs(Ftmp)**(cg_pow-2)
       do j=1,size(a)
          df(l,j)=&
               sum( dreal(Ftmp)*dreal(dg0and(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta ) + &
               sum( dimag(Ftmp)*dimag(dg0and(ispin,jspin,iorb,jorb,:,j))*Ctmp/Wdelta )
       enddo
    enddo
    !
    dchi2 = -cg_pow*sum(df,1)/Ldelta/totNso     !sum over all orbital indices
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_replica. dChi**2:",dchi2
#endif
    !
  end function grad_chi2_weiss_replica_elemental

  !> FROBENIUS NORM: global \chi^2 for all components, only i\omega are weighted
  function grad_chi2_weiss_replica_frobenius(a) result(dchi2)
    real(8),dimension(:)                                                 :: a
    real(8),dimension(size(a))                                           :: dchi2
    real(8),dimension(Ldelta,size(a))                                    :: df
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)                   :: g0and
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a))           :: dg0and
    complex(8),dimension(Ldelta)                                         :: Ftmp
    real(8),dimension(Ldelta,size(a))                                    :: dChi_freq
    integer                                                              :: i,j,idelta,iorb,jorb,ispin,jspin
    !
    !
    print*, "                                                                     "
    print*, "WARNING: the analytic gradient of the Weiss field Frobenius distance "
    print*, "         has been found giving /QUALITATIVELY WRONG/ fitted spectra! "
    print*, "         > the issue has emerged in some Nlat=Nspin=Norb=1 test runs "
    print*, "         > while the bug is investigated please switch cg_scheme to  "
    print*, "           'delta' or cg_norm to 'elemental' (or give up analytic cg)"
    print*, "                                                                     "
    !
    !
    g0and  = g0and_replica(a)
    dg0and = grad_g0and_replica(a)
    Ftmp=zero
    df=zero
    !
    do idelta=1,Ldelta
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   !
                   Ftmp(idelta) = Ftmp(idelta) + abs(g0and(ispin,jspin,iorb,jorb,idelta)-FGmatrix(ispin,jspin,iorb,jorb,idelta))**2
                   do j=1,size(a)
                      df(idelta,j) = df(idelta,j) + &
                           dreal(g0and(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                           dreal(dg0and(ispin,jspin,iorb,jorb,idelta,j)) + &
                           dimag(g0and(ispin,jspin,iorb,jorb,idelta) - FGmatrix(ispin,jspin,iorb,jorb,idelta)) * &
                           dimag(dg0and(ispin,jspin,iorb,jorb,idelta,j))
                   enddo
                enddo
             enddo
          enddo
       enddo
       Ftmp(idelta) = cg_pow * (sqrt(Ftmp(idelta))**(cg_pow-2)) / Wdelta(idelta)
       dchi_freq(idelta,:) = Ftmp(idelta) * df(idelta,:)
    enddo
    !
    dchi2 = sum(dchi_freq,1)/Ldelta/(Nspin*Norb)
    !
#ifdef _DEBUG
    if(ed_verbose>4)write(Logfile,"(A,"//str(size(a))//"ES10.2)")"DEBUG grad_chi2_weiss_replica_frobenius. dChi2:",dchi2
#endif
    !
  end function grad_chi2_weiss_replica_frobenius




  !######################################################################
  ! WEISS + GRAD WEISS: NORMAL
  !######################################################################
  function g0and_replica(a) result(G0and)
    real(8),dimension(:)                               :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
    complex(8),dimension(Norb*Nspin,Norb*Nspin)        :: FGorb
    integer                                            :: i,Nso
    !
    Nso   = Norb*Nspin
    Delta = delta_replica(a)
    do i=1,Ldelta
       FGorb = (xi*Xdelta(i)+xmu)*zeye(Nso) - nn2so_reshape(impHloc + Delta(:,:,:,:,i), Nspin,Norb)
       call inv(FGorb)
       G0and(:,:,:,:,i) = so2nn_reshape(FGorb,Nspin,Norb)
    enddo
  end function g0and_replica


  function grad_g0and_replica(a) result(dG0and)
    real(8),dimension(:)                                       :: a
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dG0and
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)         :: G0and
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta,size(a)) :: dDelta
    complex(8),dimension(Nspin*Norb,Nspin*Norb)                :: dDelta_so,dG0and_so,G0and_so
    integer                                                    :: ispin,iorb,jorb
    integer                                                    :: ik,l
    !
    G0and  = g0and_replica(a)
    dDelta = grad_delta_replica(a)
    !
    dG0and = zero
    do l=1,Ldelta
       G0and_so=nn2so_reshape(g0and(:,:,:,:,l),Nspin,Norb)
       do ik=1,size(a)
          dDelta_so=nn2so_reshape(dDelta(:,:,:,:,l,ik),Nspin,Norb)
          dG0and_so = (G0and_so .x. dDelta_so) .x. G0and_so
          dG0and(:,:,:,:,l,ik)=so2nn_reshape(dG0and_so,Nspin,Norb) !Check the sign, should it be +
       enddo
    enddo
    !
  end function grad_g0and_replica







  !######################################################################
  ! WEISS + GRAD WEISS: SUPERC
  !######################################################################
  function g0and_replica_superc(a) result(G0and)
    real(8),dimension(:)                                      :: a
    complex(8),dimension(2,Nspin,Nspin,Norb,Norb,Ldelta)      :: G0and,Delta
    complex(8),dimension(Nnambu*Norb*Nspin,Nnambu*Norb*Nspin) :: FGorb
    complex(8),dimension(Norb*Nspin,Norb*Nspin)               :: Z,zHloc
    complex(8)                                                :: iw
    integer                                                   :: i,Nson,Nso
    !
    Nso   = Norb*Nspin
    Nson  = Norb*Nspin*Nnambu
    Delta = delta_replica_superc(a) !only two components, 11,12
    !
    zHloc =  nn2so_reshape(impHloc,Nspin,Norb)
    do i=1,Ldelta
       iw = xi*Xdelta(i)
       Z = (iw+xmu)*zeye(Nso) - zHloc
       FGorb(      :Norb,      :Norb) = Z -       nn2so_reshape(Delta(1,:,:,:,:,i), Nspin,Norb)
       FGorb(      :Norb,Norb+1:    ) =   -       nn2so_reshape(Delta(2,:,:,:,:,i), Nspin,Norb)
       FGorb(Norb+1:    ,      :Norb) =   - conjg(nn2so_reshape(Delta(2,:,:,:,:,i), Nspin,Norb))
       FGorb(Norb+1:    ,Norb+1:    ) = Z + conjg(nn2so_reshape(Delta(1,:,:,:,:,i), Nspin,Norb))
       call inv(FGorb)
       G0and(1,:,:,:,:,i) = so2nn_reshape(FGorb(1:Norb,     1:Norb  ),Nspin,Norb)
       G0and(2,:,:,:,:,i) = so2nn_reshape(FGorb(1:Norb,Norb+1:2*Norb),Nspin,Norb)
    enddo
  end function g0and_replica_superc


END MODULE ED_FIT_REPLICA



















! JUNK CODE?
!
! function delta_replica(a) result(Delta)
!   real(8),dimension(:)                                :: a
!   complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta)  :: Delta
!   integer                                             :: ispin,jspin,iorb,jorb,ibath
!   integer                                             :: i,io,jo,ndx
!   complex(8),dimension(Nspin*Norb,Nspin*Norb)         :: V_k
!   complex(8),dimension(Nspin*Norb,Nspin*Norb,Ldelta)  :: invH_k
!   complex(8),dimension(Nspin,Nspin,Norb,Norb,Nbath)   :: invH_knn
!   type(effective_bath)                                :: dmft_bath_tmp
!   !
!   call allocate_dmft_bath(dmft_bath_tmp)
!   call init_dmft_bath_mask(dmft_bath_tmp)
!   call set_dmft_bath(a,dmft_bath_tmp)
!   !
!   Delta=zero
!   invH_k=zero
!   do i=1,Ldelta
!      invH_knn=zero
!      do ibath=1,Nbath
!         !
!         do ispin=1,Nspin
!            do jspin=1,Nspin
!               do iorb=1,Norb
!                  do jorb=1,Norb
!                     io = iorb + (ispin-1) * Norb
!                     jo = jorb + (jspin-1) * Norb
!                     invH_k(io,jo,i)=dmft_bath_tmp%h(ispin,jspin,iorb,jorb,ibath)
!                  enddo
!               enddo
!            enddo
!         enddo
!         !
!         invH_k(:,:,i) = zeye(Nspin*Norb) * xi * Xdelta(i) - invH_k(:,:,i)
!         call inv(invH_k(:,:,i))
!         invH_knn(:,:,:,:,ibath)=so2nn_reshape(invH_k(:,:,i),Nspin,Norb)
!         !
!      enddo
!      !
!      do ibath=1,Nbath
!         do ispin=1,Nspin
!            do jspin=1,Nspin
!               do iorb=1,Norb
!                  do jorb=1,Norb
!                     Delta(ispin,jspin,iorb,jorb,i)=Delta(ispin,jspin,iorb,jorb,i)+ &
!                          (dmft_bath_tmp%vr(ibath))*invH_knn(ispin,jspin,iorb,jorb,ibath)*dmft_bath_tmp%vr(ibath)
!                  enddo
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
!   !
!   call deallocate_dmft_bath(dmft_bath_tmp)
!   !
! end function delta_replica

! function g0and_replica(a) result(G0and)
!   real(8),dimension(:)                               :: a
!   complex(8),dimension(Nspin,Nspin,Norb,Norb,Ldelta) :: G0and,Delta
!   complex(8),dimension(Norb*Nspin,Norb*Nspin)        :: zeta,fgorb
!   integer                                            :: i,Nso
!   integer                                            :: iorb,jorb,ispin,jspin,io,jo
!   !
!   Nso = Norb*Nspin
!   !
!   Delta = delta_replica(a)
!   !
!   do i=1,Ldelta
!      fgorb=zero
!      zeta = (xi*Xdelta(i)+xmu)*zeye(Nso)
!      do ispin=1,Nspin
!         do jspin=1,Nspin
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  io = iorb + (ispin-1)*Norb
!                  jo = jorb + (jspin-1)*Norb
!                  fgorb(io,jo)   = zeta(io,jo) - impHloc(ispin,jspin,iorb,jorb) - Delta(ispin,jspin,iorb,jorb,i)
!               enddo
!            enddo
!         enddo
!      enddo
!      call inv(fgorb)
!      do ispin=1,Nspin
!         do jspin=1,Nspin
!            do iorb=1,Norb
!               do jorb=1,Norb
!                  io = iorb + (ispin-1)*Norb
!                  jo = jorb + (jspin-1)*Norb
!                  G0and(ispin,jspin,iorb,jorb,i) = fgorb(io,jo)
!               enddo
!            enddo
!         enddo
!      enddo
!   enddo
!   !
! end function g0and_replica






