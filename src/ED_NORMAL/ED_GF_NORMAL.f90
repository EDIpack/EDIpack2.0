MODULE ED_GF_NORMAL
  USE SF_CONSTANTS, only:one,xi,zero,pi
  USE SF_TIMER  
  USE SF_IOTOOLS, only: str,reg,txtfy
  USE SF_LINALG,  only: inv,eigh,eye
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  USE ED_EIGENSPACE
  USE ED_BATH
  USE ED_BATH_FUNCTIONS
  USE ED_SETUP
  USE ED_SECTOR
  USE ED_HAMILTONIAN_NORMAL
  implicit none
  private


  public :: build_gf_normal
  public :: build_sigma_normal


  integer                               :: istate
  integer                               :: isector,jsector
  integer                               :: idim,idimUP,idimDW
  !
  integer                               :: ialfa,jalfa
  integer                               :: ipos,jpos
  integer                               :: i,j,m
  integer                               :: iph,i_el
  real(8)                               :: sgn,norm2
  !
  real(8),allocatable                   :: vvinit(:)
  real(8),allocatable                   :: alfa_(:),beta_(:)
  real(8),dimension(:),pointer          :: state_dvec
  real(8)                               :: state_e




contains



  !+------------------------------------------------------------------+
  !                        NORMAL
  !+------------------------------------------------------------------+
  !PURPOSE  : Evaluate the Green's function of the impurity electrons
  ! & phonons D = -<x(\tau)x(0)> with x = (b + b^+)
  subroutine build_gf_normal()
    integer                                     :: iorb,jorb,ispin,jspin,i
    logical                                     :: MaskBool
    logical(8),dimension(Nspin,Nspin,Norb,Norb) :: Hmask
    !

    do ispin=1,Nspin
       do iorb=1,Norb
          write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_s"//str(ispin)
          if(MPIMASTER)call start_timer
          call lanc_build_gf_normal_diag(iorb,ispin)
          if(MPIMASTER)call stop_timer(unit=LOGfile)
       enddo
    enddo
    !
    if(offdiag_gf_flag)then
       write(LOGfile,"(A)")""
       write(LOGfile,"(A)")"Get mask(G):"
       Hmask= .true.
       if(.not.ed_all_g)Hmask=mask_hloc(impHloc,wdiag=.true.,uplo=.false.)
       do ispin=1,Nspin
          do iorb=1,Norb
             write(LOGfile,*)((Hmask(ispin,jspin,iorb,jorb),jorb=1,Norb),jspin=1,Nspin)
          enddo
       enddo
       do ispin=1,Nspin
          do iorb=1,Norb
             do jorb=iorb+1,Norb
                MaskBool=.true.   
                if(bath_type=="replica")MaskBool=Hmask(ispin,ispin,iorb,jorb)
                if(.not.MaskBool)cycle
                !
                write(LOGfile,"(A)")"Get G_l"//str(iorb)//"_m"//str(jorb)//"_s"//str(ispin)
                if(MPIMASTER)call start_timer
                call lanc_build_gf_normal_mix(iorb,jorb,ispin)
                if(MPIMASTER)call stop_timer(unit=LOGfile)
             enddo
          enddo
       enddo
       !
       !
       !Put here off-diagonal manipulation by symmetry:
       select case(ed_diag_type)
       case default
          do ispin=1,Nspin
             do iorb=1,Norb
                do jorb=iorb+1,Norb
                   !if(hybrid)always T; if(replica)T iff following condition is T
                   MaskBool=.true.   
                   if(bath_type=="replica")MaskBool=Hmask(ispin,ispin,iorb,jorb)
                   !
                   if(.not.MaskBool)cycle
                   impGmats(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGmats(ispin,ispin,iorb,jorb,:) &
                        - impGmats(ispin,ispin,iorb,iorb,:) - impGmats(ispin,ispin,jorb,jorb,:))
                   impGreal(ispin,ispin,iorb,jorb,:) = 0.5d0*(impGreal(ispin,ispin,iorb,jorb,:) &
                        - impGreal(ispin,ispin,iorb,iorb,:) - impGreal(ispin,ispin,jorb,jorb,:))
                   impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
                   impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
                enddo
             enddo
          enddo
       case ("full")
          !>>ACTHUNG: this relation might not be true, it depends on the value of the impHloc_ij
          ! if impHloc_ij is REAL then it is true. if CMPLX hermiticity must be ensured
          impGmats(ispin,ispin,jorb,iorb,:) = impGmats(ispin,ispin,iorb,jorb,:)
          impGreal(ispin,ispin,jorb,iorb,:) = impGreal(ispin,ispin,iorb,jorb,:)
       end select
    end if
    !
    if(DimPh>1)then
       write(LOGfile,"(A)")"Get phonon Green function:"
       if(MPIMASTER)call start_timer()
       call lanc_build_gf_phonon_main()
       if(MPIMASTER)call stop_timer(unit=LOGfile)
    endif
  end subroutine build_gf_normal






  !################################################################
  !################################################################
  !################################################################
  !################################################################





  subroutine lanc_build_gf_normal_diag(iorb,ispin)
    integer,intent(in)          :: iorb,ispin
    type(sector)                :: sectorI,sectorJ
    !
    if(ed_total_ud)then
       ialfa = 1
       ipos  = iorb
    else
       ialfa = iorb
       ipos  = 1
    endif
    !
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_dvec => es_return_dvector(MpiComm,state_list,istate) 
       else
          state_dvec => es_return_dvector(state_list,istate)
       endif
#else
       state_dvec => es_return_dvector(state_list,istate)
#endif
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'From sector  :',isector,sectorI%Nups,sectorI%Ndws
       endif
       !
       !ADD ONE PARTICLE:
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then 
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c^+_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_dvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,iorb,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       endif
       !
       !REMOVE ONE PARTICLE:
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c_a,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_dvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,iorb,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_dvec))deallocate(state_dvec)
       else
          if(associated(state_dvec))nullify(state_dvec)
       endif
#else
       if(associated(state_dvec))nullify(state_dvec)
#endif
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_diag




  !################################################################






  subroutine lanc_build_gf_normal_mix(iorb,jorb,ispin)
    integer                     :: iorb,jorb,ispin
    type(sector)                :: sectorI,sectorJ
    !
    if(ed_total_ud)then
       ialfa = 1
       jalfa = ialfa               !this is the condition to evaluate G_ab: ialfa=jalfa
       ipos  = iorb
       jpos  = jorb
    else
       write(LOGfile,"(A)")"ED_GF_NORMAL warning: can not evaluate GF_ab with ed_total_ud=F"
       return
    endif
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_dvec => es_return_dvector(MpiComm,state_list,istate)
       else
          state_dvec => es_return_dvector(state_list,istate)
       endif
#else
       state_dvec => es_return_dvector(state_list,istate)
#endif
       !
       if(MpiMaster)then
          call build_sector(isector,sectorI)
          if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
               'From sector  :',isector,sectorI%Nups,sectorI%Ndws
       endif
       !
       !EVALUATE (c^+_iorb + c^+_jorb)|gs>
       jsector = getCDGsector(ialfa,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  ' apply c^+_a,s + c^+_b,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             !c^+_iorb|gs>
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_dvec(i)
             enddo
             !+c^+_jorb|gs>
             do i=1,sectorI%Dim
                call apply_op_CDG(i,j,sgn,jpos,jalfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_dvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,1,iorb,jorb,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       !EVALUATE (c_iorb + c_jorb)|gs>
       jsector = getCsector(ialfa,ispin,isector)
       if(jsector/=0)then
          if(MpiMaster)then
             call build_sector(jsector,sectorJ)
             if(ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")&
                  '  apply c_a,s + c_b,s:',jsector,sectorJ%Nups,sectorJ%Ndws
             allocate(vvinit(sectorJ%Dim)) ; vvinit=zero
             !c_iorb|gs>
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,ipos,ialfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = sgn*state_dvec(i)
             enddo
             !+c_jorb|gs>
             do i=1,sectorI%Dim
                call apply_op_C(i,j,sgn,jpos,jalfa,ispin,sectorI,sectorJ)
                if(sgn==0d0.OR.j==0)cycle
                vvinit(j) = vvinit(j) + sgn*state_dvec(i)
             enddo
             call delete_sector(sectorJ)
          else
             allocate(vvinit(1));vvinit=0.d0
          endif
          !
          call tridiag_Hv_sector_normal(jsector,vvinit,alfa_,beta_,norm2)
          call add_to_lanczos_gf_normal(one*norm2,state_e,alfa_,beta_,-1,iorb,jorb,ispin)
          deallocate(alfa_,beta_)
          if(allocated(vvinit))deallocate(vvinit)          
       endif
       !
       if(MpiMaster)call delete_sector(sectorI)
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_dvec))deallocate(state_dvec)
       else
          if(associated(state_dvec))nullify(state_dvec)
       endif
#else
       if(associated(state_dvec))nullify(state_dvec)
#endif
       !
    enddo
    return
  end subroutine lanc_build_gf_normal_mix





  !################################################################





  subroutine add_to_lanczos_gf_normal(vnorm2,Ei,alanc,blanc,isign,iorb,jorb,ispin)
    complex(8)                                 :: vnorm2,pesoBZ,peso
    real(8)                                    :: Ei,Egs,de
    integer                                    :: nlanc,itype
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    integer                                    :: isign,iorb,jorb,ispin
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j,ierr
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    if((finiteT).and.(beta*(Ei-Egs).lt.200))then
       pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    elseif(.not.finiteT)then
       pesoBZ = vnorm2/zeta_function
    else
       pesoBZ=0.d0
    endif
    !
    !pesoBZ = vnorm2/zeta_function
    !if(finiteT)pesoBZ = vnorm2*exp(-beta*(Ei-Egs))/zeta_function
    !
    !Only the nodes in Mpi_Comm_Group did get the alanc,blanc.
    !However after delete_sectorHv MpiComm returns to be the global one
    !so we can safely Bcast the alanc,blanc (known only to the operative group)
    !to every nodes. The master is in charge of this (as a
    !participant of the operative group)
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    !
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       de = diag(j)-Ei
       peso = pesoBZ*Z(1,j)*Z(1,j)
       do i=1,Lmats
          iw=xi*wm(i)
          impGmats(ispin,ispin,iorb,jorb,i)=impGmats(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
       do i=1,Lreal
          iw=dcmplx(wr(i),eps)
          impGreal(ispin,ispin,iorb,jorb,i)=impGreal(ispin,ispin,iorb,jorb,i) + peso/(iw-isign*de)
       enddo
    enddo
  end subroutine add_to_lanczos_gf_normal






  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################



  subroutine lanc_build_gf_phonon_main()
    integer,dimension(Ns_Ud)    :: iDimUps,iDimDws
    integer                     :: Nups(Ns_Ud)
    integer                     :: Ndws(Ns_Ud)
    !
    do istate=1,state_list%size
       isector    =  es_return_sector(state_list,istate)
       state_e    =  es_return_energy(state_list,istate)
#ifdef _MPI
       if(MpiStatus)then
          state_dvec => es_return_dvector(MpiComm,state_list,istate)
       else
          state_dvec => es_return_dvector(state_list,istate)
       endif
#else
       state_dvec => es_return_dvector(state_list,istate)
#endif
       !
       call get_Nup(isector,Nups)
       call get_Ndw(isector,Ndws)
       if(MpiMaster.AND.ed_verbose>=3)write(LOGfile,"(A,I6,20I4)")'From sector:',isector,Nups,Ndws
       !
       idim = getdim(isector)
       call get_DimUp(isector,iDimUps)
       call get_DimDw(isector,iDimDws)
       iDimUp = product(iDimUps)
       iDimDw = product(iDimDws)
       !
       if(MpiMaster)then
          if(ed_verbose==3)write(LOGfile,"(A,I12)")'Apply x:',isector
          !
          allocate(vvinit(idim));vvinit=0d0
          !
          do i=1,iDim
             iph = (i-1)/(iDimUp*iDimDw) + 1
             i_el = mod(i-1,iDimUp*iDimDw) + 1
             !
             !apply destruction operator
             if(iph>1) then
                j = i_el + ((iph-1)-1)*iDimUp*iDimDw
                vvinit(j) = vvinit(j) + sqrt(dble(iph-1))*state_dvec(i)
             endif
             !
             !apply creation operator
             if(iph<DimPh) then
                j = i_el + ((iph+1)-1)*iDimUp*iDimDw
                vvinit(j) = vvinit(j) + sqrt(dble(iph))*state_dvec(i)
             endif
          enddo
       else
          allocate(vvinit(1));vvinit=0.d0
       endif
       !
       call tridiag_Hv_sector_normal(isector,vvinit,alfa_,beta_,norm2)
       call add_to_lanczos_phonon(norm2,state_e,alfa_,beta_)
       deallocate(alfa_,beta_)
       if(allocated(vvinit))deallocate(vvinit)
#ifdef _MPI
       if(MpiStatus)then
          if(associated(state_dvec))deallocate(state_dvec)
       else
          if(associated(state_dvec))nullify(state_dvec)
       endif
#else
       if(associated(state_dvec))nullify(state_dvec)
#endif
    enddo
    return
  end subroutine lanc_build_gf_phonon_main


  !################################################################


  subroutine add_to_lanczos_phonon(vnorm2,Ei,alanc,blanc)
    real(8)                                    :: vnorm2,Ei,Ej,Egs,pesoF,pesoAB,pesoBZ,de,peso
    integer                                    :: nlanc
    real(8),dimension(:)                       :: alanc
    real(8),dimension(size(alanc))             :: blanc 
    real(8),dimension(size(alanc),size(alanc)) :: Z
    real(8),dimension(size(alanc))             :: diag,subdiag
    integer                                    :: i,j
    complex(8)                                 :: iw
    !
    Egs = state_list%emin       !get the gs energy
    !
    Nlanc = size(alanc)
    !
    pesoF  = vnorm2/zeta_function 
    pesoBZ = 1d0
    if(finiteT)pesoBZ = exp(-beta*(Ei-Egs))
    !
#ifdef _MPI
    if(MpiStatus)then
       call Bcast_MPI(MpiComm,alanc)
       call Bcast_MPI(MpiComm,blanc)
    endif
#endif
    diag(1:Nlanc)    = alanc(1:Nlanc)
    subdiag(2:Nlanc) = blanc(2:Nlanc)
    call eigh(diag(1:Nlanc),subdiag(2:Nlanc),Ev=Z(:Nlanc,:Nlanc))
    !
    do j=1,nlanc
       Ej     = diag(j)
       dE     = Ej-Ei
       pesoAB = Z(1,j)*Z(1,j)
       peso   = pesoF*pesoAB*pesoBZ
       ! the correct behavior for beta*dE << 1 is recovered only by assuming that v_n is still finite
       ! beta*dE << v_n for v_n--> 0 slower. First limit beta*dE--> 0 and only then v_n -->0.
       ! This ensures that the correct null contribution is obtained.
       ! So we impose that: if (beta*dE is larger than a small qty) we sum up the contribution, else
       ! we do not include the contribution (because we are in the situation described above).
       ! For the real-axis case this problem is circumvented by the usual i*0+ = xi*eps
       if(beta*dE > 1d-3)impDmats_ph(0)=impDmats_ph(0) - peso*2*(1d0-exp(-beta*dE))/dE 
       do i=1,Lmats
          impDmats_ph(i)=impDmats_ph(i) - peso*(1d0-exp(-beta*dE))*2d0*dE/(vm(i)**2+dE**2)
       enddo
       do i=1,Lreal
          impDreal_ph(i)=impDreal_ph(i) + peso*(1d0-exp(-beta*dE))*(1d0/(dcmplx(vr(i),eps) - dE) - 1d0/(dcmplx(vr(i),eps) + dE))
       enddo
    enddo
  end subroutine add_to_lanczos_phonon




  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################
  !############################################################################################






  subroutine build_sigma_normal
    integer                                           :: i,ispin,iorb
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats) :: invG0mats,invGmats
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lreal) :: invG0real,invGreal
    complex(8),dimension(Norb,Norb)                   :: invGimp
    !
    invG0mats = zero
    invGmats  = zero
    invG0real = zero
    invGreal  = zero
    !
    !Get G0^-1
    invG0mats(:,:,:,:,:) = invg0_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    invG0real(:,:,:,:,:) = invg0_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !
    select case(bath_type)
    case default                !Diagonal in both spin and orbital
       !
       !Get Gimp^-1
       do ispin=1,Nspin
          do iorb=1,Norb
             invGmats(ispin,ispin,iorb,iorb,:) = one/impGmats(ispin,ispin,iorb,iorb,:)
             invGreal(ispin,ispin,iorb,iorb,:) = one/impGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             impSmats(ispin,ispin,iorb,iorb,:) = invG0mats(ispin,ispin,iorb,iorb,:) - invGmats(ispin,ispin,iorb,iorb,:)
             impSreal(ispin,ispin,iorb,iorb,:) = invG0real(ispin,ispin,iorb,iorb,:) - invGreal(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !
    case ("hybrid","replica")   !Diagonal in spin only. Full Orbital structure
       !
       !Get Gimp^-1
       do ispin=1,Nspin
          do i=1,Lmats
             invGimp = impGmats(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGmats(ispin,ispin,:,:,i)=invGimp
          enddo
          !
          do i=1,Lreal
             invGimp = impGreal(ispin,ispin,:,:,i)
             call inv(invGimp)
             invGreal(ispin,ispin,:,:,i)=invGimp
          enddo
       enddo
       !Get Sigma functions: Sigma= G0^-1 - G^-1
       impSmats=zero
       impSreal=zero
       do ispin=1,Nspin
          impSmats(ispin,ispin,:,:,:) = invG0mats(ispin,ispin,:,:,:) - invGmats(ispin,ispin,:,:,:)
          !
          impSreal(ispin,ispin,:,:,:) = invG0real(ispin,ispin,:,:,:) - invGreal(ispin,ispin,:,:,:)
       enddo
       !
    end select
    !
    !Get G0and:
    impG0mats(:,:,:,:,:) = g0and_bath_function(dcmplx(0d0,wm(:)),dmft_bath)
    impG0real(:,:,:,:,:) = g0and_bath_function(dcmplx(wr(:),eps),dmft_bath)
    !!
    !
  end subroutine build_sigma_normal


END MODULE ED_GF_NORMAL











