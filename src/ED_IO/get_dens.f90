subroutine ed_get_dens_n0(self,iorb)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !The density value or array of values
  integer,optional      :: iorb !the orbital index
  integer               :: iorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  if(iorb_>Norb)stop "ed_get_dens error: orbital index > N_orbital"
  self = ed_dens(iorb_)
end subroutine ed_get_dens_n0


subroutine ed_get_dens_n1(self,iorb,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self
  integer,optional     :: iorb 
  integer,optional     :: Nlat !the number of inequivalent impurity sites for real-space DMFT
  integer              :: iorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  if(iorb_>Norb)stop "ed_get_dens error: orbital index > N_orbital"
  if(present(Nlat))then
     if(.not.allocated(dens_ineq))stop "ed_get_dens error: dens_ineq not allocated"
     if(Nlat>size(dens_ineq,1))stop "ed_get_dens error: required N_sites > evaluated N_sites"
  endif
  if(present(Nlat))then
     call assert_shape(self,[Nlat],'ed_get_dens','dens')
     self = dens_ineq(:,iorb_)
  else
     call assert_shape(self,[Norb],'ed_get_dens','dens')
     self = ed_dens
  endif
end subroutine ed_get_dens_n1


subroutine ed_get_dens_n2(self,Nlat)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:) :: self
  integer                :: Nlat
  if(.not.allocated(dens_ineq))stop "ed_get_dens error: dens_ineq not allocated"
  if(Nlat>size(dens_ineq,1))stop "ed_get_dens error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,Norb],'ed_get_dens','dens')
  self = dens_ineq
end subroutine ed_get_dens_n2





subroutine ed_get_imp_info(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(2) :: self
  self = ed_imp_info
end subroutine ed_get_imp_info
