subroutine ed_get_dens_n0(self,iorb)
  real(8) :: self !The density value or array of values
  integer,optional      :: iorb !the orbital index
  integer               :: iorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  if(iorb_>Norb)stop "ed_get_dens error: orbital index > N_orbital"
  self = ed_dens(iorb_)
end subroutine ed_get_dens_n0


subroutine ed_get_dens_n1(self)
  real(8),dimension(:) :: self
  call assert_shape(self,[Norb],'ed_get_dens','dens')
  self = ed_dens
end subroutine ed_get_dens_n1



subroutine ed_get_imp_info(self)
  real(8),dimension(2) :: self
  self = ed_imp_info
end subroutine ed_get_imp_info
