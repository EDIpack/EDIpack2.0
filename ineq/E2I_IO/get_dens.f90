subroutine ed_get_dens_n1(self,Nlat,iorb)
  real(8),dimension(:) :: self
  integer              :: Nlat !the number of inequivalent impurity sites for real-space DMFT
  integer              :: iorb 
  if(.not.allocated(dens_ineq))stop "ed_get_dens error: dens_ineq not allocated"
  if(Nlat>size(dens_ineq,1))stop "ed_get_dens error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dens','dens')
  self = dens_ineq(:,iorb)
end subroutine ed_get_dens_n1


subroutine ed_get_dens_n2(self,Nlat)
  real(8),dimension(:,:) :: self
  integer                :: Nlat
  if(.not.allocated(dens_ineq))stop "ed_get_dens error: dens_ineq not allocated"
  if(Nlat>size(dens_ineq,1))stop "ed_get_dens error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,Norb],'ed_get_dens','dens')
  self = dens_ineq
end subroutine ed_get_dens_n2
