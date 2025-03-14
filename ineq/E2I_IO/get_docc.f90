
subroutine ed_get_docc_n1(self,Nlat,iorb)
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  integer              :: iorb
  if(iorb>Norb)stop "ed_get_docc error: orbital index > N_orbital"
  if(.not.allocated(docc_ineq))stop "ed_get_docc error: docc_ineq not allocated"
  if(Nlat>size(docc_ineq,1))stop "ed_get_docc error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_docc','docc')
  self = docc_ineq(:,iorb)
end subroutine ed_get_docc_n1

subroutine ed_get_docc_n2(self,Nlat)
  real(8),dimension(:,:) :: self
  integer                :: Nlat
  if(.not.allocated(docc_ineq))stop "ed_get_docc error: docc_ineq not allocated"
  if(Nlat>size(docc_ineq,1))stop "ed_get_docc error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,Norb],'ed_get_docc','docc')
  self = docc_ineq
end subroutine ed_get_docc_n2
