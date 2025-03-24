
subroutine ed_get_doubles_n2(self,Nlat)
  real(8),dimension(:,:) :: self
  integer                :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_doubles error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_doubles error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,4],'ed_get_doubles','doubles')
  self = dd_ineq
end subroutine ed_get_doubles_n2

subroutine ed_get_dust_n1(self,Nlat)
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_dust error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_dust error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dust','dust')
  self = dd_ineq(:,1)
end subroutine ed_get_dust_n1

subroutine ed_get_dund_n1(self,Nlat)
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_dund error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_dund error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dund','dund')
  self = dd_ineq(:,2)
end subroutine ed_get_dund_n1

subroutine ed_get_dse_n1(self,Nlat)
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_dse error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_dse error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dse','dse')
  self = dd_ineq(:,3)
end subroutine ed_get_dse_n1

subroutine ed_get_dph_n1(self,Nlat)
  real(8),dimension(:) :: self
  integer              :: Nlat !number of inequivalent impurity sites for real-space DMFT
  !
  if(.not.allocated(dd_ineq))stop "ed_get_dph error: dd_ineq not allocated"
  if(Nlat>size(dd_ineq,1))stop "ed_get_dph error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_dph','dph')
  self = dd_ineq(:,4)
end subroutine ed_get_dph_n1
