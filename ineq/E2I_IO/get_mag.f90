subroutine ed_get_mag_n1(self,Nlat,iorb,component)
  real(8),dimension(:)      :: self
  integer                   :: iorb
  integer                   :: Nlat !Number of inequivalent impurities for real-space DMFT
  character(len=1),optional :: component
  character(len=1)          :: char
  integer                   :: id
  !
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  if(.not.allocated(mag_ineq))stop "ed_get_mag error: mag_ineq not allocated"
  if(Nlat>size(mag_ineq,1))stop "ed_get_mag error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat],'ed_get_mag','mag')
  self = mag_ineq(:,id,iorb)
end subroutine ed_get_mag_n1

subroutine ed_get_mag_n2(self,Nlat,component)
  real(8),dimension(:,:)    :: self
  integer                   :: Nlat
  character(len=1),optional :: component
  character(len=1)          :: char
  integer                   :: id
  !
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  if(.not.allocated(mag_ineq))stop "ed_get_mag error: mag_ineq not allocated"
  if(Nlat>size(mag_ineq,1))stop "ed_get_mag error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,Norb],'ed_get_mag','mag')
  self = mag_ineq(:,id,:)
end subroutine ed_get_mag_n2

subroutine ed_get_mag_n3(self,Nlat,component)
  real(8),dimension(:,:,:)  :: self
  integer                   :: Nlat
  character(len=1),optional :: component
  character(len=1)          :: char
  integer                   :: id
  !
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  if(.not.allocated(mag_ineq))stop "ed_get_mag error: mag_ineq not allocated"
  if(Nlat>size(mag_ineq,1))stop "ed_get_mag error: required N_sites > evaluated N_sites"
  call assert_shape(self,[Nlat,3,Norb],'ed_get_mag','mag')
  self = mag_ineq
end subroutine ed_get_mag_n3
