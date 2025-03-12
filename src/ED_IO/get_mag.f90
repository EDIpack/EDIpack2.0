subroutine ed_get_mag_n0(self,component,iorb)
  real(8) :: self !Magnetization
  character(len=1),optional :: component !Component of the magnetization, can be :code:`"x"`, :code:`"y"`, :code:`"z"` (default :code:`"z"` )
  integer,optional          :: iorb !Orbital (default :code:`1`)
  !
  integer                   :: iorb_
  character(len=1)          :: char
  integer                   :: id
  !
  iorb_=1;if(present(iorb))iorb_=iorb
  char='Z';if(present(component))char=component
  select case(char)
  case default;id=3
  case('X','x');id=1
  case('Y','y');id=2
  end select
  if(iorb_>Norb)stop "ed_get_mag error: orbital index > N_orbital"
  self = ed_mag(id,iorb_)
end subroutine ed_get_mag_n0

subroutine ed_get_mag_n1(self,component)
  real(8),dimension(:)      :: self
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
  call assert_shape(self,[Norb],'ed_get_mag','mag')
  self = ed_mag(id,:)
end subroutine ed_get_mag_n1


subroutine ed_get_mag_n2(self)
  real(8),dimension(:,:)      :: self
  !
  call assert_shape(self,shape(ed_mag),'ed_get_mag','mag')
  self = ed_mag
end subroutine ed_get_mag_n2
