subroutine ed_get_exct_n0(self,component,iorb,jorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8)          :: self ! :math:`{S_0,T_x,T_y,T_z}` value 
  character(len=1) :: component
  integer,optional :: iorb ! first orbital index
  integer,optional :: jorb ! second orbital index
  integer          :: iorb_,jorb_,id
  self=0d0
  if(Norb<2)return
  iorb_=1;if(present(iorb))iorb_=iorb
  jorb_=2;if(present(jorb))jorb_=jorb
  if(iorb_>Norb.OR.jorb_>Norb)stop "ed_get_exct error: orbital index > N_orbital"
  select case(component)
  case default;stop "ed_get_exct error: wrong component, not in [0,x,y,z]"
  case('0','S','s');id=1
  case('X','x');id=2
  case('Y','y');id=3
  case('Z','z');id=4
  end select
  self = ed_exct(id,iorb_,jorb_)
end subroutine ed_get_exct_n0

subroutine ed_get_exct_n1(self,iorb,jorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(4) :: self ! :math:`[S_0,T_x,T_y,T_z]` array
  integer,optional     :: iorb ! first orbital index
  integer,optional     :: jorb ! second orbital index
  integer              :: iorb_,jorb_
  self=0d0
  if(Norb<2)return
  iorb_=1;if(present(iorb))iorb_=iorb
  jorb_=2;if(present(jorb))jorb_=jorb
  if(iorb_>Norb.OR.jorb_>Norb)stop "ed_get_exct error: orbital index > N_orbital"
  self = ed_exct(1:4,iorb_,jorb_)
end subroutine ed_get_exct_n1


!X_ab
subroutine ed_get_exct_n2(self,component)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:) :: self
  character(len=1)       :: component
  integer                :: id
  call assert_shape(self,[Norb,Norb],'ed_get_exct','exct')
  select case(component)
  case default;stop "ed_get_exct error: wrong component, not in [0,x,y,z]"
  case('0','S','s');id=1
  case('X','x');id=2
  case('Y','y');id=3
  case('Z','z');id=4
  end select
  self = ed_exct(id,:,:)
end subroutine ed_get_exct_n2


!\vecX_ab
subroutine ed_get_exct_n3(self)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:,:) :: self
  call assert_shape(self,[4,Norb,Norb],'ed_get_exct','exct')
  self = ed_exct
end subroutine ed_get_exct_n3
