subroutine ed_get_exct_n0(self,component,iorb,jorb)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8)          :: self ! :math:`[S_0,T_z]` value or array of values
  integer          :: component
  integer,optional :: iorb ! first orbital index
  integer,optional :: jorb ! second orbital index
  integer          :: iorb_,jorb_
  self=0d0
  if(Norb<2)return
  iorb_=1;if(present(iorb))iorb_=iorb
  jorb_=2;if(present(jorb))jorb_=jorb
  if(iorb_>Norb.OR.jorb_>Norb)stop "ed_get_exct error: orbital index > N_orbital"
  if(component<0.OR.component>4)stop "ed_get_exct error: component index > 4 or < 0"
  self = ed_exct(component,iorb_,jorb_)
end subroutine ed_get_exct_n0

!X_ab
subroutine ed_get_exct_n1(self,component)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:) :: self
  integer                :: component
  call assert_shape(self,[Norb,Norb],'ed_get_exct','exct')
  self = ed_exct(component,:,:)
end subroutine ed_get_exct_n1


!\vecX_ab
subroutine ed_get_exct_n2(self)
#if __INTEL_COMPILER
  use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:,:,:) :: self
  call assert_shape(self,[4,Norb,Norb],'ed_get_exct','exct')
  self = ed_exct
end subroutine ed_get_exct_n2
