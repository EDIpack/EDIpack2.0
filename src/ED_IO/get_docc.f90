subroutine ed_get_docc_n0(self,iorb)
  real(8)          :: self !double-occupation value or array of values
  integer,optional :: iorb !orbital index
  integer          :: iorb_
  iorb_=1;if(present(iorb))iorb_=iorb
  if(iorb_>Norb)stop "ed_get_docc error: orbital index > N_orbital"
  self = ed_docc(iorb_)
end subroutine ed_get_docc_n0

subroutine ed_get_docc_n1(self)
  real(8),dimension(:) :: self
  call assert_shape(self,[Norb],'ed_get_docc','docc')
  self = ed_docc
end subroutine ed_get_docc_n1
