subroutine ed_get_eimp_n1(self)
  real(8),dimension(:) :: self !energy components array
  call assert_shape(self,[4],'ed_get_eimp','eimp')
  self = [ed_Epot,ed_Eint,ed_Ehartree,ed_Eknot]
end subroutine ed_get_eimp_n1




subroutine ed_get_epot_n0(self)
  real(8) :: self !value of :f:var:`ed_epot`
  self = ed_Epot
end subroutine ed_get_epot_n0




subroutine ed_get_eint_n0(self)
  real(8) :: self !value of :f:var:`ed_int`
  self = ed_Eint
end subroutine ed_get_eint_n0





subroutine ed_get_ehartree_n0(self)
  real(8) :: self !value of :f:var:`ed_ehartree`
  self = ed_Ehartree
end subroutine ed_get_ehartree_n0




subroutine ed_get_eknot_n0(self)
  real(8)          :: self !value of :f:var:`ed_eknot`
  self = ed_Eknot
end subroutine ed_get_eknot_n0


