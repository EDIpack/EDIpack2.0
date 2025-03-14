subroutine ed_get_doubles_n1(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8),dimension(:) :: self !array of two-body terms expectation values
  call assert_shape(self,[4],'ed_get_doubles','doubles')
  self = [ed_Dust,ed_Dund,ed_Dse,ed_Dph]
end subroutine ed_get_doubles_n1


subroutine ed_get_dust_n0(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !value of :f:var:`dust`
  self = ed_Dust
end subroutine ed_get_dust_n0


subroutine ed_get_dund_n0(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !value of :f:var:`dund`
  self = ed_Dund
end subroutine ed_get_dund_n0


subroutine ed_get_dse_n0(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !value of :f:var:`dse`
  self = ed_Dse
end subroutine ed_get_dse_n0


subroutine ed_get_dph_n0(self)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  real(8) :: self !value of :f:var:`dph`
  self = ed_Dph
end subroutine ed_get_dph_n0
