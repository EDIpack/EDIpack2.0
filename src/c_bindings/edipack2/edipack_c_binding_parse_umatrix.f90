  !Inspect interaction parameters
  subroutine inspect_uparams_c(which_param, coeffmatrix, ioflag, ierrflag) bind(c, name='inspect_uparams')
    use, intrinsic :: iso_c_binding
    use ED_PARSE_UMATRIX, only: inspect_uparams
    integer(c_int),value                :: which_param, ioflag
    integer(c_int),dimension(1)         :: ierrflag
    real(c_double),dimension(Norb,Norb) :: coeffmatrix
    !
    call inspect_uparams(which_param, coeffmatrix, ioflag, ierrflag(1))
  end subroutine inspect_uparams_c
  
  
  subroutine reset_umatrix_c() bind(c, name='reset_umatrix')
    use, intrinsic :: iso_c_binding
    use ED_PARSE_UMATRIX, only: reset_umatrix
    call reset_umatrix()
  end subroutine reset_umatrix_c
    

  subroutine add_twobody_operator_c(o1,s1,o2,s2,o3,s3,o4,s4,U) bind(c, name='add_twobody_operator')
    use, intrinsic :: iso_c_binding
    use ED_PARSE_UMATRIX, only: add_twobody_operator
    integer(c_int),value          :: o1, o2, o3, o4, s1, s2, s3, s4
    character(len=1)              :: s1_, s2_, s3_, s4_
    real(c_double),value          :: U
    
    s1_ = merge("u", "d", s1 == 1)
    s2_ = merge("u", "d", s2 == 1)
    s3_ = merge("u", "d", s3 == 1)
    s4_ = merge("u", "d", s4 == 1)
    
    call add_twobody_operator(o1,s1_,o2,s2_,o3,s3_,o4,s4_,U)

  end subroutine add_twobody_operator_c
