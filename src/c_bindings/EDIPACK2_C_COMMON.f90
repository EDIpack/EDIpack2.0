MODULE EDIPACK2_C_COMMON
!:synopsis: Fortran-C bindings: common routines and variables
  USE EDIPACK2
  USE SCIFOR
  USE ISO_C_BINDING
  implicit none

#ifdef _WINEQ
  integer(c_int),bind(c, name="has_ineq") :: has_ineq=1 
#else
  integer(c_int),bind(c, name="has_ineq") :: has_ineq=0
#endif

contains

  !integer to logical
  function i2l(var_integer) result (var_logical)
    use, intrinsic :: iso_c_binding
    integer        :: var_integer
    logical        :: var_logical   

    if (var_integer == 1) then
       var_logical = .true.
    else
       var_logical = .false.
    endif
  end function i2l

  !logical to integer
  function l2i(var_logical) result (var_integer)
    use, intrinsic :: iso_c_binding
    integer    :: var_integer
    logical    :: var_logical   

    if (var_logical) then
       var_integer = 1
    else
       var_integer = 0
    endif
  end function l2i

  !c string to fortran string
  subroutine c2f(c_str)
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*),intent(IN) :: c_str
    character(len=120), allocatable                 :: f_str
    integer                                         :: length
    integer                                         :: i

    length=0
    f_str=" "
    do
       if (c_str(length+1) == C_NULL_CHAR) exit
       length = length + 1
    end do
    do i = 1, length
       f_str(i:i) = c_str(i)
    enddo
    f_str=trim(f_str)
  end subroutine c2f

  !Get Bath Type
  integer(c_int) function get_bath_type_c() result(bt) bind(c, name='get_bath_type')
    use, intrinsic :: iso_c_binding
    select case(bath_type)
    case("normal");  bt = 1
    case("hybrid");  bt = 2
    case("replica"); bt = 3
    case("general"); bt = 4
    end select
  end function get_bath_type_c

  !Get Bath Type
  integer(c_int) function get_ed_mode_c() result(edm) bind(c, name='get_ed_mode')
    use, intrinsic :: iso_c_binding
    select case(ed_mode)
    case("normal");  edm = 1
    case("superc");  edm = 2
    case("nonsu2");  edm = 3
    end select
  end function get_ed_mode_c
  
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

END MODULE EDIPACK2_C_COMMON
