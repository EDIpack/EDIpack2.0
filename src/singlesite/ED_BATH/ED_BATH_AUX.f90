MODULE ED_BATH_AUX
  !:synopsis: Auxiliary routines for bath creation and manipulation
  !Implements a number of auxiliary procedures used in the bath handling
  !
  USE SF_CONSTANTS, only: zero
  USE SF_IOTOOLS, only:free_unit,reg,file_length,str
  USE SF_LINALG, only: eye,inv
  USE SF_ARRAYS, only:linspace
  USE SF_MISC, only: assert_shape
  USE ED_INPUT_VARS
  USE ED_VARS_GLOBAL
  USE ED_AUX_FUNX
  implicit none



  interface get_Whyb_matrix
     ! This subroutine build up the hybridization matrix :math:`W_{\sigma\sigma^{'}}` used in the :f:var:`ed_mode` = nonsu2 with :f:var:`bath_type` = hybrid. The input can have different shape and type:
     !
     !   * :f:var:`u` , :f:var:`v` with dimensions [ |Nspin| , |Norb| ]  
     !   * :f:var:`u` , :f:var:`v` with dimensions [ |Nspin| , |Norb| , :f:var:`nbath`]  
     !   * :f:var:`dmft_bath` 
     !
     module procedure get_Whyb_matrix_1orb
     module procedure get_Whyb_matrix_Aorb
     module procedure get_Whyb_matrix_dmft_bath
  end interface get_Whyb_matrix


  interface is_identity
     !
     ! This subroutine checks if a matrix :math:`\hat{O}`  in the basis of the :code:`replica` or :code:`general` baths is the identity.
     ! 
     ! The input matrix can have different shapes:
     !    *  [ |Nnambu| . |Nspin| . |Norb| , |Nnambu| . |Nspin| . |Norb| ]
     !    *  [ |Nnambu| . |Nspin| , |Nnambu| . |Nspin| , |Norb| , |Norb| ]
     !
     module procedure ::  is_identity_so
     module procedure ::  is_identity_nn
  end interface is_identity

  interface is_diagonal
     !
     ! This subroutine checks if a matrix :math:`\hat{O}`  in the basis of the :code:`replica` or :code:`general` baths is diagonal.
     ! 
     ! The input matrix can have different shapes:
     !    *  [ |Nnambu| . |Nspin| . |Norb| , |Nnambu| . |Nspin| . |Norb| ]
     !    *  [ |Nnambu| . |Nspin| , |Nnambu| . |Nspin| , |Norb| , |Norb| ]
     !
     module procedure ::  is_diagonal_so
     module procedure ::  is_diagonal_nn
  end interface is_diagonal




contains



  function get_Whyb_matrix_1orb(v,u) result(w)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    real(8),dimension(Nspin,Nbath)       :: v,u
    real(8),dimension(Nspin,Nspin,Nbath) :: w
    integer                              :: ispin
    do ispin=1,Nspin
       w(ispin,ispin,:) = v(ispin,:)
    enddo
    w(1,Nspin,:) = u(1,:)
    w(Nspin,1,:) = u(2,:)
  end function get_Whyb_matrix_1orb

  function get_Whyb_matrix_Aorb(v,u) result(w)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    real(8),dimension(Nspin,Norb,Nbath)       :: v,u
    real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
    integer                                   :: ispin
    do ispin=1,Nspin
       w(ispin,ispin,:,:) = v(ispin,:,:)
    enddo
    w(1,Nspin,:,:) = u(1,:,:)
    w(Nspin,1,:,:) = u(2,:,:)
  end function get_Whyb_matrix_Aorb

  function get_Whyb_matrix_dmft_bath(dmft_bath_) result(w)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    type(effective_bath)                      :: dmft_bath_
    real(8),dimension(Nspin,Nspin,Norb,Nbath) :: w
    integer                                   :: ispin
    !
    do ispin=1,Nspin
       w(ispin,ispin,:,:) = dmft_bath_%v(ispin,:,:)
    enddo
    w(1,Nspin,:,:) = dmft_bath_%u(1,:,:)
    w(Nspin,1,:,:) = dmft_bath_%u(2,:,:)
  end function get_Whyb_matrix_dmft_bath




  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is the identity
  !+-------------------------------------------------------------------+
  function is_identity_nn(mnnn) result(flag)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: mnnn
    real(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=nn2so_reshape(mnnn,Nnambu*Nspin,Norb)
    !
    do i=1,Nnambu*Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_nn

  function is_identity_so(mlso) result(flag)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mlso
    real(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                  :: i,j
    logical                                  :: flag
    !
    flag=.true.
    !
    mtmp=mlso
    !
    do i=1,Nnambu*Nspin*Norb-1
       if((mtmp(i,i).ne.mtmp(i+1,i+1)).or.(mtmp(i,i).lt.1.d-6))flag=.false.
    enddo
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(mtmp(i,j).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_identity_so



  !+-------------------------------------------------------------------+
  !PURPOSE  : Check if a matrix is diagonal
  !+-------------------------------------------------------------------+
  function is_diagonal_nn(mnnn) result(flag)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(Nnambu*Nspin,Nnambu*Nspin,Norb,Norb) :: mnnn
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                     :: i,j
    logical                                     :: flag
    !
    flag=.true.
    !
    mtmp=abs((nn2so_reshape(mnnn,Nnambu*Nspin,Norb)))
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(abs(mtmp(i,j)).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_nn

  function is_diagonal_so(mlso) result(flag)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mlso
    complex(8),dimension(Nnambu*Nspin*Norb,Nnambu*Nspin*Norb) :: mtmp
    integer                                     :: i,j
    logical                                     :: flag
    !
    flag=.true.
    !
    mtmp=abs((mlso))
    !
    do i=1,Nnambu*Nspin*Norb
       do j=1,Nnambu*Nspin*Norb
          if((i.ne.j).and.(abs(mtmp(i,j)).gt.1.d-6))flag=.false.
       enddo
    enddo
    !
  end function is_diagonal_so






  function check_herm(A,N,error) result(bool)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    integer,intent(in)                   :: N
    complex(8),dimension(N,N),intent(in) :: A
    logical                              :: bool
    real(8),optional                     :: error
    real(8)                              :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    bool   = all(abs(A - conjg(transpose(A)))<error_)
  end function check_herm


  function check_nambu(A,N,error) result(bool)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb,Nbath
#endif
    integer,intent(in)                       :: N
    complex(8),dimension(2*N,2*N),intent(in) :: A
    complex(8),dimension(N,N)                :: h11,h22
    logical                                  :: bool
    real(8),optional                         :: error
    real(8)                                  :: error_
    error_ = 1d-6 ; if(present(error))error_=error
    h11    = A(1:N    ,1:N)
    h22    = A(N+1:2*N,N+1:2*N)
    bool   = check_herm(A,2*N,error_) !this checks also for F = A_12, s.t. A_21=herm(A_12)
    bool   = bool.AND.( all(abs(h22 + conjg(h11))<error_) )
  end function check_nambu





END MODULE ED_BATH_AUX

