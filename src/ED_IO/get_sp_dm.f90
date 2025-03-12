subroutine get_density_matrix_n2(dm,iprint)!,custom_rot,dm_eig_,dm_rot_)
  !passed
  complex(8),intent(out) :: dm(:,:)
  logical,optional       :: iprint
  logical                :: iprint_
  integer                :: unit
  integer                :: iorb,jorb,ispin,jspin,io,jo
  complex(8)             :: Tr
  complex(8),allocatable :: dm_custom_rot(:,:)
  real(8)                :: soc
  !
  iprint_=.false.;if(present(iprint))iprint_=iprint
  !
  if (.not.allocated(single_particle_density_matrix)) then
     write(LOGfile,"(A)") "single_particle_density_matrix is not allocated"
     stop
  endif
  !
  ! if(allocated(dm_)) deallocate(dm_)
  ! allocate(dm_(Nspin*Norb,Nspin*Norb)); dm_ = zero
  call assert_shape(dm,[Nspin*Norb,Nspin*Norb],"get_density_matrix_single","dm")
  ! dm in the impurity problem basis
  dm = nn2so_reshape(single_particle_density_matrix,Nspin,Norb)
  if(iprint_)call print_dm(dm,1)
end subroutine get_density_matrix_n2

subroutine get_density_matrix_n4(dm,iprint)!,custom_rot,dm_eig_,dm_rot_)
  !passed
  complex(8),intent(out) :: dm(:,:,:,:)
  logical,optional       :: iprint
  logical                :: iprint_
  integer                :: unit
  integer                :: iorb,jorb,ispin,jspin,io,jo
  complex(8)             :: Tr
  complex(8),allocatable :: dm_custom_rot(:,:)
  real(8)                :: soc
  !
  iprint_=.false.;if(present(iprint))iprint_=iprint
  !
  if (.not.allocated(single_particle_density_matrix)) then
     write(LOGfile,"(A)") "single_particle_density_matrix is not allocated"
     stop
  endif
  !
  ! if(allocated(dm_)) deallocate(dm_)
  ! allocate(dm_(Nspin*Norb,Nspin*Norb)); dm_ = zero
  call assert_shape(dm,[Nspin,Nspin,Norb,Norb],"get_density_matrix_single","dm")
  ! dm in the impurity problem basis
  dm = single_particle_density_matrix
  if(iprint_)call print_dm(dm,1)
end subroutine get_density_matrix_n4


subroutine print_dm(dm_,ndx)
  implicit none
  integer               ,intent(in) :: ndx
  complex(8),intent(in)             :: dm_(Nspin*Norb,Nspin*Norb)
  integer                           :: unit
  character(len=24)                 :: suffix
  integer                           :: iorb,jorb,ispin,jspin,io,jo
  !
  suffix="single_particle_density_matrix_"//reg(str(ndx))//".dat"
  !
  unit = free_unit()
  open(unit,file=suffix,action="write",position="rewind",status='unknown')
  !
  do io=1,Nspin*Norb
     write(unit,"(*(F15.9,1X))") (dreal(dm_(io,jo)),jo=1,Nspin*Norb)
  enddo
  write(unit,*)
  !
  do io=1,Nspin*Norb
     write(unit,"(*(F15.9,1X))") (dimag(dm_(io,jo)),jo=1,Nspin*Norb)
  enddo
  write(unit,*)
  close(unit)
  !
end subroutine print_dm
