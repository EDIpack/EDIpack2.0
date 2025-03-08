subroutine get_density_matrix_single(dm_)!,custom_rot,dm_eig_,dm_rot_)
  !passed
  complex(8),allocatable,intent(out)           :: dm_(:,:)
  ! complex(8),allocatable,intent(in) ,optional  :: custom_rot(:,:)
  ! real(8),allocatable,intent(out)   ,optional  :: dm_eig_(:)
  ! complex(8),allocatable,intent(out),optional  :: dm_rot_(:,:)
  !internal
  integer                                      :: unit
  integer                                      :: iorb,jorb,ispin,jspin,io,jo
  complex(8)                                   :: Tr
  complex(8),allocatable                       :: dm_custom_rot(:,:)
  real(8)                                      :: soc
  !
  if (.not.allocated(single_particle_density_matrix)) then
     write(LOGfile,"(A)") "single_particle_density_matrix is not allocated"
     stop
  endif
  !
  if(allocated(dm_))                         deallocate(dm_)          ;allocate(dm_(Nspin*Norb,Nspin*Norb))          ;dm_ = zero
  ! if(allocated(dm_custom_rot))               deallocate(dm_custom_rot);allocate(dm_custom_rot(Nspin*Norb,Nspin*Norb));dm_custom_rot = zero
  ! if(present(dm_eig_).and.allocated(dm_eig_))deallocate(dm_eig_)      ;allocate(dm_eig_(Nspin*Norb))                 ;dm_eig_ = 0.0d0
  ! if(present(dm_rot_).and.allocated(dm_rot_))deallocate(dm_rot_)      ;allocate(dm_rot_(Nspin*Norb,Nspin*Norb))      ;dm_rot_ = zero
  !
  ! dm in the impurity problem basis
  dm_ = nn2so_reshape(single_particle_density_matrix,Nspin,Norb)
  !
  ! if(bath_type=="replica")then
  !    !
  !    ! dm in her diagonal basis
  !    if(present(dm_eig_).and.present(dm_rot_))then
  !       dm_rot_=dm_
  !       call eigh(dm_rot_,dm_eig_,'V','U')
  !    endif
  !    !
  !    ! dm in the basis defined by custom_rot
  !    if(present(custom_rot))then
  !       dm_custom_rot=matmul(transpose(conjg(custom_rot)),matmul(dm_,custom_rot))
  !    endif
  !    !
  ! elseif(bath_type=="normal")then !.and.SOC/=0.d0
  !    !
  !    ! here I assume that custom_rot is: {J,jz}-->{t2g,Sz} / {Lz,Sz}
  !    dm_custom_rot = nn2so_reshape(single_particle_density_matrix,Nspin,Norb)
  !    dm_=matmul(custom_rot,matmul(dm_custom_rot,transpose(conjg(custom_rot))))
  !    !
  !    ! dm in her diagonal basis
  !    if(present(dm_eig_).and.present(dm_rot_))then
  !       dm_rot_=dm_
  !       call eigh(dm_rot_,dm_eig_,'V','U')
  !    endif
  !    !
  ! endif
  !
  call print_dm(dm_,1)!,dm_rot_,dm_eig_,dm_custom_rot,1)
  !
end subroutine get_density_matrix_single


subroutine print_dm(dm_,ndx)!,dm_rot_,dm_eig_,dm_custom_rot,ndx)
  implicit none
  integer               ,intent(in)            :: ndx
  complex(8),allocatable,intent(in)            :: dm_(:,:)
  ! complex(8),allocatable,intent(in)            :: dm_custom_rot(:,:)
  ! real(8),allocatable   ,intent(in),optional   :: dm_eig_(:)
  ! complex(8),allocatable,intent(in),optional   :: dm_rot_(:,:)
  !internal
  integer                                      :: unit
  character(len=24)                            :: suffix
  integer                                      :: iorb,jorb,ispin,jspin,io,jo
  !
  suffix="single_particle_density_matrix_"//reg(str(ndx))//".dat"
  !
  unit = free_unit()
  open(unit,file=suffix,action="write",position="rewind",status='unknown')
  !
  write(unit,"(A90)")"# density matrix in the impurity problem basis REAL part:"
  do io=1,Nspin*Norb
     write(unit,"(90(F15.9,1X))") (dreal(dm_(io,jo)),jo=1,Nspin*Norb)
  enddo
  write(unit,*)
  !
  write(unit,"(A90)")"# density matrix in the impurity problem basis IMAGINARY part:"
  do io=1,Nspin*Norb
     write(unit,"(90(F15.9,1X))") (dimag(dm_(io,jo)),jo=1,Nspin*Norb)
  enddo
  write(unit,*)
  close(unit)
  !
end subroutine print_dm
