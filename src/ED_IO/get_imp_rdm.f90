
subroutine get_rdm_single(rdm,doprint)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:),intent(inout) :: rdm
  logical,intent(in),optional             :: doprint
  logical                                 :: doprint_
  !
  doprint_=.false.; if(present(doprint)) doprint_=doprint
  !
  if(.not.allocated(impurity_density_matrix))stop "ERROR: impurity_density_matrix is not allocated"
  !
  ! if(allocated(rdm))deallocate(rdm)
  ! allocate(rdm, source=impurity_density_matrix)
  call assert_shape(rdm,shape(impurity_density_matrix),"get_rdm_single","rdm")
  rdm = impurity_density_matrix
  !
  if(doprint_)call print_rdm(rdm,4**Norb)
  !
contains
  !
  subroutine print_rdm(dm,Nrdm)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
    integer                  ,intent(in)            :: Nrdm
    complex(8),dimension(:,:),intent(in)            :: dm
    integer                                         :: unit
    character(len=64)                               :: fname
    integer                                         :: io,jo
    !
    if(size(dm,1)/=Nrdm.OR.size(dm,2)/=Nrdm)then
       stop "ERROR: actual dm argument has incogruent size wrt explicitly passed Nrdm"
    endif
    !
    !
    fname = "reduced_density_matrix"//str(ed_file_suffix)//".ed"
    !
    unit = free_unit()
    open(unit,file=fname,action="write",position="rewind",status='unknown')
    do io=1,Nrdm
       write(unit,"(*(F20.16,1X))") (dreal(dm(io,jo)),jo=1,Nrdm)
    enddo
    if(any(dimag(dm)/=0d0))then
       write(unit,*)
       do io=1,Nrdm
          write(unit,"(*(F20.16,1X))") (dimag(dm(io,jo)),jo=1,Nrdm)
       enddo
    endif
    close(unit)
    !
  end subroutine print_rdm
  !
end subroutine get_rdm_single

