!ED_MAIN:
subroutine ed_read_input_c(instr) bind(c, name='read_input')
    use, intrinsic :: iso_c_binding
    character(kind=c_char), dimension(*), intent(IN) :: instr
    character(len=:), allocatable :: INPUTunit_fortran
    integer :: length
    integer :: i
    length=0
    do
       if (instr(length+1) == C_NULL_CHAR) exit
       length = length + 1
    end do
    if(allocated(INPUTunit_fortran))deallocate(INPUTunit_fortran)
    allocate(character(len=length) :: INPUTunit_fortran)
    do i = 1, length
      INPUTunit_fortran(i:i) = instr(i)
    enddo
    INPUTunit_fortran=trim(INPUTunit_fortran)
    call ed_read_input(INPUTunit_fortran) 
    deallocate(INPUTunit_fortran)
end subroutine ed_read_input_c
