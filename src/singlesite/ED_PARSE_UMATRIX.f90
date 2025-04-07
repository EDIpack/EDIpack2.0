MODULE ED_PARSE_UMATRIX
  !:synopsis: Routines and types for bath-impurity sparse maps
  USE SF_IOTOOLS, only: str,free_unit,file_length
  USE ED_VARS_GLOBAL
  USE ED_INPUT_VARS
  implicit none
  private
 
 
contains

  subroutine allocate_coulomb_matrix(length)
    !This subroutine allocates the list of two-body operators
    integer         :: length  !The number of two-body operators
    !
    if(coulombmatrix%status) stop "Coulomb matrix already allocated"
    allocate(coulombmatrix%oplist(length))
    coulombmatrix%status = .true.
  
  end subroutine allocate_coulomb_matrix
  
  subroutine deallocate_coulomb_matrix
    !This subroutine deallocates the list of two-body operators
    !
    if (.not. coulombmatrix%status) stop "Coulomb matrix not allocated"
    deallocate(coulombmatrix%oplist)
    coulombmatrix%status = .false.  
  end subroutine deallocate_coulomb_matrix

  
  subroutine read_umatrix_file(ufile)
    !This subroutine reads the interaction Hamiltonian from a user-specified file.
    character(len=*) :: ufile  !File containing a properly formatted interaction Hamiltonian
    logical          :: master=.true.,ufile_exists, verbose
    integer          :: iline,flen,unit_umatrix,rank
    integer          :: o1, o2, o3, o4
    character(len=1) :: s1, s2, s3, s4
    !
#ifdef _DEBUG
    if(ed_verbose>2)write(Logfile,"(A,A)")"DEBUG ed_read_input: read input from",trim(INPUTunit)
#endif
#ifdef _MPI    
    if(check_MPI())then
       master=get_Master_MPI()
       rank  =get_Rank_MPI()
    endif
#endif
    inquire(file=str(ufile),exist=ufile_exists)
    if(.not.ufile_exists)stop "read_dmft_bath ERROR: indicated file does not exist" !#FIXME: change this to make it default back to Uloc&co.
    if(ed_verbose>2)write(LOGfile,"(A)")'Reading interaction Hamiltonian from file '//str(ufile)
    if(ed_verbose>2)write(LOGfile,"(A)")'There are '//str(flen)//' two-body operators.'
    !
    flen = file_length(str(ufile),verbose=ed_verbose>2)
    !
    call  allocate_coulomb_matrix(flen)
    !
    open(free_unit(unit_umatrix),file=str(ufile))
    !
    !Take care of the comment in the preamble:
    read(unit_umatrix,*)
    !
    !Parse lines
    do iline=1,flen
       read(unit_umatrix,*) o1,s1,o2,s2,o3,s3,o4,s4,coulombmatrix%oplist(iline)%U
       if(max(o1, o2, o3, o4)>Norb) stop "read_umatrix_file: at line "//str(iline)//" too many orbitals" 
       if(min(o1, o2, o3, o4) <  1) stop "read_umatrix_file: at line "//str(iline)//" orbital index < 1" 
       if (any(( [s1, s2, s3, s4] /= "u" ) .AND. ( [s1, s2, s3, s4] /= "d" ))) stop "read_umatrix_file: at line "//str(iline)//" spin index malformed" 
       coulombmatrix%oplist(iline)%cd_i = [o1, merge(1, 2, s1 == "u")]
       coulombmatrix%oplist(iline)%cd_j = [o2, merge(1, 2, s2 == "u")]
       coulombmatrix%oplist(iline)%c_k  = [o3, merge(1, 2, s3 == "u")]
       coulombmatrix%oplist(iline)%c_l  = [o4, merge(1, 2, s4 == "u")]
       if(ed_verbose>4)write(LOGfile,"(A)")'Two-body operator '//str(iline)//' '//&
                                             str(coulombmatrix%oplist(iline)%U)//&
                                           ' cdg_['//str(o1)//str(s1)//']'&
                                           ' cdg_['//str(o2)//str(s2)//']'&
                                           ' cdg_['//str(o3)//str(s3)//']'&
                                           ' cdg_['//str(o4)//str(s4)//']'
    enddo
    !
    close(unit_umatrix)
    !
  end subroutine read_umatrix_file
  
  
  
  !subroutine normal_order
  
  
  !end subroutine normal_order
  
  
  !subroutine update_interaction_params
  
  
  !end subroutine update_interaction_params



end module ED_PARSE_UMATRIX






