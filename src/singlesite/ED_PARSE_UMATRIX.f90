MODULE ED_PARSE_UMATRIX
  !:synopsis: Routines and types for bath-impurity sparse maps
  USE SF_IOTOOLS, only: str,free_unit,file_length
  USE ED_VARS_GLOBAL
  USE ED_INPUT_VARS
  implicit none
  private
 
 
contains


  
  subroutine read_umatrix_file(ufile)
    !This subroutine reads the interaction Hamiltonian from a user-specified file.
    !It would make sense for the elements to be ordered in this way:
    !cd_s cd_sprime c_s c_sprime.
    character(len=*)              :: ufile  !File containing a properly formatted interaction Hamiltonian
    type(coulomb_matrix_element)  :: opline
    logical                       :: master=.true.,ufile_exists, verbose
    integer                       :: iline,flen,unit_umatrix,rank
    integer                       :: o1, o2, o3, o4
    character(len=1)              :: s1, s2, s3, s4
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
    !Set internal interaction coefficient matrices to zero
    mfHloc        = zero
    Uloc_internal = zero
    Ust_internal  = zero
    Jh_internal   = zero
    Jx_internal   = zero
    Jp_internal   = zero
    !
    flen = file_length(str(ufile),verbose=ed_verbose>2)
    !
    !
    open(free_unit(unit_umatrix),file=str(ufile))
    !
    !Take care of the comment in the preamble:
    read(unit_umatrix,*)
    !
    !Parse lines
    do iline=1,flen
       read(unit_umatrix,*) o1, s1, o2, s2, o3, s3, o4, s4, opline%U
       if(max(o1, o2, o3, o4)>Norb) stop "read_umatrix_file: at line "//str(iline)//" too many orbitals" 
       if(min(o1, o2, o3, o4) <  1) stop "read_umatrix_file: at line "//str(iline)//" orbital index < 1" 
       if (any(( [s1, s2, s3, s4] /= "u" ) .AND. ( [s1, s2, s3, s4] /= "d" ))) stop "read_umatrix_file: at line "//str(iline)//" spin index malformed" 
       opline%cd_i = [o1, merge(1, 2, s1 == "u")]
       opline%cd_j = [o2, merge(1, 2, s2 == "u")]
       opline%c_k  = [o3, merge(1, 2, s3 == "u")]
       opline%c_l  = [o4, merge(1, 2, s4 == "u")]
       if(ed_verbose>4)write(LOGfile,"(A)")'Two-body operator '//str(iline)//' '//&
                                             str(opline%U)//&
                                           ' cd_['//str(o1)//str(s1)//']'&
                                           ' cd_['//str(o2)//str(s2)//']'&
                                           ' c_['//str(o3)//str(s3)//']'&
                                           ' c_['//str(o4)//str(s4)//']'
      call parse_umatrix_line(opline)
    enddo
    !
    !Here we need to operate on the various Uloc, Ust, Jh, Jx, Jp matrices
    !-Uloc needs the 1/2 to be dealt with
    !-Ust, Jh, Jx, Jp have to keep in mind that the summations are constrained jorb > iorb
    !-Jh needs to be rescaled, because the user has inputted Ust - Jh in the umatrix file
    !-Remember to add mfHloc to the hamiltonian constructor!
    !
    close(unit_umatrix)
    !
  end subroutine read_umatrix_file
  
  
  subroutine parse_umatrix_line(line)
    !This funcion switches the second and third elements of the coulomb matrix applying
    !commutation relations. If the resulting series of operators is of the type that would
    !couple to Uloc, Ust, Ust-Jh, Jx, Jp, store the resulting coefficient in the appropriate
    !matrix .
    !If some mean-field terms result from the anticommutation, store them in a mean-field 
    !matrix to be added along with impHloc at Fock space H creation time.
    !If a Coulomb matrix element comes along that would be inconsistent with the previously
    !determined value of Uloc, Ust, Ust-Jh, Jx, Jp then set the previously determined coeff
    !to zero and add this to the list of "everything_else".
    type(coulomb_matrix_element)        :: line
    character(len=10)                   :: operator_type
    integer,dimension(2)                :: dummy
    !
    !First: order the two creation operators so that they
    !are set in an increasing order of (first) spin and (second) orbital from left to right
    !
    !spin
    if (line%cd_i(2) > line%cd_j(2)) then
      dummy = line%cd_i
      line%cd_i = line%cd_j
      line%cd_j = dummy
      line%U    = -1.0*line%U
    !orbital
    elseif (line%cd_i(1) > line%cd_j(1)) then
      dummy = line%cd_i
      line%cd_i = line%cd_j
      line%cd_j = dummy
      line%U    = -1.0*line%U
    endif
    !
    !Second: order the two annihilation operators so that they
    !are set in an increasing order of (first) spin and (second) orbital from left to right
    !
    !spin
    if (line%c_k(2) > line%c_l(2)) then
      dummy = line%c_k
      line%c_k = line%c_l
      line%c_l = dummy
      line%U    = -1.0*line%U
    !orbital
    elseif (line%c_k(1) > line%c_l(1)) then
      dummy = line%c_k
      line%c_k = line%c_l
      line%c_l = dummy
      line%U    = -1.0*line%U
    endif
    !
    !Third: are the indices of the swapped operators the same? If so, there is a mean-field
    !term coming out of the anticommutator
    !
    if(all((line%cd_j==line%c_k)))then
      mfHloc(line%cd_i(2),line%c_k(2),line%cd_i(1),line%c_k(1)) = mfHloc(line%cd_i(2),line%c_k(2),line%cd_i(1),line%c_k(1)) + line%U
    endif  
    !
    !Fourth: multiply U by (-1) because operators will always be applied from right to left
    !as c->cd->c->cd, so second and third element are to be swapped
    !
    line%U = -1.0 * line%U
    !
    !Fifth: look at the four-operator term. Does it look like any of the ones we would put 
    !in an Hubbard-Kanamori density-density interaction?
    !
    if(line%cd_i(1)==line%c_k(1))then          !c_dag c: First two operators form a density
      if(line%cd_j(1)==line%c_l(1))then        !c_dag c: Second two operators form a density
        if(line%cd_i(2)/=line%cd_j(2))then          !c_dag spin are opposite: It is either Uloc or Ust
          if(line%cd_i(1)==line%cd_j(1))then        !It is Uloc
            Uloc_internal(line%cd_i(1)) = Uloc_internal(line%cd_i(1)) + line%U
            return
          else    !It is Ust
            Ust_internal(line%cd_i(1),line%cd_j(1)) = Ust_internal(line%cd_i(1),line%cd_j(1)) + line%U
            return
          endif
        else
          if(line%cd_i(1)/=line%cd_j(1))then        !It is Ust - Jh. Store this coefficient in Jh for now. Will need to do Jh = Ust - this in the end.
            Jh_internal(line%cd_i(1),line%cd_j(1)) = Jh_internal(line%cd_i(1),line%cd_j(1)) + line%U
            return
          endif
        endif
      endif
    endif
    !
    !Sixth: is it spin-exchange?
    !S-E: -J c^+_a_up c^+_b_dw c_b_up c_a_dw (i.ne.j)
    !
    if(line%cd_i(1) /= line%cd_j(1) .and.& !iorb != jorb
       line%cd_i(1) /= line%c_k(1)  .and.& 
       line%cd_i(2) == line%c_k(2)  .and.&
       line%cd_j(1) /= line%c_l(1)  .and.&
       line%cd_j(2) == line%c_l(2))  then         
       Jx_internal(line%cd_i(1),line%cd_j(1)) = Jx_internal(line%cd_i(1),line%cd_j(1)) - line%U
       return
    endif
    !
    !Seventh: is it pair-hopping?
    !P-H: -J c^+_iorb_up c^+_iorb_dw   c_jorb_up  c_jorb_dw (i.ne.j)
    !
    if(line%cd_i(1) /= line%c_k(1)  .and.& !iorb != jorb
       line%cd_i(1) == line%cd_j(1) .and.& 
       line%cd_i(2) /= line%cd_j(2) .and.&
       line%c_k(1)  == line%c_l(1)  .and.&
       line%c_k(2)  /= line%c_l(2))  then          
       Jx_internal(line%cd_i(1),line%cd_j(1)) = Jp_internal(line%cd_i(1),line%cd_j(1)) - line%U
       return
    endif
    !
    !Eight: if it is none of the above, put this into coulomb_everything_else
    !
    call grow_sundry_array(line)
      
    
  end subroutine parse_umatrix_line
  

  subroutine grow_sundry_array(new_element)
    type(coulomb_matrix_element),intent(in)                 :: new_element
		type(coulomb_matrix_element),dimension(:),allocatable   :: temp
		integer                                                 :: dim_old
    !
    if(.not.allocated(coulomb_sundry))then
      allocate(coulomb_sundry(1))
      coulomb_sundry(1) = new_element
    else     
      dim_old = size(coulomb_sundry)
      allocate(temp(dim_old+1))
		  temp(1:dim_old) = coulomb_sundry
      temp(dim_old+1) = new_element
		  call move_alloc(temp,coulomb_sundry)
		endif
    !
  end subroutine grow_sundry_array    


end module ED_PARSE_UMATRIX






