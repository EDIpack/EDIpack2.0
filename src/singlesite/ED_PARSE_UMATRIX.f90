MODULE ED_PARSE_UMATRIX
  !:synopsis: Routines and types for bath-impurity sparse maps
  USE SF_IOTOOLS, only: str,free_unit,file_length,to_lower,txtfy
  USE ED_AUX_FUNX, only: print_hloc
  USE ED_VARS_GLOBAL
  USE ED_INPUT_VARS
  implicit none
  private
  
  public :: read_umatrix_file
 
contains


  
  subroutine read_umatrix_file(ufile)
    !This subroutine reads the interaction Hamiltonian from a user-specified file.
    !It would make sense for the elements to be ordered in this way:
    !cd_s cd_sprime c_s c_sprime.
    character(len=*)              :: ufile  !File containing a properly formatted interaction Hamiltonian
    character(len=300)            :: dummy
    type(coulomb_matrix_element)  :: opline
    logical                       :: master=.true.,ufile_exists, verbose, preamble
    integer                       :: iline,flen,unit_umatrix,rank, nops, iorb, jorb, ierr
    integer                       :: o1, o2, o3, o4
    character(len=1)              :: s1, s2, s3, s4
    !
#ifdef _MPI    
    if(check_MPI())then
       master=get_Master_MPI()
       rank  =get_Rank_MPI()
    endif
#endif
    inquire(file=trim(ufile)//reg(ed_file_suffix)//".restart",exist=ufile_exists)
    if(.not.ufile_exists)stop "read_umatrix_file ERROR: indicated file does not exist" !#FIXME: change this to make it default back to Uloc&co.
    if(ed_verbose>0)write(LOGfile,"(A)")'Reading interaction Hamiltonian from file '//trim(ufile)//reg(ed_file_suffix)//".restart"
    !
    !Set internal interaction coefficient matrices to zero
    mfHloc        = zero
    Uloc_internal = zero
    Ust_internal  = zero
    Jh_internal   = zero
    Jx_internal   = zero
    Jp_internal   = zero
    !
    !
    open(free_unit(unit_umatrix),file=trim(ufile)//reg(ed_file_suffix)//".restart")
    !
    !Take care of the comment in the preamble:
    preamble = .true.
    do while(preamble)
      read(unit_umatrix,*) dummy
      if(trim(dummy)/="#" .and. trim(dummy)/="!" .and. trim(dummy)/="%")then
        preamble = .false.
        read(dummy, *) o1
      endif
    enddo
    !
    if(o1/=Norb)STOP "Wrong number of orbitals in umatrix file header"
    !
    !Parse lines
    nops = 0
    do
         read(unit_umatrix, '(A)', iostat=ierr) dummy
         if (ierr /= 0) then
            exit  ! end-of-file error
         else
           read(dummy, *, iostat=ierr) o1, s1, o2, s2, o3, s3, o4, s4, opline%U
           if (ierr /= 0) then
               cycle  ! not a valid data line
           else
             nops = nops + 1
             if(max(o1, o2, o3, o4)>Norb) stop "read_umatrix_file: at line "//str(iline)//" too many orbitals" 
             if(min(o1, o2, o3, o4) <  1) stop "read_umatrix_file: at line "//str(iline)//" orbital index < 1" 
             if (any(( [s1, s2, s3, s4] /= "u" ) .AND. ( [s1, s2, s3, s4] /= "d" ))) stop "read_umatrix_file: at line "//str(iline)//" spin index malformed" 
             opline%cd_i = [o1, merge(1, 2, s1 == "u")]
             opline%cd_j = [o2, merge(1, 2, s2 == "u")]
             opline%c_k  = [o3, merge(1, 2, s3 == "u")]
             opline%c_l  = [o4, merge(1, 2, s4 == "u")]
             if(ed_verbose>3)then
               write(dummy, '(F10.6)') opline%U 
               write(LOGfile,"(A)")'Two-body operator '//txtfy(nops,3)//':     '//&
                                                        trim(dummy)//&
                                                      ' cd_['//str(o1)//str(s1)//']'//&
                                                      ' cd_['//str(o2)//str(s2)//']'//&
                                                      ' c_['//str(o3)//str(s3)//']'//&
                                                      ' c_['//str(o4)//str(s4)//']'
           endif
         endif
         call parse_umatrix_line(opline)
        endif
    enddo
    if(ed_verbose>2)write(LOGfile,"(A)")str(nops)//" two-body operators parsed."
    !
    close(unit_umatrix)
    !    
    if(ed_verbose>3)then
      write(LOGfile,"(A)")''
      write(LOGfile,"(A)")'Interaction coefficients:'
      write(LOGfile,"(A)")''
    endif
    !
    !Here we need to operate on the various Uloc, Ust, Jh, Jx, Jp matrices
    !-the Hubbard terms are passed with an 1/2, but if the user made things
    !correctly they already have the two nup/ndw and ndw/nup terms. So nothing to do here.
    if(ed_verbose>3)then
      write(LOGfile,"(A)")'ULOC:'
      write(LOGfile,"(90(F15.9,1X))") (Uloc_internal(iorb),iorb=1,Norb)
      write(LOGfile,"(A)")''
    endif
    !
    !-Ust is symmetric with respect to orbital exchange: the elements we got should be ordered
    !in increasing orbital order. So the matrix should be upper triangular. In the routine that
    !creates Hloc, it does nup*ndw + ndw*nup, so here we divide by 2.0
    Ust_internal = (Ust_internal + transpose(Ust_internal))/2.0
    if(ed_verbose>3)then
      write(LOGfile,"(A)")'UST:'
      do iorb=1,Norb
        write(LOGfile,"(90(F15.9,1X))") (Ust_internal(iorb,jorb),jorb=1,Norb)
      enddo
      write(LOGfile,"(A)")''
    endif

    !-Jh needs to be rescaled. First, it also is symmetric w.r.t. orbital exchange.
    !Then, the coefficient of the terms in the H constructor is actually Ust-Jh. 
    !So, if the user passed this as Jh, we need to recast it as Ust - what the user passed
    Jh_internal = (Jh_internal + transpose(Jh_internal))/2.0
    Jh_internal = Ust_internal - Jh_internal
    if(ed_verbose>3)then
      write(LOGfile,"(A)")'JH:'
      do iorb=1,Norb
        write(LOGfile,"(90(F15.9,1X))") (Jh_internal(iorb,jorb),jorb=1,Norb)
      enddo
      write(LOGfile,"(A)")''
    endif
    
    !Jx and Jp have a summation that goes from 1 to Norb for both orbital indices, so no change there
    if(ed_verbose>3)then
      write(LOGfile,"(A)")'JX:'
      do iorb=1,Norb
        write(LOGfile,"(90(F15.9,1X))") (Jx_internal(iorb,jorb),jorb=1,Norb)
      enddo
      write(LOGfile,"(A)")'JP:'
      do iorb=1,Norb
        write(LOGfile,"(90(F15.9,1X))") (Jp_internal(iorb,jorb),jorb=1,Norb)
      enddo
      write(LOGfile,"(A)")''
    endif
    
     !Print mean-field terms
     if(ed_verbose>3)then
      write(LOGfile,"(A)")'Mean-field terms from anticommutators:'
      call print_hloc(mfHloc)
     endif
    
    !Is there anything else?
    if(ed_verbose>3 .and. allocated(coulomb_sundry))then
      write(LOGfile,"(A)")'There are '//str(size(coulomb_sundry))//' sundry terms.'
      do iline=1,size(coulomb_sundry)
         o1 = coulomb_sundry(iline)%cd_i(1)
         o2 = coulomb_sundry(iline)%cd_j(1)
         o3 = coulomb_sundry(iline)%c_k(1)
         o4 = coulomb_sundry(iline)%c_l(1)
         s1 = merge("u", "d", coulomb_sundry(iline)%cd_i(2)==1)
         s1 = merge("u", "d", coulomb_sundry(iline)%cd_j(2)==1)
         s1 = merge("u", "d", coulomb_sundry(iline)%c_k(2)==1)
         s1 = merge("u", "d", coulomb_sundry(iline)%c_l(2)==1)
         write(dummy, '(F10.6)') coulomb_sundry(iline)%U
         write(LOGfile,"(A)")'Sundry operator '//txtfy(iline,3)//':     '//&
                                            trim(dummy)//&
                                          ' cd_['//str(o1)//str(s1)//']'//&
                                          ' c_['//str(o2)//str(s2)//']'//&
                                          ' cd_['//str(o3)//str(s3)//']'//&
                                          ' c_['//str(o4)//str(s4)//']'
      enddo
      write(LOGfile,"(A)")''
    endif
    !
    !
  end subroutine read_umatrix_file
  
  
  subroutine parse_umatrix_line(line)
    !This subroutine switches the second and third elements of the coulomb matrix applying
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
    !Zero: onsistently with w2dynamics, 1/2 prefactor is applied by the code.
    !It is also multiplied by -1, because the line in the w2dynamics file is the following:
    ![i j k l U_ijkl], but the operator is defined as cd_i cd_j U_ijkl c_l c_k.
    !Note that l and k are inverted.
    line%U = -0.5d0 * line%U
    !
    !First: order the two creation operators so that they
    !are set in an increasing order of (first) spin and (second) orbital from left to right
    !
    !orbital
    if (line%cd_i(1) > line%cd_j(1)) then
      dummy = line%cd_i
     line%cd_i = line%cd_j
      line%cd_j = dummy
      line%U    = -1.0*line%U
    endif
    !spin (overrides orbital)
    if (line%cd_i(2) > line%cd_j(2)) then
      dummy = line%cd_i
      line%cd_i = line%cd_j
      line%cd_j = dummy
      line%U    = -1.0*line%U
    endif
    !
    !Second: order the two annihilation operators so that they
    !are set in an increasing order of (first) spin and (second) orbital from left to right
    !
    !orbital
    if (line%c_k(1) > line%c_l(1)) then
      dummy = line%c_k
      line%c_k = line%c_l
      line%c_l = dummy
      line%U    = -1.0*line%U
    endif
    !spin (overrides orbital)
    if (line%c_k(2) > line%c_l(2)) then
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
    !S-E: -J c^+_a_up c_a_dw c^+_b_dw c_b_up (i.ne.j)
    !Note that the sign change was already done at step 0
    !
    if(line%cd_i(1) /= line%cd_j(1) .and.& !iorb != jorb
       line%cd_i(1) /= line%c_k(1)  .and.& 
       line%cd_i(2) == line%c_k(2)  .and.&
       line%cd_j(1) /= line%c_l(1)  .and.&
       line%cd_j(2) == line%c_l(2))  then         
       Jx_internal(line%cd_i(1),line%c_k(1)) = Jx_internal(line%cd_i(1),line%c_k(1)) + line%U
       return
    endif
    !
    !Seventh: is it pair-hopping?
    !P-H: -J c^+_iorb_up c_jorb_up c^+_iorb_dw  c_jorb_dw (i.ne.j)
    !Note that the sign change was already done at step 0
    !
    if(line%cd_i(1) /= line%c_k(1)  .and.& !iorb != jorb
       line%cd_i(1) == line%cd_j(1) .and.& 
       line%cd_i(2) /= line%cd_j(2) .and.&
       line%c_k(1)  == line%c_l(1)  .and.&
       line%c_k(2)  /= line%c_l(2))  then          
       Jp_internal(line%cd_i(1),line%c_k(1)) = Jp_internal(line%cd_i(1),line%c_k(1)) + line%U
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






