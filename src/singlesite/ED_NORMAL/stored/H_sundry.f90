  !We build the transposed H_non_local here (symmetric)
  !to comply with the MPI decomposition of the matrix.
  !A better MPI handling might be necessary here...
  do i=MpiIstart,MpiIend
     iup = iup_index(i,DimUp)
     idw = idw_index(i,DimUp)
     !
     mup = Hsector%H(1)%map(iup)
     mdw = Hsector%H(2)%map(idw)
     !
     nup = bdecomp(mup,Ns)
     ndw = bdecomp(mdw,Ns)
     !
     !
     if(allocated(coulomb_sundry))then
        do iline=1,size(coulomb_sundry)
           o1=coulomb_sundry(iline)%cd_i(1) 
           o2=coulomb_sundry(iline)%cd_j(1)
           o3=coulomb_sundry(iline)%c_k(1) 
           o4=coulomb_sundry(iline)%c_l(1)
           s1=coulomb_sundry(iline)%cd_i(2) 
           s2=coulomb_sundry(iline)%cd_j(2)
           s3=coulomb_sundry(iline)%c_k(2) 
           s4=coulomb_sundry(iline)%c_l(2) 
           
           !This needs to happen:
           !- this is normal, so no terms that unbalance the spin count can exist (e.g. no c_up c_up cd_dw cd_dw)
           !- we need to operate on either m_up or m_dw depending on the spins
           !- the operators need to be applied from right to left as c -> cd -> c -> cd to be consistent with the
           !  density terms originating from the anticommutation relations that are included in mfHloc
           
           !do jorb=1,Norb
           !   Jcondition=(&
           !        (iorb/=jorb).AND.&
           !        (nup(jorb)==1).AND.&
           !        (ndw(iorb)==1).AND.&
           !        (ndw(jorb)==0).AND.&
           !        (nup(iorb)==0))
           !   if(Jcondition)then
           !      call c(iorb,mdw,k1,sg1)  !DW
           !      call cdg(jorb,k1,k2,sg2) !DW
           !      jdw=binary_search(Hsector%H(2)%map,k2)
           !      call c(jorb,mup,k3,sg3)  !UP
           !      call cdg(iorb,k3,k4,sg4) !UP
           !      jup=binary_search(Hsector%H(1)%map,k4)
           !      htmp = Jx_internal(iorb,jorb)*sg1*sg2*sg3*sg4
           !      j = jup + (jdw-1)*DimUp
           !      !
           !      select case(MpiStatus)
           !      case (.true.)
           !         call sp_insert_element(MpiComm,spH0nd,htmp,i,j)
           !      case (.false.)
           !         call sp_insert_element(spH0nd,htmp,i,j)
           !      end select
           !      !
           !   endif
           !enddo
        enddo
     endif
     !
     !
  enddo
