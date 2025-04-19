  !Diagonal Elements, i.e. local part
  htmp = zero
  htmp = htmp - xmu*(sum(nup)+sum(ndw))
  !
  do iorb=1,Norb
    htmp = htmp + (impHloc(1,1,iorb,iorb) + mfHloc(1,1,iorb,iorb))*nup(iorb)
    htmp = htmp + impHloc(Nspin,Nspin,iorb,iorb)*ndw(iorb)
    htmp = htmp + mfHloc(2,2,iorb,iorb) * ndw(iorb)
  enddo
  !
  hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(i)
  !
  !Off-diagonal elements, i.e. non-local part
  !1. same spin:
  do iorb=1,Norb
     do jorb=1,Norb
        !UP
       Jcondition = &
            ((impHloc(1,1,iorb,jorb)/=zero) .OR. &
             (mfHloc(1,1,iorb,jorb)/=zero)).AND. &
            (ib(jorb)==1)                  .AND. &
            (ib(iorb)==0)
        if (Jcondition) then
           call c(jorb,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp = conjg(impHloc(1,1,iorb,jorb)+mfHloc(1,1,iorb,jorb))*sg1*sg2
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
        !DW
       Jcondition = &
            ((impHloc(Nspin,Nspin,iorb,jorb)/=zero) .OR. &
             (mfHloc(2,2,iorb,jorb)/=zero)).AND. &
            (ib(jorb+Ns)==1)                       .AND. &
            (ib(iorb+Ns)==0)
        if (Jcondition) then
           call c(jorb+Ns,m,k1,sg1)
           call cdg(iorb+Ns,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp = conjg(impHloc(Nspin,Nspin,iorb,jorb))*sg1*sg2
           htmp = htmp + conjg(mfHloc(2,2,iorb,jorb))*sg1*sg2
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
     enddo
  enddo
  !
  !Evaluate: Fd . D = Fd . (C^+_{a,up}C^+_{a,dw} + C_{a,dw}C_{a,up})
  if(any(pair_field/=0d0))then
     do iorb=1,Norb
        !
        Jcondition = (ib(iorb)==1) .AND. (ib(iorb+Ns)==1)
        if(Jcondition)then
           call c(iorb,m,k1,sg1)
           call c(iorb+Ns,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp=one*pair_field(iorb)*sg1*sg2
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
        !
        Jcondition = (ib(iorb)==0) .AND. (ib(iorb+Ns)==0)
        if(Jcondition)then
           call cdg(iorb+Ns,m,k1,sg1)
           call cdg(iorb,k1,k2,sg2)
           j_el = binary_search(Hsector%H(1)%map,k2)
           j    = j_el + (iph-1)*DimEl
           htmp=one*pair_field(iorb)*sg1*sg2 !
           !
           hv(i-MpiIshift) = hv(i-MpiIshift) + htmp*vin(j)
           !
        endif
     enddo
  endif
