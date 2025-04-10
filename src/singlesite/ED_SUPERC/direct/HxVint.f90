  !density-density interaction: same orbital, opposite spins:
  ! = \sum_\a U_\a*(n_{\a,up}*n_{\a,dw})
  htmp = zero
  do iorb=1,Norb
     htmp = htmp + Uloc_internal(iorb)*nup(iorb)*ndw(iorb)
  enddo
  if(Norb>1)then
     !density-density interaction: different orbitals, opposite spins:
     ! =   U'   *     sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
     ! =  (Uloc_internal-2*Jh_internal(iorb,jorb))*sum_{i/=j} [ n_{i,up}*n_{j,dw} + n_{j,up}*n_{i,dw} ]
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           htmp = htmp + Ust_internal(iorb,jorb)*(nup(iorb)*ndw(jorb) + nup(jorb)*ndw(iorb))
        enddo
     enddo
     !density-density interaction: different orbitals, parallel spins
     ! = \sum_{i<j}    U''     *[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
     ! = \sum_{i<j} (Uloc_internal-3*Jh_internal(iorb,jorb))*[ n_{i,up}*n_{j,up} + n_{i,dw}*n_{j,dw} ]
     do iorb=1,Norb
        do jorb=iorb+1,Norb
           htmp = htmp + (Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*(nup(iorb)*nup(jorb) + ndw(iorb)*ndw(jorb))
        enddo
     enddo
  endif
  !if using the Hartree-shifted chemical potential: mu=0 for half-filling
  !sum up the contributions of hartree terms:
  if(hfmode)then
     do iorb=1,Norb
        htmp = htmp - 0.5d0*Uloc_internal(iorb)*(nup(iorb)+ndw(iorb)) + 0.25d0*Uloc_internal(iorb)
     enddo
     if(Norb>1)then
        do iorb=1,Norb
           do jorb=iorb+1,Norb
              htmp=htmp-0.5d0*Ust_internal(iorb,jorb)*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.5d0*Ust_internal(iorb,jorb)
              htmp=htmp-0.5d0*(Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))*(nup(iorb)+ndw(iorb)+nup(jorb)+ndw(jorb))+0.5d0*(Ust_internal(iorb,jorb)-Jh_internal(iorb,jorb))
           enddo
        enddo
     endif
  endif
  !
  j = i
  hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)
  !


  ! SPIN-EXCHANGE (S-E) and PAIR-HOPPING TERMS
  !    S-E: J c^+_iorb_up c^+_jorb_dw c_iorb_dw c_jorb_up  (i.ne.j) 
  !    S-E: J c^+_{iorb} c^+_{jorb+Ns} c_{iorb+Ns} c_{jorb}
  if(Norb>1.AND.(Jx_internal(iorb,jorb)/=0d0.OR.Jp_internal(iorb,jorb)/=0d0))then
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition=(&
                (iorb/=jorb).AND.&
                (ib(jorb)==1).AND.&
                (ib(iorb+Ns)==1).AND.&
                (ib(jorb+Ns)==0).AND.&
                (ib(iorb)==0))
           if(Jcondition)then
              call c(jorb,m,k1,sg1)
              call c(iorb+Ns,k1,k2,sg2)
              call cdg(jorb+Ns,k2,k3,sg3)
              call cdg(iorb,k3,k4,sg4)
              j_el = binary_search(Hsector%H(1)%map,k2)
              j    = j_el + (iph-1)*DimEl
              htmp = one*Jx_internal(iorb,jorb)*sg1*sg2*sg3*sg4
              !
              if(j/=0)hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)
              !
           endif
        enddo
     enddo
  endif
  !
  ! PAIR-HOPPING (P-H) TERMS
  !    P-H: J c^+_iorb_up c^+_iorb_dw   c_jorb_dw   c_jorb_up  (i.ne.j) 
  !    P-H: J c^+_{iorb}  c^+_{iorb+Ns} c_{jorb+Ns} c_{jorb}
  if(Norb>1.AND.(Jx_internal(iorb,jorb)/=0d0.OR.Jp_internal(iorb,jorb)/=0d0))then
     do iorb=1,Norb
        do jorb=1,Norb
           Jcondition=(&
                (iorb/=jorb).AND.&
                (ib(jorb)==1).AND.&
                (ib(jorb+Ns)==1).AND.&
                (ib(iorb+Ns)==0).AND.&
                (ib(iorb)==0))
           if(Jcondition)then
              call c(jorb,m,k1,sg1)
              call c(jorb+Ns,k1,k2,sg2)
              call cdg(iorb+Ns,k2,k3,sg3)
              call cdg(iorb,k3,k4,sg4)
              j_el = binary_search(Hsector%H(1)%map,k2)
              j    = j_el + (iph-1)*DimEl
              htmp = one*Jp_internal(iorb,jorb)*sg1*sg2*sg3*sg4
              !
              if(j/=0)hv(j-MpiIshift) = hv(j-MpiIshift) + htmp*vin(i)
              !
           endif
        enddo
     enddo
  endif
