!##################################################################
!LATTICE EXTENSION: need to read to the 
!##################################################################


subroutine ed_get_dimp_lattice_n2(self,nlat,axis,z)
#if __INTEL_COMPILER
    use ED_INPUT_VARS, only: Nspin,Norb
#endif
  complex(8),dimension(:,:),intent(inout)       :: self !! [Nlso,Nlso,:]
  integer,intent(in)                            :: nlat  ! Number of inequivalent impurity sites for real-space DMFT
  character(len=*),optional                     :: axis ! Can be :f:var:`"m"` for Matsubara (default), :f:var:`"r"` for real
  complex(8),dimension(:),optional              :: z    ! User provided array of complex frequency where to evaluate Self
  character(len=1)                              :: axis_
  complex(8),dimension(:),allocatable           :: z_
  integer                                       :: ilat
  !
  axis_='m';if(present(axis))axis_=trim(axis)
  !
  call allocate_grids
  !
  if(present(z))then
     allocate(z_, source=z)
  else
     select case(axis_)
     case default;stop "ed_get_dimp ERROR: axis is neither Matsubara, nor Realaxis"
     case ('m','M')
        z_ = dcmplx(0d0,wm)
     case ('r','R')
        z_ = dcmplx(wr,eps)
     end select
  endif
  !
  L = size(z_)
  !
  call assert_shape(self,[Nlat,L],'ed_get_dimp','self')
  !
  do ilat=1,Nlat
     call ed_set_suffix(ilat)
     call ed_read_impDmatrix()
     self(ilat,:) = ed_build_impD(z_,axis_)
  enddo
  !
  call ed_reset_suffix()
  call deallocate_grids()
  ! call deallocate_GFmatrix(impDmatrix)
  !
end subroutine ed_get_dimp_lattice_n2

