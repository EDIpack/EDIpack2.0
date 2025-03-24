MODULE ED_VERSION
  implicit none
  !GIT VERSION
  character(len=41),parameter,public :: version = "d12b53fcc9ae32bbcf2ff299a786b3203e0ec144"
  contains			 !this is a silly trick to avoid the ar no symbol issues in OSX
  subroutine foo_ed_version_foo()
  end subroutine foo_ed_version_foo
END MODULE ED_VERSION
