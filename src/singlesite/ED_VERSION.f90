MODULE ED_VERSION
  implicit none
  !GIT VERSION
  character(len=41),parameter,public :: version = "ec5e02071231acaa1de09c591f76bf225f3be9e6"
  contains			 !this is a silly trick to avoid the ar no symbol issues in OSX
  subroutine foo_ed_version_foo()
  end subroutine foo_ed_version_foo
END MODULE ED_VERSION
