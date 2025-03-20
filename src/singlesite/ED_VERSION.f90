MODULE ED_VERSION
  implicit none
  !GIT VERSION
  character(len=41),parameter,public :: version = "9d2699879630e789092bc6db917d4c89556f7ee4"
  contains			 !this is a silly trick to avoid the ar no symbol issues in OSX
  subroutine foo_ed_version_foo()
  end subroutine foo_ed_version_foo
END MODULE ED_VERSION
