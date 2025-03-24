MODULE ED_VERSION
  implicit none
  !GIT VERSION
  character(len=41),parameter,public :: version = "17a7eee9950127e69a661ce6170479df1c018599"
  contains			 !this is a silly trick to avoid the ar no symbol issues in OSX
  subroutine foo_ed_version_foo()
  end subroutine foo_ed_version_foo
END MODULE ED_VERSION
