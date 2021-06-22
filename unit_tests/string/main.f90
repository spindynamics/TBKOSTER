program main
	print *
	!call test_clean_str()
	!print *
	call test_int2str()
	print *
	call test_is_str_int()
	print *
	call test_lower()
	print *
	call test_real2str()
	print *
	call test_str2int()
	print *
	call test_unique_str()
	print *
	call test_upper()
	print *
end program main

subroutine test_clean_str()
	use string_mod
	implicit none

end subroutine test_clean_str

subroutine test_int2str()
	use string_mod
	implicit none

	print *, "test_int2str()"

	print *, "int2str(-123) ="
	print *, int2str(-123)
end subroutine test_int2str

subroutine test_is_str_int()
	use string_mod
	implicit none

	print *, "test_is_str_int()"

	print *, "is_str_int('-123') ="
	print *, is_str_int('-123')

	print *, "is_str_int('123a') ="
	print *, is_str_int('123a')
end subroutine test_is_str_int

subroutine test_lower()
	use string_mod
	implicit none

	print *, "test_lower()"
	print *, "lower('Hello') ="
	print *, "'", lower('Hello'), "'"
end subroutine test_lower

subroutine test_real2str()
	use string_mod
	implicit none

	print *, "test_real2str()"

	print *, "real2str(-12.3) ="
	print *, real2str(-12.3_8)
end subroutine test_real2str

subroutine test_str2int()
	use string_mod
	implicit none

	print *, "test_str2int()"

	print *, "str2int(-123) ="
	print *, str2int('-123')

	!print *, "str2int(123a) ="
	!print *, str2int('123a')
end subroutine test_str2int

subroutine test_unique_str()
	use string_mod
	implicit none

	!character(len=2),dimension(*),parameter :: a = (/'H ','He','He','Li','Li','Li'/)
	character(len=2),dimension(*),parameter :: a = (/'Li','Li','Li','He','He','H '/)
	character(len=2),dimension(:),allocatable :: c
	integer,dimension(:),allocatable :: ia
	integer,dimension(size(a)) :: ic

	print *, "test_unique_str()"
	!print *, "unique_str((/'H ','He','He','Li','Li','Li'/),c,ia,ic)"
	print *, "unique_str((/'Li','Li','Li','He','He','H '/),c,ia,ic)"
	call unique_str(a,c,ia,ic)
	print *, "c = ", c
	print *, "ia = ", ia
	print *, "ic = ", ic
end subroutine test_unique_str

subroutine test_upper()
	use string_mod
	implicit none

	print *, "test_upper()"
	print *, "upper('Hello') ="
	print *, "'", upper('Hello'), "'"
end subroutine test_upper

