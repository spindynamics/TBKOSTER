program main
	print *
	call test_zeros1()
	print *
	call test_zeros2()
	print *
	call test_ones1()
	print *
	call test_ones2()
	print *
	call test_diag1()
	print *
	call test_diag2()
	print *
	call test_eye()
	print *
	call test_cross_product()
	print *
	call test_determinant()
	print *
	call test_inverse()
	print *
	call test_unique_int()
	print *
end program main

subroutine test_zeros1()
	use math_mod
	implicit none

	print *, "test_zeros1()"
	print *, "zeros1(3) ="
	print *, zeros1(3)
end subroutine test_zeros1

subroutine test_zeros2()
	use math_mod
	use precision_mod
	implicit none

	real(rp), dimension(3,3) :: m

	print *, "test_zeros2()"
	print *, "zeros2(3,3) ="
	m = zeros2(3,3)
	print *, m(1,:)
	print *, m(2,:)
	print *, m(3,:)
end subroutine test_zeros2

subroutine test_ones1()
	use math_mod
	implicit none

	print *, "test_ones1()"
	print *, "ones1(3) ="
	print *, ones1(3)
end subroutine test_ones1

subroutine test_ones2()
	use math_mod
	use precision_mod
	implicit none

	real(rp), dimension(3,3) :: m

	print *, "test_ones2()"	
	print *, "ones2(3,3) ="
	m = ones2(3,3)
	print *, m(1,:)
	print *, m(2,:)
	print *, m(3,:)
end subroutine test_ones2

subroutine test_diag1()
	use math_mod
	use precision_mod
	implicit none

	real(rp), dimension(3,3) :: m

	print *, "test_diag1()"
	print *, "diag1((/1.0,2.0,3.0/)) ="
	m = diag1((/1.0_rp,2.0_rp,3.0_rp/))
	print *, m(1,:)
	print *, m(2,:)
	print *, m(3,:)
end subroutine test_diag1

subroutine test_diag2()
	use math_mod
	use precision_mod
	implicit none

	print *, "test_diag2()"
	print *, "diag2((/1.0,0.0,0.0; 0.0,2.0,0.0; 0.0,0.0,3.0/)) ="
	print *, diag2(reshape(&
	 (/1.0_rp,0.0_rp,0.0_rp, &
	   0.0_rp,2.0_rp,0.0_rp, &
	   0.0_rp,0.0_rp,3.0_rp/), &
	 (/3,3/)))
end subroutine test_diag2

subroutine test_eye()
	use math_mod
	use precision_mod
	implicit none

	real(rp), dimension(3,3) :: m

	print *, "test_eye()"
	print *, "eye(3) ="
	m = eye(3)
	print *, m(1,:)
	print *, m(2,:)
	print *, m(3,:)
end subroutine test_eye

subroutine test_cross_product()
	use math_mod
	use precision_mod
	implicit none

	print *, "test_cross_product()"
	print *, "cross_product((/1.0,2.0,3.0/),(/2.0,3.0,4.0/)) ="
	print *, cross_product((/1.0_rp,2.0_rp,3.0_rp/),(/2.0_rp,3.0_rp,4.0_rp/))
end subroutine test_cross_product

subroutine test_determinant()
	use math_mod
	use precision_mod
	implicit none

	print *, "test_determinant()"
	print *, "determinant((/1.0,0.0,0.0; 0.0,2.0,0.0; 0.0,0.0,3.0/)) ="
	print *, determinant(reshape(&
	 (/1.0_rp,0.0_rp,0.0_rp, &
	   0.0_rp,2.0_rp,0.0_rp, &
	   0.0_rp,0.0_rp,3.0_rp/), &
	 (/3,3/)))
end subroutine test_determinant

subroutine test_inverse()
	use math_mod
	use precision_mod
	implicit none

	real(rp), dimension(3,3) :: m
	
	print *, "test_inverse()"
	print *, "inverse(diag1((/1.0,2.0,3.0/))) ="
	m = inverse(diag1((/1.0_rp,2.0_rp,3.0_rp/)))
	print *, m(1,:)
	print *, m(2,:)
	print *, m(3,:)
end subroutine test_inverse

subroutine test_unique_int()
	use math_mod
	implicit none

	!integer,dimension(*),parameter :: a = (/1,2,2,3,3,3/)
	integer,dimension(*),parameter :: a = (/3,3,3,2,2,1/)
	integer,dimension(:),allocatable :: c,ia
	integer,dimension(size(a)) :: ic

	print *, "test_unique_int()"
	!print *, "unique_int((/1,2,2,3,3,3/),c,ia,ic)"
	print *, "unique_int((/3,3,3,2,2,1/),c,ia,ic)"
	call unique_int(a,c,ia,ic)
	print *, "c = ", c
	print *, "ia = ", ia
	print *, "ic = ", ic
end subroutine test_unique_int

