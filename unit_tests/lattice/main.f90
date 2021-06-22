program main
	print *
	call test_constructor()
	print *
	call test_accessors()
	print *
	call test_cart2dir()
	print *
	call test_construct_reciprocal()
	print *
	call test_dir2cart()
	print *
	call test_read_namelist()
	print *
	call test_write_namelist()
	print *
end program main

subroutine test_constructor()
	use lattice_mod
	use math_mod
	implicit none

	type(lattice) :: obj

	print *, "test_constructor()"
	print *, "obj = lattice(zeros2(3,3))"
	print *, "obj%print()"
	obj = lattice(zeros2(3,3))
	call obj%write_namelist()
end subroutine test_constructor

subroutine test_accessors()
	use lattice_mod
	use math_mod
	use precision_mod
	implicit none

	type(lattice) :: obj
	real(rp), dimension(3,3) :: v

	print *, "test_accessors()"

	print *, "obj%set_vector(zeros2(3,3))"
	print *, "obj%print()"
	!call obj%set_vector(zeros2(3,3))
	call obj%write_namelist()

	print *, "obj%get_vector() ="
	v = obj%get_v()
	print *, v(1,:)
	print *, v(2,:)
	print *, v(3,:)
end subroutine test_accessors

subroutine test_cart2dir()
	use lattice_mod
	use math_mod
	use precision_mod
	implicit none

	type(lattice) :: obj
	real(rp), dimension(3,3) :: v

	print *, "test_cart2dir()"
	print *, "obj = lattice(diag1((/1.0,2.0,3.0/)))"
	print *, "obj%cart2dir(eye(3)) ="
	obj = lattice(diag1((/1.0_rp,2.0_rp,3.0_rp/)))
	v = obj%cart2dir(eye(3))
	print *, v(1,:)
	print *, v(2,:)
	print *, v(3,:)
end subroutine test_cart2dir

subroutine test_construct_reciprocal()
	use lattice_mod
	use math_mod
	use precision_mod
	implicit none

	type(lattice) :: obj1,obj2
	real(rp), dimension(3,3) :: v

	print *, "test_construct_reciprocal()"
	print *, "obj1 = lattice(diag1((/1.0,2.0,3.0/)))"
	print *, "obj2 = obj1%construct_reciprocal()"
	print *, "obj2%print()"
	obj1 = lattice(diag1((/1.0_rp,2.0_rp,3.0_rp/)))
	obj2 = obj1%construct_reciprocal()
	call obj2%write_namelist()
end subroutine test_construct_reciprocal

subroutine test_dir2cart()
	use lattice_mod
	use math_mod
	use precision_mod
	implicit none

	type(lattice) :: obj
	real(rp), dimension(3,3) :: v

	print *, "test_dir2cart()"
	print *, "obj = lattice(diag1((/1.0,2.0,3.0/)))"
	print *, "obj%dir2cart(eye(3)) ="
	obj = lattice(diag1((/1.0_rp,2.0_rp,3.0_rp/)))
	v = obj%dir2cart(eye(3))
	print *, v(1,:)
	print *, v(2,:)
	print *, v(3,:)
end subroutine test_dir2cart

subroutine test_read_namelist()
	use lattice_mod
	implicit none

	type(lattice) :: obj

	print *, "test_read_namelist()"
	print *, "obj%read_namelist()"
	print *, "obj%print()"
	call obj%read_namelist()
	call obj%write_namelist()
end subroutine test_read_namelist

subroutine test_write_namelist()
	use lattice_mod
	use math_mod
	use precision_mod
	implicit none

	type(lattice) :: obj

	print *, "test_write_namelist()"
	print *, "obj = lattice(diag1((/1.0,2.0,3.0/)))"
	print *, "obj%write_namelist()"
	obj = lattice(diag1((/1.0_rp,2.0_rp,3.0_rp/)))
	call obj%write_namelist()
end subroutine test_write_namelist
