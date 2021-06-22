program main
	print *
	!call test_constructor()
	!print *
	!call test_accessors()
	!print *
	!call test_build_ia2ie()
	!print *
	call test_read_namelist()
	print *
end program main

subroutine test_constructor()
	!use atom_mod
	!implicit none

	!type(atom) :: obj

	!print *, "test_constructor()"
end subroutine test_constructor

subroutine test_accessors()
	!use atom_mod
	!implicit none

	!type(atom) :: obj

	!print *, "test_accessors()"
end subroutine test_accessors

subroutine test_build_ia2ie()
	use atom_mod
	implicit none

	type(atom) :: obj

	print *, "test_build_ia2ie()"
	print *, "print *, build_ia2ie((/'Fe_atom1','Fe_atom2','Fe_atom3','Fe_atom4'/))"
	!print *, build_ia2ie((/'Fe_atom1','Fe_atom2','Fe_atom3','Fe_atom4'/))
	print *
	print *, "print *, build_ia2ie((/'Fe_atom1','C _atom1','Fe_atom2','C _atom2'/))"
	!print *, build_ia2ie((/'Fe_atom1','C _atom1','Fe_atom2','C _atom2'/))
end subroutine test_build_ia2ie

subroutine test_read_namelist()
	use, intrinsic :: iso_fortran_env, only: output_unit
	use atom_mod
	implicit none

	type(atom) :: obj

	print *, "test_read_namelist()"
	print *, "obj%read_namelist()"
	print *, "obj%write_namelist_formatted()"
	call obj%read_namelist()
	call obj%write_namelist_formatted(unit=output_unit)
end subroutine test_read_namelist
