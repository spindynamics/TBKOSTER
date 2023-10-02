program main
	print *
	call test_constructor()
	print *
	!call test_accessors()
	!print *
	call test_read_namelist()
	print *
	call test_write_namelist()
	print *
end program main

subroutine test_constructor()
	use element_tb_mod
	use string_mod, only: sl
	implicit none

	character(len= 2),dimension(2) :: symbol = (/'Fe','C '/)
	character(len=sl),dimension(2) :: filename = [character(len=sl) :: &
     "/mnt/Data/Dropbox/github.com/TBKOSTER/Codes/TB/TB_PARAM/fe_par_fcc_bcc_sc_gga", &
     "/mnt/Data/Dropbox/github.com/TBKOSTER/Codes/TB/TB_PARAM/c_par_105"]
	type(element_tb) :: obj

	print *, "test_constructor()"
	print *, "obj = element_tb((/'Fe','C '/),2,filename)"
	print *, "obj%print()"
	print *
	obj = element_tb(symbol,2,filename)
	call obj%write_namelist()
end subroutine test_constructor

subroutine test_accessors()
	!use element_tb_mod
	!implicit none

	!type(element_tb) :: obj

	!print *, "test_accessors()"

	!print *, "obj = element_tb('C')"
	!print *, "obj%set_symbol('carbon')"
	!print *, "obj%print()"
	!obj = element_tb('C')
	!call obj%set_symbol('carbon')
	!call obj%print()

	!print *, "obj = element_tb('C')"
	!print *, "obj%get_symbol() ="
	!obj = element_tb('C')
	!print *, obj%get_symbol()
end subroutine test_accessors

subroutine test_read_namelist()
	use element_tb_mod
	implicit none

	type(element_tb) :: obj

	print *, "test_read_namelist()"
	print *, "obj%read_namelist()"
	print *, "obj%print()"
	print *
	call obj%read_namelist()
	call obj%write_namelist()
end subroutine test_read_namelist

subroutine test_write_namelist()
	use element_tb_mod
	use string_mod, only: sl
	implicit none

	character(len= 2),dimension(2) :: symbol = (/'Fe','C '/)
	character(len=sl),dimension(2) :: filename = [character(len=sl) :: &
     "/mnt/Data/Dropbox/github.com/TBKOSTER/Codes/TB/TB_PARAM/fe_par_fcc_bcc_sc_gga", &
     "/mnt/Data/Dropbox/github.com/TBKOSTER/Codes/TB/TB_PARAM/c_par_105"]
	type(element_tb) :: obj

	print *, "test_write_namelist()"
	print *, "obj = element_tb((/'Fe','C '/),2,filename)"
	print *, "obj%write_namelist()"
	obj = element_tb(symbol,2,filename)
	call obj%write_namelist()
end subroutine test_write_namelist
