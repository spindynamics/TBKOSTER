!
! Copyright (C) 2017
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Mathieu Cesar <mailto:mathieu.cesar@cea.fr>,
! Pascal Thibaudeau <mailto:pascal.thibaudeau@cea.fr>.
!
! This software is a computer program whose purpose is DyNaMol.
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/ or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and inRIA at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeCILL license and that you accept its terms.
!
!  string.f90
!  DyNaMol
module string_mod
  use iso_fortran_env, only: error_unit
  use precision_mod, only: il, ip, rl, rp
  implicit none

  !> Working length of character string variables (arbitrary)
  !> 132 is the maximum length of a line free-form Fortran source line
  integer,parameter :: sl = 132

contains
  !> Eliminate the tab characters and adjust to the left side
  function clean_str(str)
    character(len=*) :: str
    character(len=len(str)) :: clean_str
    integer :: i

    do i=1,len(str)
      if(ichar(str(i:i)) == 9) then
        clean_str(i:i) = ' '
      else
        clean_str(i:i) = str(i:i)
      end if
    end do

    clean_str = adjustl(clean_str)
  end function clean_str

  function cmplx2str(c) result(str)
    complex(rp),intent(in) :: c
    character(len=:),allocatable :: str

    str = '(' // real2str(real(c)) // ', ' // real2str(aimag(c)) // ')'
  end function cmplx2str

 !> Custom flush subroutine to activate with a flag
  subroutine dynamol_flush(unit)
    integer, intent(in)  :: unit
#if defined(Enable_flush)
    call flush(unit=unit)
#endif
  end subroutine dynamol_flush

  !> Convert integer to character string
  function int2str(i) result(str)
    integer(ip),intent(in) :: i
    character(len=il) :: str_temp
    character(len=:),allocatable :: str

    write(str_temp,*) i
    str = trim(adjustl(str_temp))

    !integer,intent(in) :: i
    !character(len=floor(log(real(abs(i),wp))/log(10.0_rp))+2) :: str

    !print *, floor(log(real(abs(i),wp))/log(10.0_rp))+2
    !write(str,*) i
  end function int2str

  !> Test if a character string is an integer
  function is_str_int(str) result(l)
    character(len=*),intent(in) :: str
    logical :: l
    character,parameter :: char_minus = '-'
    character(len=10),parameter :: char_digit = '0123456789'
    integer :: is,id
    integer :: ls

    ls = len(str)
    if(ls == 0) then
      l = .false.
    elseif(ls == 1) then
      l = .false.
      do id=1,10
        l = l .or. (str(1:1) == char_digit(id:id))
        if(l) exit
      end do
    else
      ! Test first character
      l = .false.
      l = l .or. (str(1:1) == char_minus)
      if(.not. l) then
        do id=1,10
          l = l .or. (str(1:1) == char_digit(id:id))
          if(l) exit
        end do
      end if
      if(.not. l) return
      ! Test remaining characters
      do is=2,ls
        l = .false.
        do id=1,10
          l = l .or. (str(is:is) == char_digit(id:id))
          if(l) exit
        end do
        if(.not. l) return
      end do
    end if
  end function is_str_int

  !> Test if a character string is an integer
  function is_str_log(str) result(l)
    character(len=*),intent(in) :: str
    logical :: l

    if(str=='t' .or. str=='.true.' .or. str=='f' .or. str=='.false.') then
      l = .true.
    else
      l = .false.
    end if
  end function is_str_log

  ! Check if a character string is a real
  !function is_str_real(string)

  !> Convert logical to character string
  function log2str(l) result(str)
    logical,intent(in) :: l
    character(len=:),allocatable :: str

    if(l) then
      str = '.true.'
    else
      str = '.false.'
    end if
  end function log2str

  !> Convert lowercase characters to uppercase, copy other characters
  function lower(str1) result(str2)
    character(len=*),intent(in) :: str1
    character(len=len(str1))    :: str2
    character(len=26),parameter :: char_lower = 'abcdefghijklmnopqrstuvwxyz'
    character(len=26),parameter :: char_upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: is,ic

    str2 = str1
    do is=1,len(str1)
      do ic=1,26
        if(str1(is:is) == char_upper(ic:ic)) then
          str2(is:is) = char_lower(ic:ic)
          exit
        end if
      end do
    end do
  end function lower

  ! Convert real to character string
  function real2str(r) result(str)
    real(rp),intent(in) :: r
    character(len=rl) :: str_temp
    character(len=:),allocatable :: str

    write(str_temp,*) r
    str = trim(adjustl(str_temp))
  end function real2str

  ! Convert real to character string in a fixed form
  function fixedreal2str(r) result(str)
    real(rp),intent(in) :: r
    character(len=rl) :: str_temp
    character(len=:),allocatable :: str

    write(str_temp,"(F18.14)") r
    str = trim(adjustl(str_temp))
  end function fixedreal2str

  ! Convert character string to complex
  !function str2cmplx(str) result(r)

  !> Convert character string to integer
  function str2int(str) result(i)
    character(len=*),intent(in) :: str
    integer :: i

    if(.not. is_str_int(trim(adjustl(str)))) then
      write(error_unit,*) 'string%str2int(): character string ', str, &
       ' is not integer'
      error stop
    end if
    read(str,*) i
  end function str2int

  !> Convert lcharacter string to logical
  function str2log(str) result(l)
    character(len=*),intent(in) :: str
    character(len=len(str)) :: str_temp
    logical :: l

    str_temp = lower(str)
    if(.not. is_str_log(str_temp)) then
      write(error_unit,*) 'string%str2log(): character string ', str, &
       ' is not logical'
      error stop
    end if

    if(str_temp=='t' .or. str_temp=='.true.') then
      l = .true.
    else
      l = .false.
    end if
  end function str2log

  ! Convert character string to real
  !function str2real(str) result(r)

  !> Unique values c in a string array a, index ia is such that c = a(ia) and
  !> index ic is such that a = c(ic)
  subroutine unique_str(a, c, ia, ic)
    character(len=*),dimension(:),intent(in) :: a
    !character(len=len(a)),dimension(:),allocatable,intent(out) :: c
    character(len=*),dimension(:),allocatable,intent(out) :: c
    integer,dimension(:),allocatable,intent(out) :: ia
    integer,dimension(size(a)),intent(out) :: ic
    integer :: sa,sc,ia1,ia2,ic1
    integer,dimension(size(a)) :: index_a
    logical,dimension(size(a)) :: mask

    sa = size(a)
    index_a = (/(ia1, ia1=1, sa)/)
    mask = .false.
    do ia1=1,sa-1
      if(.not. mask(ia1)) then
        do ia2=ia1+1,sa
          mask(ia2) = a(ia1) == a(ia2)
          if(mask(ia2)) index_a(ia2) = index_a(ia1)
        end do
      end if
    end do
    sc = sa - count(mask)
    allocate(c(sc),ia(sc))
    c = pack(a, .not. mask)
    ia = pack(index_a, .not. mask)
    do ia1=1,sa
      do ic1=1,sc
        if(index_a(ia1) == ia(ic1)) then
          ic(ia1) = ic1
          exit
        end if
      end do
    end do
  end subroutine unique_str

  !> Convert uppercase characters to lowercase, copy other characters
  function upper(str1) result(str2)
    character(len=*),intent(in) :: str1
    character(len=len(str1))    :: str2
    character(len=26),parameter :: char_lower = 'abcdefghijklmnopqrstuvwxyz'
    character(len=26),parameter :: char_upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer :: is,ic

    str2 = str1
    do is=1,len(str1)
      do ic=1,26
        if(str1(is:is) == char_lower(ic:ic)) then
          str2(is:is) = char_upper(ic:ic)
          exit
        end if
      end do
    end do
  end function upper
end module string_mod
