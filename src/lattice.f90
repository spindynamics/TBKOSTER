!
! Copyright (C) 2017
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Mathieu Cesar <mailto:mathieu.cesar@cea.fr>,
! Pascal Thibaudeau <mailto:pascal.thibaudeau@cea.fr>.
!
! This software is a computer program whose purpose is TBKOSTER.
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
!  lattice.f90
!  TBKOSTER
module lattice_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
!#if defined(BLAS95_FOUND)
!  use blas95, only: gemm
!#endif
  use math_mod, only: two_pi, cross_product, determinant, inverse_3x3
  use precision_mod, only: rp
  use string_mod, only: sl, int2str, lower, real2str, TBKOSTER_flush
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'v' &
   ]

  type,public :: lattice
    !> Units
    class(units),pointer :: u

    !> Vectors
    real(rp),dimension(3,3) :: v
    !> Vectors inverse
    real(rp),dimension(3,3) :: vi

  contains
    ! Procedures
    procedure :: construct_reciprocal
    procedure :: cart2dir
    procedure :: dir2cart
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type lattice

  ! Constructor
  interface lattice
    procedure :: constructor
  end interface lattice

contains
  function constructor(u) result(obj)
    class(units),target,intent(in) :: u
    type(lattice) :: obj

    obj%u => u
  end function constructor

  !> Parse vectors from cartesian to direct coordinates
  function cart2dir(obj,v_car) result(v_dir)
    class(lattice),intent(in) :: obj
    real(rp),dimension(:,:),intent(in) :: v_car
    real(rp),dimension(size(v_car,1),size(v_car,2)) :: v_dir
! #if !defined(BLAS95_FOUND)
!     integer :: m, n, k
!     integer :: lda, ldb, ldc
! #endif
!
! #if defined(BLAS95_FOUND)
!     call gemm(v_car,obj%vi,v_dir)
! #else
!     m = size(v_car,1)
!     n = size(obj%vi,2)
!     k = size(obj%vi,1)
!     lda = max(1,m)
!     ldb = max(1,n)
!     ldc = max(1,m)
!     call dgemm('N','N',m,n,k,1.0_rp,v_car,lda,obj%vi,ldb,0.0_rp,v_dir,ldc)
! #endif

    v_dir = matmul(v_car,obj%vi)
  end function cart2dir

  !> Construct reciprocal lattice
  function construct_reciprocal(obj1) result(obj2)
    class(lattice),intent(in) :: obj1
    type(lattice)  :: obj2
    real(rp),dimension(3,3) :: v

    real(rp), dimension(3) :: v1,v2,v3

    v1=obj1%v(1,:)
    v2=obj1%v(2,:)
    v3=obj1%v(3,:)

    v(1,:) = cross_product(v2,v3)
    v(2,:) = cross_product(v3,v1)
    v(3,:) = cross_product(v1,v2)

    v = v/abs(determinant(obj1%v))
    !v = two_pi*v/abs(determinant(obj%v))

    obj2%v = v
    obj2%vi = inverse_3x3(v)
    obj2%u => obj1%u
  end function construct_reciprocal

  !> Parse vectors from direct to cartesian coordinates
  function dir2cart(obj,v_dir) result(v_car)
    class(lattice),intent(in) :: obj
    real(rp),dimension(:,:),intent(in) :: v_dir
    real(rp),dimension(size(v_dir,1),size(v_dir,2)) :: v_car
    ! #if !defined(BLAS95_FOUND)
    !     integer :: m, n, k
    !     integer :: lda, ldb, ldc
    ! #endif
    !
    ! #if defined(BLAS95_FOUND)
    !     call gemm(v_dir,obj%v,v_car)
    ! #else
    !     m = size(v_dir,1)
    !     n = size(obj%v,2)
    !     k = size(obj%v,1)
    !     lda = max(1,m)
    !     ldb = max(1,n)
    !     ldc = max(1,m)
    !     call dgemm('N','N',m,n,k,1.0_rp,v_dir,lda,obj%v,ldb,0.0_rp,v_car,ldc)
    ! #endif

    v_car = matmul(v_dir,obj%v)
  end function dir2cart

  !> Read object in text format from file (default: 'in_lattice.txt')
  subroutine read_txt(obj,file)
    class(lattice),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    real(rp) :: v_factor
    real(rp),dimension(3,3) :: v
    ! Namelist
    namelist /lattice/ v_factor, v

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_lattice.txt'
    end if

    inquire(unit=10, opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'lattice%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'lattice%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    v_factor = 1.0_rp
    read(10,nml=lattice)
    v = v_factor*v

    obj%v = v * obj%u%convert_length('to','hau')
    obj%vi = inverse_3x3(obj%v)

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_lattice.txt')
  subroutine write_txt(obj,file,unit)
    class(lattice),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    real(rp),dimension(3,3) :: v
    ! Namelist
    namelist /lattice/ v

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_lattice.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    v = obj%v * obj%u%convert_length('from','hau')

    write(unit_rt,nml=lattice)
    call TBKOSTER_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: 'all') in text format to unit (default: 10),
  !> if it's a file its name is set to file (default: 'out_lattice.txt'), if
  !> status (default: .true.) the file is created otherwise it is continued
  subroutine write_txt_formatted(obj,file,property,status,unit)
    class(lattice),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: status
    logical                     :: status_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    real(rp),dimension(3,3) :: v
    real(rp),dimension(3,3) :: vi
    ! Local variables
    integer :: iv, ip

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_lattice.txt'
    end if
    if(present(property)) then
      property_rt = property
    else
      property_rt = property_list
    end if
    if(present(status)) then
      status_rt = status
    else
      status_rt = .true.
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(status_rt) then
      if(.not. present(unit)) then
        open(unit=unit_rt,file=file_rt,action='write')
      end if
      write(unit_rt,'(a)') '&lattice'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('v')
        v = obj%v * obj%u%convert_length('from','hau')
        do iv=1,3
          write(unit_rt,'(a)') ' v(' // int2str(iv) // ',:) = ' &
           // real2str(v(iv,1)) // ', ' &
           // real2str(v(iv,2)) // ', ' &
           // real2str(v(iv,3))
        end do
      case('vi')
        vi = obj%vi / obj%u%convert_length('from','hau')
        do iv=1,3
          write(unit_rt,'(a)') ' vi(' // int2str(iv) // ',:) = ' &
           // real2str(vi(iv,1)) // ', ' &
           // real2str(vi(iv,2)) // ', ' &
           // real2str(vi(iv,3))
        end do
      end select
    end do

    if(status_rt) then
      write(unit_rt,'(a)') ' /'
      call TBKOSTER_flush(unit_rt)
      if(.not. present(unit)) then
        close(unit_rt)
      end if
    end if
    !deallocate(file_rt,property_rt)
  end subroutine write_txt_formatted
end module lattice_mod
