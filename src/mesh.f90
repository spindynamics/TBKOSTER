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
!  mesh.f90
!  TBKOSTER
module mesh_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use lattice_mod
  use precision_mod, only: rp
  use string_mod, only: TBKOSTER_flush, int2str, lower, real2str, sl
  use units_mod
  implicit none
  private

  type,public :: mesh
    !> Units
    class(units),pointer :: u
    !> Lattice
    class(lattice),pointer :: l

    !> Mesh type ; options:
    !>  'mp'  : Monkhorst-Pack scheme
    !>  'list': explicit list of points
    !>  'path': path of high-symmetry points
    character(len=4) :: type
    !> Divisions between reciprocal space unit vectors (type='mp')
    integer,dimension(3) :: gx
    !> Shift (type='mp')
    integer,dimension(3) :: dx
    !> Divisions between high-symmetry points (type='path')
    integer :: gxs
    !> Number of high-symmetry points (type='path')
    integer :: nxs
    !> High-symmetry point labels (type='path')
    character(len=2),dimension(:),allocatable :: xs_label
    !> High-symmetry points coordinates (type='path'); options:
    !>  'cartesian'
    !>  'direct' (default)
    character(len=9) :: xs_coord
    !> High-symmetry points (type='path')
    real(rp),dimension(:,:),allocatable :: xs
    !> Number of points
    integer :: nx
    !> Points coordinates; options:
    !>  'cartesian'
    !>  'direct' (default)
    character(len=9) :: x_coord
    !> Points
    real(rp),dimension(:,:),allocatable :: x
    !> Weights
    real(rp),dimension(:),allocatable :: w

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: cart2dir
    procedure :: dir2cart
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type mesh

  interface mesh
    procedure :: constructor
  end interface mesh

contains
  function constructor(l) result(obj)
    class(lattice),target,intent(in) :: l
    type(mesh) :: obj

    obj%u => l%u
    obj%l => l
  end function constructor

  subroutine destructor(obj)
    type(mesh) :: obj

    if(allocated(obj%xs_label)) deallocate(obj%xs_label)
    if(allocated(obj%xs))       deallocate(obj%xs)
    if(allocated(obj%x))        deallocate(obj%x)
    if(allocated(obj%w))        deallocate(obj%w)
  end subroutine destructor

  function build_monkhorst_pack(gx,dx) result(x)
    integer,dimension(3),intent(in) :: gx,dx
    real(rp),dimension(:,:),allocatable :: x
    integer :: ix,igx1,igx2,igx3
    real(rp) :: x1,x2,x3

    allocate(x(product(gx),3))

    ix = 1
    do igx1=1,gx(1)
      x1 = (2*igx1-gx(1)-1.0_rp+dx(1))/(2*gx(1))
      do igx2=1,gx(2)
        x2 = (2*igx2-gx(2)-1.0_rp+dx(2))/(2*gx(2))
        do igx3=1,gx(3)
          x3 = (2*igx3-gx(3)-1.0_rp+dx(3))/(2*gx(3))
          x(ix,:) = (/x1,x2,x3/)
          ix = ix+1
        end do
      end do
    end do
  end function build_monkhorst_pack

  function build_path(gxs,xs) result(x)
    integer,intent(in) :: gxs
    real(rp),dimension(:,:),intent(in) :: xs
    real(rp),dimension(:,:),allocatable :: x
    integer :: nxs,ix,ixs,igxs

    allocate(x((size(xs,1)-1)*gxs+1,3))

    nxs = size(xs,1)
    ix = 1
    do ixs=1,nxs-1
      do igxs=0,gxs-1
        x(ix,:) = xs(ixs,:) + (xs(ixs+1,:)-xs(ixs,:))*igxs/gxs
        ix = ix+1
      end do
    end do
    x(ix,:) = xs(nxs,:)
  end function build_path

  subroutine cart2dir(obj)
    class(mesh),intent(inout) :: obj

    if(obj%nxs > 0 .and. obj%xs_coord == 'cartesian') then
      obj%xs_coord = 'direct'
      obj%xs = obj%l%cart2dir(obj%xs)
    end if
    if(obj%x_coord == 'cartesian') then
      obj%x_coord = 'direct'
      obj%x = obj%l%cart2dir(obj%x)
    end if
  end subroutine cart2dir

  !> Check the validity of mesh type
  subroutine check_type(type)
    character(len=*) :: type

    if(type /= 'list' &
     .and. type /= 'mp' &
     .and. type /= 'path') then
      write(error_unit,*) 'mesh%check_type(): mesh%type must be one of: &
       ''list'', ''mp'', ''path'''
      error stop
    end if
  end subroutine check_type

  !> Check the validity of the position coordinates
  subroutine check_x_coord(x_coord)
    character(len=*),intent(in) :: x_coord

    if(x_coord /= 'cartesian' .and. x_coord /= 'direct') then
      write(error_unit,*) 'atom%check_x_coord(): x_coord must be one of: &
       &''cartesian'', ''direct'''
      error stop
    end if
  end subroutine check_x_coord

  subroutine dir2cart(obj)
    class(mesh),intent(inout) :: obj

    if(obj%nxs > 0 .and. obj%xs_coord == 'direct') then
      obj%xs_coord = 'cartesian'
      obj%xs = obj%l%dir2cart(obj%xs)
    end if
    if(obj%x_coord == 'direct') then
      obj%x_coord = 'cartesian'
      obj%x = obj%l%dir2cart(obj%x)
    end if
  end subroutine dir2cart

  subroutine initialize_monkhorst_pack(gx,dx)
    integer,dimension(3),intent(out) :: gx,dx

    gx = 0
    dx = 0
  end subroutine initialize_monkhorst_pack

  subroutine initialize_path(gxs,nxs)
    integer,intent(out) :: gxs,nxs

    gxs = 0
    nxs = 0
  end subroutine initialize_path

  !> Read object in text format from file (default: 'in_mesh.txt')
  subroutine read_txt(obj,file)
    class(mesh),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    character(len=4) :: type
    integer,dimension(3) :: gx, dx
    integer :: gxs, nxs
    character(len=9) :: xs_coord
    character(len=2),dimension(:),allocatable :: xs_label
    real(rp),dimension(:,:),allocatable :: xs
    integer :: nx
    character(len=9) :: x_coord
    real(rp),dimension(:,:),allocatable :: x
    real(rp),dimension(:),allocatable :: w
    ! Namelist
    namelist /mesh/ type, gx, dx, gxs, nxs, xs_coord, xs_label, xs, nx, &
     x_coord, x, w

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_mesh.txt'
    end if
    
    inquire(unit=10, opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'mesh%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'mesh%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    call initialize_monkhorst_pack(gx,dx)
    call initialize_path(gxs,nxs)
    allocate(xs_label(0))
    read(10,nml=mesh,iostat=iostatus)
    type = lower(type)
    call check_type(type)
    select case(type)
    case('list')
      x_coord = 'direct'
      allocate(x(nx,3))
      allocate(w(nx))
      w = 1.0_rp/nx
      rewind(10)
      read(10,nml=mesh)
      x_coord = lower(x_coord)
      call check_x_coord(trim(x_coord))
      if(trim(x_coord) == 'cartesian') then
        x_coord = 'direct'
        x = obj%l%cart2dir(x)
      end if
    case('mp')
      nx = product(gx)
      allocate(x(nx,3))
      allocate(w(nx))
      x_coord = 'direct'
      x = build_monkhorst_pack(gx,dx)
      w = 1.0_rp/nx
    case('path')
      deallocate(xs_label)
      allocate(xs_label(nxs))
      xs_coord = 'direct'
      allocate(xs(nxs,3))
      rewind(10)
      read(10,nml=mesh)
      xs_coord = lower(xs_coord)
      call check_x_coord(trim(xs_coord))
      if(trim(xs_coord) == 'cartesian') then
        xs_coord = 'direct'
        xs = obj%l%cart2dir(xs)
      end if
      nx = (nxs-1)*gxs+1
      allocate(x(nx,3))
      allocate(w(nx))
      x_coord = 'direct'
      x = build_path(gxs,xs)
      w = 1.0_rp/nx
    end select

    obj%type = type
    obj%gx = gx
    obj%dx = dx
    obj%gxs = gxs
    obj%nxs = nxs
    if(nxs > 0) then
      call move_alloc(xs_label,obj%xs_label)
    end if
    obj%xs_coord = xs_coord
    if(nxs > 0) then
      if(xs_coord == 'cartesian') then
        xs = xs / obj%u%convert_length('to','hau')
      end if
      call move_alloc(xs,obj%xs)
    end if
    obj%nx = nx
    obj%x_coord = x_coord
    if(x_coord == 'cartesian') then
      x = x / obj%u%convert_length('to','hau')
    end if
    call move_alloc(x,obj%x)
    call move_alloc(w,obj%w)

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_mesh.txt')
  subroutine write_txt(obj,file,unit)
    class(mesh),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=len(obj%type)) :: type
    integer,dimension(3) :: gx, dx
    integer :: gxs, nxs
    character(len=2),dimension(obj%nxs) :: xs_label
    character(len=9) :: xs_coord
    real(rp),dimension(obj%nxs,3) :: xs
    integer :: nx
    character(len=9) :: x_coord
    real(rp),dimension(obj%nx,3) :: x
    real(rp),dimension(obj%nx) :: w
    ! Namelist
    namelist /mesh/ type, gx, dx, gxs, nxs, xs_coord, xs_label, xs, nx, &
     x_coord, x, w

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_mesh.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    type = obj%type
    gx = obj%gx
    dx = obj%dx
    gxs = obj%gxs
    nxs = obj%nxs
    xs_label = obj%xs_label
    xs_coord = obj%xs_coord
    if(obj%nxs > 0 .and. obj%xs_coord == 'cartesian') then
      xs = obj%xs / obj%u%convert_length('from','hau')
    else
      xs = obj%xs
    end if
    nx = obj%nx
    x_coord = obj%x_coord
    if(obj%x_coord == 'cartesian') then
      x = obj%x / obj%u%convert_length('from','hau')
    else
      x = obj%x
    end if
    w = obj%w

    write(unit_rt,nml=mesh)
    call TBKOSTER_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: see source code) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_mesh.txt'), if tag (default: .true.) the namelist opening and closing
  !> tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(mesh),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    real(rp),dimension(obj%nxs,3) :: xs
    real(rp),dimension(obj%nx ,3) :: x
    ! Local variables
    integer :: ip, ix, ixs

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_mesh.txt'
    end if
    if(present(property)) then
      property_rt = property
    else
      select case(obj%type)
      case('list')
        property_rt = [character(len=sl) :: 'type','nx','x_coord','x','w']
      case('mp')
        property_rt = [character(len=sl) :: 'type','gx','dx','nx']
      case('path')
        property_rt = [character(len=sl) :: 'type','gxs','xs_label','xs_coord',&
         'xs','nx']
      end select
    end if
    if(present(tag)) then
      tag_rt = tag
    else
      tag_rt = .true.
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if
    if(tag_rt) then
      write(unit_rt,'(a)') '&mesh'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('type')
        write(unit_rt,'(a)') ' type = ''' // trim(obj%type) // ''''
      case('gx')
        write(unit_rt,'(a)') ' gx = ' &
         // int2str(obj%gx(1)) // ', ' &
         // int2str(obj%gx(2)) // ', ' &
         // int2str(obj%gx(3))
      case('dx')
        write(unit_rt,'(a)') ' dx = ' &
         // int2str(obj%dx(1)) // ', ' &
         // int2str(obj%dx(2)) // ', ' &
         // int2str(obj%dx(3))
      case('gxs')
        write(unit_rt,'(a)') ' gxs = ' // int2str(obj%gxs)
      case('nxs')
        write(unit_rt,'(a)') ' nxs = ' // int2str(obj%nxs)
      case('xs_label')
        do ixs=1,obj%nxs
          write(unit_rt,'(a)') ' xs_label(' // int2str(ixs) // ') = ''' &
           // trim(obj%xs_label(ixs)) // ''''
        end do
      case('xs_coord')
        write(unit_rt,'(a)') ' xs_coord = ''' // trim(obj%xs_coord) // ''''
      case('xs')
        if(obj%nxs > 0 .and. obj%xs_coord == 'cartesian') then
          xs = obj%xs / obj%u%convert_length('from','hau')
        else
          xs = obj%xs
        end if
        do ixs=1,obj%nxs
          write(unit_rt,'(a)') ' xs(' // int2str(ixs) // ',:) = ' &
           // real2str(xs(ixs,1)) // ', ' &
           // real2str(xs(ixs,2)) // ', ' &
           // real2str(xs(ixs,3))
        end do
      case('nx')
        write(unit_rt,'(a)') ' nx = ' // int2str(obj%nx)
      case('x_coord')
        write(unit_rt,'(a)') ' x_coord = ''' // trim(obj%x_coord) // ''''
      case('x')
        if(obj%x_coord == 'cartesian') then
          x = obj%x / obj%u%convert_length('from','hau')
        else
          x = obj%x
        end if
        do ix=1,obj%nx
          write(unit_rt,'(a)') ' x(' // int2str(ix) // ',:) = ' &
           // real2str(x(ix,1)) // ', ' &
           // real2str(x(ix,2)) // ', ' &
           // real2str(x(ix,3))
        end do
      case('w')
        do ix=1,obj%nx
          write(unit_rt,'(a)') ' w(' // int2str(ix) // ') = ' &
           // real2str(obj%w(ix))
        end do
      end select
    end do

    if(tag_rt) then
      write(unit_rt,'(a)') ' /'
    end if
    call TBKOSTER_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt,property_rt)
  end subroutine write_txt_formatted
end module mesh_mod
