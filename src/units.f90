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
!  units.f90
!  TBKOSTER
module units_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use constant_mod, only: a_0, m_e, e_ha, e_ry, t_ha, t_ry
  use precision_mod, only: rp
  use string_mod, only: TBKOSTER_flush, lower
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=6),dimension(*),parameter :: property_list = &
   [character(len=6) :: &
   'energy', &
   'length', &
   'time', &
   'mass' &
   ]

  type,public :: units
    !> Energy units ; options:
    !>  'hau': Hartree atomic units eq. Hartree energy (default)
    !>  'rau': Rydberg atomic units eq. Rydberg energy
    !>  'ev' : electronvolts
    character(len=3) :: energy
    !> Length units ; options:
    !>  'hau': Hartree atomic units eq. Bohr radius (default)
    !>  'rau': Rydberg atomic units eq. Bohr radius
    !>  'nm' : nanometers
    !>  'ang': angstroms
    character(len=3) :: length
    !> Time units ; options:
    !>  'hau': Hartree atomic units (default)
    !>  'rau': Rydberg atomic units
    !>  'fs' : femtoseconds
    character(len=3) :: time
    !> Mass units ; options:
    !>  'hau': Hartree atomic units (default)
    !>  'rau': Rydberg atomic units
    !>  'g/mol' : g/mol
    character(len=5) :: mass

  contains
    ! Procedures
    procedure :: convert_energy
    procedure :: convert_length
    procedure :: convert_mass
    procedure :: convert_time
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type units

contains
  subroutine check_direction(direction)
    character(len=*), intent(in) :: direction

    if(direction /= 'from' &
     .and. direction /= 'to') then
      write(error_unit,*) 'units%check_direction(): direction must be one of: &
       &''from'', ''to'''
      error stop
    end if
  end subroutine check_direction

  subroutine check_energy(energy)
    character(len=*), intent(in) :: energy

    if(energy /= 'hau' &
     .and. energy /= 'rau' &
     .and. energy /= 'ev') then
      write(error_unit,*) 'units%check_energy(): energy must be one of: &
       &''hau'', ''rau'', ''ev'''
      error stop
    end if
  end subroutine check_energy

  subroutine check_length(length)
    character(len=*), intent(in) :: length

    if(length /= 'hau' &
     .and. length /= 'rau' &
     .and. length /= 'nm' &
     .and. length /= 'ang') then
      write(error_unit,*) 'units%check_length(): length must be one of: &
       &''hau'', ''rau'', ''nm'', ''ang'''
      error stop
    end if
  end subroutine check_length

  subroutine check_mass(mass)
    character(len=*), intent(in) :: mass

    if(mass /= 'hau' &
     .and. mass /= 'rau' &
     .and. mass /= 'g/mol') then
      write(error_unit,*) 'units%check_mass(): mass must be one of: &
       &''hau'', ''rau'', ''g/mol'''
      error stop
    end if
  end subroutine check_mass

  subroutine check_time(time)
    character(len=*), intent(in) :: time

    if(time /= 'hau' &
     .and. time /= 'rau' &
     .and. time /= 'fs') then
      write(error_unit,*) 'units%check_time(): time must be one of: &
       &''hau'', ''rau'', ''fs'''
      error stop
    end if
  end subroutine check_time

  !> Returns the energy conversion factor from/to the given energy units
  function convert_energy(obj,direction,energy) result(f)
    class(units),intent(in) :: obj
    character(len=*),intent(in) :: direction, energy
    character(len=3) :: energy_1, energy_2
    real(rp) :: f

    call check_direction(direction)
    call check_energy(energy)
    select case(direction)
    case('from')
      energy_1 = energy
      energy_2 = obj%energy
    case('to')
      energy_1 = obj%energy
      energy_2 = energy
    end select

    select case(trim(energy_1))
    case('hau')
      select case(trim(energy_2))
      case('hau')
        f = 1.0_rp
      case('rau')
        f = 2.0_rp
      case('ev')
        f = e_ha
      end select
    case('rau')
      select case(trim(energy_2))
      case('hau')
        f = 0.5_rp
      case('rau')
        f = 1.0_rp
      case('ev')
        f = e_ry
      end select
    case('ev')
      select case(trim(energy_2))
      case('hau')
        f = 1.0_rp/e_ha
      case('rau')
        f = 1.0_rp/e_ry
      case('ev')
        f = 1.0_rp
      end select
    end select
  end function convert_energy

  !> Returns the length conversion factor from/to the given length units
  function convert_length(obj,direction,length) result(f)
    class(units),intent(in) :: obj
    character(len=*),intent(in) :: direction, length
    character(len=3) :: length_1, length_2
    real(rp) :: f

    call check_direction(direction)
    call check_length(length)
    select case(direction)
    case('from')
      length_1 = trim(length)
      length_2 = trim(obj%length)
    case('to')
      length_1 = trim(obj%length)
      length_2 = trim(length)
    end select

    select case(length_1)
    case('hau')
      select case(length_2)
      case('hau')
        f = 1.0_rp
      case('rau')
        f = 1.0_rp
      case('nm')
        f = a_0/10.0_rp
      case('ang')
        f = a_0
      end select
    case('rau')
      select case(length_2)
      case('hau')
        f = 1.0_rp
      case('rau')
        f = 1.0_rp
      case('nm')
        f = a_0/10.0_rp
      case('ang')
        f = a_0
      end select
    case('nm')
      select case(length_2)
      case('hau')
        f = 10.0_rp/a_0
      case('rau')
        f = 10.0_rp/a_0
      case('nm')
        f = 1.0_rp
      case('ang')
        f = 10.0_rp
      end select
    case('ang')
      select case(length_2)
      case('hau')
        f = 1.0_rp/a_0
      case('rau')
        f = 1.0_rp/a_0
      case('nm')
        f = 0.1_rp
      case('ang')
        f = 1.0_rp
      end select
    end select
  end function convert_length

  !> Returns the mass conversion factor from/to the given mass units
  function convert_mass(obj,direction,mass) result(f)
    class(units),intent(in) :: obj
    character(len=*),intent(in) :: direction, mass
    character(len=5) :: mass_1, mass_2
    real(rp) :: f

    call check_direction(direction)
    call check_mass(mass)
    select case(direction)
    case('from')
      mass_1 = trim(mass)
      mass_2 = trim(obj%mass)
    case('to')
      mass_1 = trim(obj%mass)
      mass_2 = trim(mass)
    end select

    select case(mass_1)
    case('hau')
      select case(mass_2)
      case('hau')
        f = 1.0_rp
      case('rau')
        f = 1.0_rp
      case('g/mol')
        f = 1.0_rp
      end select
    case('rau')
      select case(mass_2)
      case('hau')
        f = 1.0_rp
      case('rau')
        f = 1.0_rp
      case('g/mol')
        f = 1.0_rp
      end select
    case('g/mol')
      select case(mass_2)
      case('hau')
        f = 1.0_rp
      case('rau')
        f = 1.0_rp
      case('g/mol')
        f = 1.0_rp
      end select
    end select
  end function convert_mass

  !> Returns the time conversion factor from/to the given time units
  function convert_time(obj,direction,time) result(f)
    class(units),intent(in) :: obj
    character(len=*),intent(in) :: direction, time
    character(len=3) :: time_1, time_2
    real(rp) :: f

    call check_direction(direction)
    call check_time(time)
    select case(direction)
    case('from')
      time_1 = time
      time_2 = obj%time
    case('to')
      time_1 = obj%time
      time_2 = time
    end select

    select case(trim(time_1))
    case('hau')
      select case(trim(time_2))
      case('hau')
        f = 1.0_rp
      case('rau')
        f = 0.5_rp
      case('fs')
        f = t_ha
      end select
    case('rau')
      select case(trim(time_2))
      case('hau')
        f = 2.0_rp
      case('rau')
        f = 1.0_rp
      case('fs')
        f = t_ry
      end select
    case('fs')
      select case(trim(time_2))
      case('hau')
        f = 1.0_rp/t_ha
      case('rau')
        f = 1.0_rp/t_ry
      case('fs')
        f = 1.0_rp
      end select
    end select
  end function convert_time

  subroutine initialize(energy,length,time,mass)
    character(len=3),intent(out) :: energy, length, time, mass

    energy = 'hau'
    length = 'hau'
    time = 'hau'
    mass = 'hau'
  end subroutine initialize

  !> Read object in text format from file (default: 'in_units.txt')
  subroutine read_txt(obj,file)
    class(units),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    character(len=5) :: energy,length,time,mass
    ! Namelist
    namelist /units/ energy,length,time,mass

    write(output_unit,*) "DEBUG == Entering units & read_txt"
    call TBKOSTER_flush(output_unit)
    
    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_units.txt'
    end if

    inquire(unit=10, opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'units%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'units%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    call initialize(energy,length,time,mass)
    read(10,nml=units)
    energy = lower(energy)
    call check_energy(trim(energy))
    length = lower(length)
    call check_length(trim(length))
    time = lower(time)
    call check_time(trim(time))
    mass = lower(mass)
    call check_mass(trim(mass))

    obj%energy = energy
    obj%length = length
    obj%time = time
    obj%mass = mass

    close(unit=10)
    !deallocate(file_rt)
    write(output_unit,*) "DEBUG == Exiting units & read_txt"
    call TBKOSTER_flush(output_unit)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_units.txt')
  subroutine write_txt(obj,file,unit)
    class(units),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=len(obj%energy)) :: energy
    character(len=len(obj%length)) :: length
    character(len=len(obj%time))   :: time
    character(len=len(obj%mass))   :: mass
    ! Namelist
    namelist /units/ energy, length, time, mass

    if (present(file)) then
      file_rt = file
    else
      file_rt = 'out_units.txt'
    end if

    if (present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    energy = obj%energy
    length = obj%length
    time   = obj%time
    mass   = obj%mass

    write(unit_rt,nml=units)
    
    if(.not. present(unit)) then
      call TBKOSTER_flush(unit_rt)
      close(unit_rt)
    else
      call TBKOSTER_flush(unit)
      close(unit)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_units.txt'), if tag (default: .true.) the namelist opening and closing
  !> tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(units),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Local variables
    integer :: ip

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_units.txt'
    end if
    if(present(property)) then
      property_rt = property
    else
      property_rt = property_list
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
      write(unit_rt,'(a)') '&units'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('energy')
        write(unit_rt,'(a)') ' energy = ''' // trim(obj%energy) // ''''
      case('length')
        write(unit_rt,'(a)') ' length = ''' // trim(obj%length) // ''''
      case('time')
        write(unit_rt,'(a)') ' time = ''' // trim(obj%time) // ''''
      case('mass')
        write(unit_rt,'(a)') ' mass = ''' // trim(obj%mass) // ''''
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
end module units_mod
