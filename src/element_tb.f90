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
!  element_tb.f90
!  TBKOSTER
module element_tb_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use element_mod
  use precision_mod, only: rp
  use string_mod, only: sl, int2str, lower, TBKOSTER_flush
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'filename' &
   ]

  !> Maximum number of NRL TB parameters
  integer,parameter :: np_max = 97
  
  type,public,extends(element) :: element_tb
    !> @defgroup TB_parameters TB parameter variables
    !> See \cite Barreteau2016 page 7-8, \cite Mehl1996 page 1-4.\n
    !> A TB calculation is completely defined by a finite set of parameters.
    !> These parameters are distributed between the on-site energies, the
    !> hopping integrals and the overlap integrals. The on-site energies are
    !> written as:\n
    !> \f$ \epsilon_{i\lambda} = a_{\lambda}
    !> + b_{\lambda}\rho_i^{1/3} + c_{\lambda}\rho_i^{2/3}
    !> + d_{\lambda}\rho_i^{4/3} + e_{\lambda}\rho_i^2 \f$\n
    !> with the pseudo-density:\n
    !> \f$ \rho_i = \sum_{j \neq i} \mathrm{exp}[-\Lambda^2 R_{ij}]
    !> F_c(R_{ij}) \f$\n
    !> where the parameters are \f$ \Lambda \f$ and \f$ a_{\lambda},b_{\lambda},
    !> c_{\lambda},d_{\lambda},e_{\lambda} \f$ with \f$ \lambda = s,p,d \f$.
    !> The hopping and overlap integrals are written as:\n
    !> \f$ \beta_{\gamma}(R) = (p_{\gamma} + f_{\gamma}R + g_{\gamma}R^2)
    !> \mathrm{exp}[-h_{\gamma}^2 R] F_c(R) \f$\n
    !> where the parameters are \f$ p_{\gamma},f_{\gamma},g_{\gamma},h_{\gamma}
    !> \f$ with \f$ \gamma = ss\sigma, sp\sigma, sd\sigma, pp\sigma, pp\pi,
    !> pd\sigma, pd\pi, dd\sigma, dd\pi, dd\delta \f$.
    !> @{

    !> TB parameter filenames
    character(len=sl),dimension(:),allocatable :: filename
    !> TB model types ; options: 'nrl', 'mod', 'wan'
    character(len=3) :: tb_type
    !> TB parameter types ; options: 'old', 'new', 'cyr', 'pow'
    character(len=3),dimension(:),allocatable :: nrl_type
    !> TB parameters of NRL type
    !> TB parameter \f$ \Lambda \f$
    !real(rp),dimension(:),allocatable :: lambda = abs(p(:,1)) (absolute value) ?
    !> TB parameters:\n
    !> (1): \f$ \Lambda \f$\n
    !> (2:16): \f$ a_{\lambda},b_{\lambda},c_{\lambda},d_{\lambda},e_{\lambda}
    !> \f$ with \f$ \lambda = s,p,d \f$\n
    !> (17): unused\n
    !> (18:97): \f$ p_{\gamma},f_{\gamma},g_{\gamma},h_{\gamma} \f$ with
    !> \f$ \gamma = ss\sigma, sp\sigma, sd\sigma, pp\sigma, pp\pi, pd\sigma,
    !> pd\pi, dd\sigma, dd\pi, dd\delta \f$ for hopping integrals, then for
    !> overlap integrals
    real(rp),dimension(:,:),allocatable :: p
    !> @}

    !> @defgroup Cutoff Cutoff function variables
    !> See \cite Barreteau2016 page 7-8.\n
    !> Both the on-site energies:\n
    !> \f$ \epsilon_{i\lambda} = a_{\lambda}
    !> + b_{\lambda}\rho_i^{1/3} + c_{\lambda}\rho_i^{2/3}
    !> + d_{\lambda}\rho_i^{4/3} + e_{\lambda}\rho_i^2 \f$\n
    !> with the pseudo-density:\n
    !> \f$ \rho_i = \sum_{j \neq i} \mathrm{exp}[-\Lambda^2 R_{ij}] F_c(R_{ij})
    !> \f$\n
    !> and the overlap/hopping integrals:\n
    !> \f$ \beta_{\gamma}(R) = (p_{\gamma} + f_{\gamma}R + g_{\gamma}R^2)
    !> \mathrm{exp}[-h_{\gamma}^2 R] F_c(R) \f$\n
    !> involve a cutoff function combining an exponential decay and a step
    !> function:\n
    !> \f$ F_c(R) = \frac{1-\theta(R-r_c)}{1+\mathrm{exp}[(R-r_0)/r_l]} \f$\n
    !> where \f$ r_0 \f$ is the exponential offset, \f$ r_l \f$ is the
    !> exponential lifetime, \f$ \theta(R) \f$ is the Heaviside step function
    !> and \f$ r_c = r_0 + 5r_l \f$ is the step function offset
    !> @{

    !> Exponential offsets \f$ r_0 \f$ of the cutoff function \f$ F_c(R) \f$
    real(rp),dimension(:),allocatable :: r_0
    !> Step function offsets \f$ r_c \f$ of the cutoff function \f$ F_c(R) \f$
    real(rp),dimension(:),allocatable :: r_c
    !> Maximum step function offset
    real(rp) :: r_c_max
    !> Exponential lifetimes \f$ r_l \f$ of the cutoff function \f$ F_c(R) \f$
    real(rp),dimension(:),allocatable :: r_l
    !> @}

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: read_file_nrl
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type element_tb

  ! Constructor
  interface element_tb
    procedure :: constructor
  end interface element_tb

contains
  function constructor(u) result(obj)
    class(units),target,intent(in) :: u
    type(element_tb) :: obj

    ! Parent type constructor
    obj%element = element(u)
  end function constructor

  subroutine destructor(obj)
    type(element_tb) :: obj

    if(allocated(obj%filename)) deallocate(obj%filename)
    if(allocated(obj%nrl_type)) deallocate(obj%nrl_type)
    if(allocated(obj%p))        deallocate(obj%p)
    if(allocated(obj%r_0))      deallocate(obj%r_0)
    if(allocated(obj%r_c))      deallocate(obj%r_c)
    if(allocated(obj%r_l))      deallocate(obj%r_l)

  end subroutine destructor

  !> Check compatibility of TB parameter types
  subroutine check_nrl_type(type)
    character(len=*),dimension(:),intent(in) :: type
    integer :: ie

    do ie=1,size(type)-1
      if(type(ie) /= type(ie+1)) then
        write(error_unit,*) 'element_tb%check_type(): element_tb%nrl_type of &
         &elements ', ie, ' and ', ie+1, ' incompatible'
        error stop
      end if
    end do
  end subroutine check_nrl_type

  !> Read object from formatted TB nrlconstantopoulos file
  subroutine read_file_nrl(obj)
    class(element_tb),intent(inout) :: obj
    integer :: iostatus
    logical :: isopen
    character(len=7) :: type
    integer :: ie,ip

    write(output_unit,*) 'DEBUG == Entering element_tb & read_file_nrl'
    call TBKOSTER_flush(output_unit)

    if(allocated(obj%nrl_type)) deallocate(obj%nrl_type)
    allocate(obj%nrl_type(obj%ne))
    if(allocated(obj%p)) deallocate(obj%p)
    allocate(obj%p(obj%ne,np_max))
    if(allocated(obj%r_0)) deallocate(obj%r_0)
    allocate(obj%r_0(obj%ne))
    if(allocated(obj%r_c)) deallocate(obj%r_c)
    allocate(obj%r_c(obj%ne))
    if(allocated(obj%r_l)) deallocate(obj%r_l)
    allocate(obj%r_l(obj%ne))

    do ie=1,obj%ne
      inquire(unit=10,opened=isopen)
      if (isopen) then
        write(error_unit,'(a)') 'element_tb%read_file_nrl() : Unit 10 is already open'
        error stop
      else
        open(unit=10,file=obj%filename(ie),action='read',iostat=iostatus, status='old')
        if(iostatus /= 0) then
          write(error_unit,*) 'element_tb%read_file_nrl(): file ', obj%filename(ie), ' not found'
          error stop
        end if
      end if
    
      read(10,*) obj%nrl_type(ie)
      read(10,*)
      read(10,*)
      read(10,*) obj%r_0(ie), obj%r_c(ie), obj%r_l(ie)
      if(abs(obj%r_c(ie) - obj%r_0(ie) - 5*obj%r_l(ie)) > 0.00001_rp) then
        write(error_unit,*) 'element_tb%read_file_nrl(): bad cut-off'
        error stop
      end if
      read(10,*)
      read(10,*)
      read(10,*)
      do ip=1,np_max
        read(10,*) obj%p(ie,ip)
      end do

      close(unit=10)

      select case(trim(type))
      case('NN00000')
        obj%nrl_type(ie) = 'old'
      case('NN00001')
        obj%nrl_type(ie) = 'new'
      case('NN00002')
        obj%nrl_type(ie) = 'cyr'
      case('power')
        obj%nrl_type(ie) = 'pow'
      case default
        obj%nrl_type(ie) = 'cyr'
      end select
    end do

    obj%r_c_max = maxval(obj%r_c)

    ! Unit conversion from Rydberg atomic units to Hartree atomic units
    obj%p(:, 2:16) = 0.5_rp * obj%p(:, 2:16)
    do ip=18,54,4
      obj%p(:,ip:ip+2) = 0.5_rp * obj%p(:,ip:ip+2)
    end do

    write(output_unit,*) 'DEBUG == Exiting element_tb & read_file_nrl'
    call TBKOSTER_flush(output_unit)
  end subroutine read_file_nrl

  !> Read object in text format from file (default: 'in_element_tb.txt')
  subroutine read_txt(obj,file)
    class(element_tb),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    character(len=512) :: msg
    ! Namelist variables
    character(len=3) :: tb_type
    character(len=sl),dimension(:),allocatable :: filename
    ! Namelist
    namelist /element_tb/ tb_type, filename

    write(output_unit,*) 'DEBUG == Entering element_tb & read_txt'
    call TBKOSTER_flush(output_unit)

    tb_type='nrl' ! default type is nrl 
    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_element_tb.txt'
    end if

    ! Parent type procedure
    call obj%element%read_txt(file_rt)
    ! Derived type procedure
    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'element_tb%read_txt() : Unit 10 already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
      if(iostatus /= 0) then
        write(error_unit,*) 'element_tb%read_txt(): file ', file_rt, ' not found'
        error stop
      end if
    end if

    allocate(filename(obj%ne))
    read(unit=10,nml=element_tb,iostat=iostatus,iomsg=msg)
    obj%tb_type=tb_type
    close(unit=10)

    select case (lower(trim(tb_type)))
      case ("nrl")
        call move_alloc(filename,obj%filename)
        call obj%read_file_nrl()
        call check_nrl_type(obj%nrl_type)
      case ('wan')
        if (allocated(filename)) deallocate(filename)
        allocate(filename(1))
        filename(1)='hr.dat'
        call move_alloc(filename,obj%filename)
        write(output_unit,*) 'will read TB parameters from hr.dat file in build_b_r function'
      case ('mod') 
        if (allocated(filename)) deallocate(filename)
        allocate(filename(1))
        filename(1)='mod.dat'
        call move_alloc(filename,obj%filename)
        write(output_unit,*) 'will read TB parameters from mod.dat file in build_b_r function'
    end select

    write(output_unit,*) 'DEBUG == Exiting element_tb & read_txt'
    call TBKOSTER_flush(output_unit)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_element_tb.txt')
  subroutine write_txt(obj,file,unit)
    class(element_tb),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=len(obj%filename)),dimension(obj%ne) :: filename
    ! Namelist
    namelist /element_tb/ filename

    if (present(file)) then
      file_rt = file
    else
      file_rt = 'out_element_tb.txt'
    end if

    if (present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    ! Parent type procedure
    call obj%element%write_txt(file_rt,unit_rt)

    ! Derived type procedure
    if(.not. present(unit)) then
      !open(unit=unit_rt,file=file,action='write',status='old', &
      ! position='append')
      open(unit=unit_rt,file=file,action='write')
    end if

    filename = obj%filename
    
    write(unit_rt,nml=element_tb)

    if(.not. present(unit)) then
      call TBKOSTER_flush(unit_rt)
      close(unit_rt)
    else
      call TBKOSTER_flush(unit)
    end if

    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_element_tb.txt'), if tag (default: .true.) the namelist opening and
  !> closing tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(element_tb),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Local variables
    integer :: ie, ip

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_element_tb.txt'
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

    ! Parent type procedure
    if(present(property)) then
      if(present(unit)) then
        call obj%element%write_txt_formatted(file=file_rt, &
         property=property_rt,tag=tag_rt,unit=unit_rt)
      else
        call obj%element%write_txt_formatted(file=file_rt, &
         property=property_rt,tag=tag_rt)
       end if
    else !if(.not. present(property))
      if(present(unit)) then
        call obj%element%write_txt_formatted(file=file_rt,tag=tag_rt, &
         unit=unit_rt)
      else
        call obj%element%write_txt_formatted(file=file_rt,tag=tag_rt)
      end if
    end if

    ! Derived type procedure
    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if
    if(tag_rt) then
      write(unit_rt,'(a)') '&element_tb'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('filename')
        select case(obj%tb_type)
           case('nrl')
        do ie=1,obj%ne
        ! filenames with quotes
          write(unit_rt,'(a)') ' filename(' // int2str(ie) // ') = ' &
           // "'" // trim(obj%filename(ie)) // "'"
        ! filenames without quotes
        !          write(unit_rt,'(a)') ' filename(' // int2str(ie) // ') = ' &
        !           // trim(obj%filename(ie))
        end do
           case('mod','wan')
            write(unit_rt,'(a)') ' filename= '// "'" // trim(obj%filename(1)) // "'"
           end select
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
end module element_tb_mod
