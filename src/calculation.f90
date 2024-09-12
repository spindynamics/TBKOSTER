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
!  calculation.f90
!  TBKOSTER
module calculation_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_tb_mod
  use charge_mod
  use density_of_states_mod
  use magnetic_anisotropy_mod
  use band_structure_mod
  use element_tb_mod
  use energy_mod
  use forces_mod
  use hamiltonian_tb_mod
  use lattice_mod
  !use magnetic_molecular_dynamics_mod
  use mesh_mod
  use mixing_mod
  use molecular_dynamics_mod
  use precision_mod, only: rp
  use self_consistent_field_mod
  use spin_dynamics_mod
  use string_mod, only: TBKOSTER_flush, int2str, lower, sl
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'pre_processing', &
   'pre_processing_dir', &
   'processing', &
   'post_processing', &
   'post_processing_dir' &
   ]

  type,public :: calculation
    !> TBKOSTER executable directory
    character(len=sl) :: TBKOSTER_dir
    !> Pre-processing ; options:
    !>  'none' (default)
    !>  'atom': plot the atomic structure
    !>  'mesh': plot the reciprocal space mesh
    !>  'txt2xyz': convert an "atom.txt" file into an "atom.xyz" file
    character(len=7)  :: pre_processing
    !> Pre-processing directory
    character(len=sl) :: pre_processing_dir
    !> Processing ; options:
    !>  'none' (default)
    !>  'md' : molecular dynamics
    !>  'mmd': magnetic molecular dynamics
    !>  'scf': self-consistent field
    !>  'sd' : spin dynamics
    character(len=4)  :: processing
    !> Post-processing ; options:
    !>  'none' (default)
    !>  'band': calculate the band structure
    !>  'dos' : calculate the density of states
    !>  'forces' : calculate the forces on atoms
    !>  'mae' : calculate the magnetic anisotropy energy
    !>  'txt2xyz': convert an "atom.txt" file into an "atom.xyz" file
    character(len=7)  :: post_processing
    !> Post-processing directory
    character(len=sl) :: post_processing_dir

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: post_process_band
    procedure :: post_process_dos
    procedure :: post_process_mae
    procedure :: post_process_forces
    procedure :: post_process_txt2xyz
    procedure :: pre_process_txt2xyz
    procedure :: process
    procedure :: process_md
    !procedure :: process_mmd
    procedure :: process_scf
    procedure :: process_sd
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type calculation

  ! Constructor
  interface calculation
    procedure :: constructor
  end interface calculation

contains
  function constructor() result(obj)
    type(calculation) :: obj
    integer :: i

    call get_command_argument(0,obj%TBKOSTER_dir)
    if(obj%TBKOSTER_dir(1:1) /= '/') then
      write(error_unit,*) 'calculation%constructor(): TBKOSTER must be started &
       &using its absolute path'
      error stop
    end if
    i = index(obj%TBKOSTER_dir,'/',.true.)
    obj%TBKOSTER_dir = obj%TBKOSTER_dir(1:i-1)
  end function constructor

  subroutine check_post_processing(post_processing)
    character(len=*),intent(in) :: post_processing

    if(post_processing /= 'none' &
     .and. post_processing /= 'band' &
     .and. post_processing /= 'dos' &
     .and. post_processing /= 'forces' &
     .and. post_processing /= 'mae' &
     .and. post_processing /= 'txt2xyz') then
      write(error_unit,*) 'calculation%check_post_processing(): &
       &calculation%post_processing must be one of: ''none'', ''band'', &
       &''dos'', ''mae'',''forces'', ''txt2xyz'''
      error stop
    end if

  end subroutine check_post_processing

  subroutine check_pre_processing(pre_processing)
    character(len=*),intent(in) :: pre_processing

    if(pre_processing /= 'none' &
     .and. pre_processing /= 'atom' &
     .and. pre_processing /= 'mesh' &
     .and. pre_processing /= 'txt2xyz') then
      write(error_unit,*) 'calculation%check_pre_processing(): &
       &calculation%pre_processing must be one of: ''none'', ''atom'', &
       &''mesh'', ''txt2xyz'''
      error stop
    end if
  end subroutine check_pre_processing

  subroutine check_processing(processing)
    character(len=*),intent(in) :: processing

    if(processing /= 'none' &
     .and. processing /= 'md' &
     .and. processing /= 'mmd' &
     .and. processing /= 'scf' &
     .and. processing /= 'sd') then
      write(error_unit,*) 'calculation%check_processing(): &
       &calculation%processing must be one of: ''none'', ''md'', ''mmd'', &
       &''scf'', ''sd'''
      error stop
    end if
    
  end subroutine check_processing

  subroutine destructor(obj)
		type(calculation) :: obj
  end subroutine destructor

  subroutine initialize(pre_processing, pre_processing_dir, processing, &
   post_processing, post_processing_dir)
    character(len=7) ,intent(out) :: pre_processing
    character(len=sl),intent(out) :: pre_processing_dir
    character(len=4) ,intent(out) :: processing
    character(len=7) ,intent(out) :: post_processing
    character(len=sl),intent(out) :: post_processing_dir

    pre_processing      = 'none'
    pre_processing_dir  = '.'
    processing          = 'none'
    post_processing     = 'none'
    post_processing_dir = '.'
  end subroutine initialize

  subroutine post_process_band(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target :: atom_tb_obj
    type(charge),target :: charge_obj
    type(density_of_states),target :: dos_obj
    type(band_structure),target :: band_obj
    type(element_tb),target :: element_tb_obj
    type(energy),target :: energy_obj
    type(hamiltonian_tb),target :: hamiltonian_tb_obj
    type(lattice),target :: lattice_r, lattice_k
    type(mesh),target :: mesh_k
    type(self_consistent_field),target :: scf_obj
    type(units),target :: units_obj

    character(len=:),allocatable :: dir
    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence,file_existence2
    character(len=4) :: old_mesh_type
    integer :: unit
    logical :: isopen

    dir = trim(obj%post_processing_dir)
    ! Read input
    inquire(file=master_file,exist=file_existence)
    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
    end if
    mesh_k = mesh(lattice_k)
    call mesh_k%read_txt(file=dir // '/in_mesh.txt')

    charge_obj = charge(atom_tb_obj)
    call charge_obj%read_txt(file='out_charge.txt')
    if(file_existence) then
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt(file=master_file)
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt(file=master_file)
      dos_obj = density_of_states(energy_obj)
      call dos_obj%read_txt(file=master_file)
      band_obj = band_structure(energy_obj)
      call band_obj%read_txt(file=dir // '/in_band.txt')
    else
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt()
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt(file=dir // '/in_energy.txt')
      dos_obj = density_of_states(energy_obj)
      call dos_obj%read_txt(file=dir // '/in_dos.txt')
      band_obj = band_structure(energy_obj)
      call band_obj%read_txt(file=dir // '/in_band.txt')   
    end if
    !scf_obj = self_consistent_field(en=energy_obj)
    ! scf_obj = self_consistent_field(dos=dos_obj)
       scf_obj = self_consistent_field(band=band_obj)
    call scf_obj%initialize('nscf')
    ! Write mesh in list format
    old_mesh_type = mesh_k%type
    mesh_k%type = 'list'
    call mesh_k%dir2cart()
    call mesh_k%write_txt_formatted(file=dir // '/out_mesh.txt')
    call mesh_k%cart2dir()
    mesh_k%type = old_mesh_type
    ! Write out_hamiltonian_tb.txt
    call hamiltonian_tb_obj%write_txt_formatted(file=dir // '/out_hamiltonian_tb.txt',&
     property= [character(len=sl) :: 'nh','ns'])

    !unit = output_unit
    ! Open log
    unit = 11
    inquire(unit=unit,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'calculation%post_process_band() : Unit 11 is already open'
      error stop
    else
      ! Write log
      open(unit=11,file=dir // '/out_log.txt',action='write')
    end if

    call units_obj%write_txt_formatted(unit=unit)
    call element_tb_obj%write_txt_formatted(unit=unit)
    call lattice_r%write_txt_formatted(unit=unit)
    call atom_tb_obj%write_txt_formatted(unit=unit)
    call mesh_k%write_txt_formatted(unit=unit)
    call hamiltonian_tb_obj%write_txt_formatted(unit=unit)
    call energy_obj%write_txt_formatted(unit=unit)
    call band_obj%write_txt_formatted(unit=unit)
    ! Run
    write(unit,'(a)') ''
    write(unit,'(a)') 'calculation%run(): Post-processing band'
    call TBKOSTER_flush(unit)
    call hamiltonian_tb_obj%calculate_h_r()
    call hamiltonian_tb_obj%calculate_s_r()
    call scf_obj%run(unit,obj%post_processing)
    ! Close log
    close(unit)

    ! Write band structure
    call band_obj%write_txt_formatted(file=dir // '/out_band.txt', &
     property=[character(len=sl) :: 'en_k'])
  end subroutine post_process_band

  subroutine post_process_dos(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target :: atom_tb_obj
    type(charge),target :: charge_obj
    type(density_of_states),target :: dos_obj
    type(element_tb),target :: element_tb_obj
    type(energy),target :: energy_obj
    type(hamiltonian_tb),target :: hamiltonian_tb_obj
    type(lattice),target :: lattice_r, lattice_k
    type(mesh),target :: mesh_k
    type(self_consistent_field),target :: scf_obj
    type(units),target :: units_obj

    character(len=:),allocatable :: dir
    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence
    integer :: unit

    dir = trim(obj%post_processing_dir)

    ! Read input
    inquire(file=master_file,exist=file_existence)
    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
      !mesh_k = mesh(lattice_k)
      !call mesh_k%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
      !mesh_k = mesh(lattice_k)
      !call mesh_k%read_txt()
    end if
    mesh_k = mesh(lattice_k)
    call mesh_k%read_txt(file=dir // '/in_mesh.txt')
    charge_obj = charge(atom_tb_obj)
    call charge_obj%read_txt(file='out_charge.txt')
    if(file_existence) then
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt(file=master_file)
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt(file=dir // '/in_energy.txt')
      dos_obj = density_of_states(energy_obj)
      call dos_obj%read_txt(file=dir // '/in_dos.txt')
    else
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt()
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt()
      dos_obj = density_of_states(energy_obj)
      call dos_obj%read_txt()
    end if
    scf_obj = self_consistent_field(dos=dos_obj)
    call scf_obj%initialize('nscf')

    !unit = output_unit
    ! Open log
    unit = 11
    open(unit=unit,file=dir // '/out_log.txt',action='write')
    ! Write log
    call units_obj%write_txt_formatted(unit=unit)
    call element_tb_obj%write_txt_formatted(unit=unit)
    call lattice_r%write_txt_formatted(unit=unit)
    call atom_tb_obj%write_txt_formatted(unit=unit)
    call mesh_k%write_txt_formatted(unit=unit)
    call hamiltonian_tb_obj%write_txt_formatted(unit=unit)
    call energy_obj%write_txt_formatted(unit=unit)
    call dos_obj%write_txt_formatted(unit=unit)
    ! Write out_hamiltonian_tb.txt
    call hamiltonian_tb_obj%write_txt_formatted(file=dir // '/out_hamiltonian_tb.txt',&
     property= [character(len=sl) :: 'nh','ns'])
    ! Run
    write(unit,'(a)') ''
    write(unit,'(a)') 'calculation%run(): Post-processing dos'
    call TBKOSTER_flush(unit)
    call hamiltonian_tb_obj%calculate_h_r()
    call hamiltonian_tb_obj%calculate_s_r()
    call scf_obj%run(unit,obj%post_processing)
    ! Close log
    close(unit)

    ! Write density of states
    if(dos_obj%na_dos==0) then
      call dos_obj%write_txt_formatted(file=dir // '/out_dos.txt', &
       property=[character(len=sl) :: 'nen','en_min','en_max','na_dos','dos'])
    else !obj%na_dos>0
      call dos_obj%write_txt_formatted(file=dir // '/out_dos.txt', &
       property=[character(len=sl) :: 'nen','en_min','en_max','na_dos','ia', &
       'dos','dos_s','dos_p','dos_px','dos_py','dos_pz', &
       'dos_d','dos_dxy','dos_dyz','dos_dzx','dos_dx2y2','dos_dz2r2'])
    end if

    !    Plot density of states
    !    call execute_command_line('python ' // trim(obj%TBKOSTER_dir) &
    !     // '/../python/plot_dos.py')

    call execute_command_line(trim(obj%TBKOSTER_dir) // '/pdos.x ')
  end subroutine post_process_dos

  subroutine post_process_mae(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target :: atom_tb_obj
    type(charge),target :: charge_obj
    type(magnetic_anisotropy),target :: mae_obj
    type(element_tb),target :: element_tb_obj
    type(energy),target :: energy_obj
    type(hamiltonian_tb),target :: hamiltonian_tb_obj
    type(lattice),target :: lattice_r, lattice_k
    type(mesh),target :: mesh_k
    type(self_consistent_field),target :: scf_obj
    type(units),target :: units_obj
    character(len=:),allocatable :: dir
    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence
    integer :: unit, iangle
    integer :: ia, l
    real(rp),dimension(2) :: angle

    dir = trim(obj%post_processing_dir)

    ! Read input
    inquire(file=master_file,exist=file_existence)
    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
      !mesh_k = mesh(lattice_k)
      !call mesh_k%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
      !mesh_k = mesh(lattice_k)
      !call mesh_k%read_txt()
    end if
    mesh_k = mesh(lattice_k)
    call mesh_k%read_txt(file=dir // '/in_mesh.txt')
    charge_obj = charge(atom_tb_obj)
    charge_obj%a%ns=4
    call atom_tb_obj%calculate_spin()
    call charge_obj%initialize()
  !  call charge_obj%read_charge_col_to_ncol()
    if(file_existence) then
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt(file=master_file)
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt(file=dir // '/in_energy.txt')
      mae_obj = magnetic_anisotropy(energy_obj)
      call mae_obj%read_txt(file=dir // '/in_mae.txt')
    else
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt()
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt()
      mae_obj = magnetic_anisotropy(energy_obj)
      call mae_obj%read_txt()
    end if
    scf_obj = self_consistent_field(mae=mae_obj)
    call scf_obj%initialize('nscf')
    do iangle=1,mae_obj%nangle
    !unit = output_unit
    ! Open log
    unit = 11
    open(unit=unit,file=dir // '/out_log.txt',action='write')
    ! Write log
    call units_obj%write_txt_formatted(unit=unit)
    call element_tb_obj%write_txt_formatted(unit=unit)
    call lattice_r%write_txt_formatted(unit=unit)
    !call atom_tb_obj%write_txt_formatted(unit=unit)
    call atom_tb_obj%write_txt_formatted(unit=unit,property= [character(len=sl) :: 'ns','na','ia2ie','r'])
    call mesh_k%write_txt_formatted(unit=unit)
    call hamiltonian_tb_obj%write_txt_formatted(unit=unit)
    call energy_obj%write_txt_formatted(unit=unit)
    call mae_obj%write_txt_formatted(unit=unit)
    ! Write out_hamiltonian_tb.txt
    call hamiltonian_tb_obj%write_txt_formatted(file=dir // '/out_hamiltonian_tb.txt',&
     property= [character(len=sl) :: 'nh','ns'])
    ! Run
    !do iangle=1,mae_obj%nangle
    angle(:)=mae_obj%angle1(:)+(iangle-1)*(mae_obj%angle2(:)-mae_obj%angle1(:))/(mae_obj%nangle-1)
    call charge_obj%set_angle_charge_ncol(angle)
    call mae_obj%initialize()
    write(unit,'(a)') ''
    write(unit,'(a)') 'calculation%run(): Post-processing mae'//int2str(iangle)
    call TBKOSTER_flush(unit)
    call hamiltonian_tb_obj%calculate_h_r()
    call hamiltonian_tb_obj%calculate_s_r()
    call scf_obj%run(unit,obj%post_processing)
    ! Close log
    !close(unit)

    ! Write magnetic anisotropy
  
      call mae_obj%write_txt_formatted(file=dir // '/out_mae_'//int2str(iangle)//'.txt', &
       property=[character(len=sl) :: 'mae'])
    end do
    close(unit)

    call execute_command_line(trim(obj%TBKOSTER_dir) // '/mae.x ')
  end subroutine post_process_mae

  subroutine post_process_forces(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target :: atom_tb_obj
    type(charge),target :: charge_obj
    type(element_tb),target :: element_tb_obj
    type(energy),target :: energy_obj
    type(forces),target :: forces_obj
    type(hamiltonian_tb),target :: hamiltonian_tb_obj
    type(lattice),target :: lattice_r, lattice_k
    type(mesh),target :: mesh_k
    type(self_consistent_field),target :: scf_obj
    type(units),target :: units_obj

    character(len=:),allocatable :: dir
    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence
    integer :: unit
    logical :: isopen

    dir = trim(obj%post_processing_dir)

    ! Read input
    inquire(file=master_file,exist=file_existence)
    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
      mesh_k = mesh(lattice_k)
      call mesh_k%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
      mesh_k = mesh(lattice_k)
      call mesh_k%read_txt()
    end if
    charge_obj = charge(atom_tb_obj)
    call charge_obj%read_txt(file='out_charge.txt')
    if(file_existence) then
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt(file=master_file)
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt(file=master_file)
    else
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt()
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt()
    end if

    scf_obj = self_consistent_field(en=energy_obj)
    call scf_obj%initialize('fscf')

    forces_obj = forces(scf_obj)

    call forces_obj%read_txt(file=master_file)

    !unit = output_unit
    ! Open log
    unit = 11
    inquire(unit=unit,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'calculation%post_process_forces() : Unit 11 is already open'
      error stop
    else
      open(unit=unit,file=dir // '/out_log.txt',action='write')
    end if
    ! Write log
    call units_obj%write_txt_formatted(unit=unit)
    call element_tb_obj%write_txt_formatted(unit=unit)
    call lattice_r%write_txt_formatted(unit=unit)
    call atom_tb_obj%write_txt_formatted(unit=unit)
    call mesh_k%write_txt_formatted(unit=unit)
    call hamiltonian_tb_obj%write_txt_formatted(unit=unit)
    call energy_obj%write_txt_formatted(unit=unit)
    call forces_obj%write_txt_formatted(unit=unit)

    ! Run
    write(unit,'(a)') ''
    write(unit,'(a)') 'calculation%run(): Post-processing forces'
    call TBKOSTER_flush(unit)

    call hamiltonian_tb_obj%calculate_h_r()
    call hamiltonian_tb_obj%calculate_s_r()

    call scf_obj%run(unit)

    call forces_obj%calculate_forces()
    call forces_obj%print_forces(unit=unit)

    ! Close log
    close(unit)

  end subroutine post_process_forces

  subroutine post_process_txt2xyz(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target    :: atom_tb_obj
    type(element_tb),target :: element_tb_obj
    type(lattice),target    :: lattice_r, lattice_k
    type(units),target      :: units_obj
    character(len=:),allocatable :: dir, file
    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence
    integer :: ifile

    dir = trim(obj%post_processing_dir)

    ! Read input
    inquire(file=master_file,exist=file_existence)

    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
    end if

    ifile = 0
    file = dir // '/out_atom_tb_0'
    inquire(file=file // '.txt', exist=file_existence)
    do while(file_existence)
      ! Read txt
      call atom_tb_obj%read_txt(file=file // '.txt')
      ! Write xyz
      call atom_tb_obj%write_xyz(file=file // '.xyz')

      ifile = ifile+1
      deallocate(file)
      file = dir // '/out_atom_tb_' // int2str(ifile)
      inquire(file=file // '.txt', exist=file_existence)
    end do
  end subroutine post_process_txt2xyz

  subroutine pre_process_txt2xyz(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target :: atom_tb_obj
    type(element_tb),target :: element_tb_obj
    type(lattice),target :: lattice_r, lattice_k
    type(units),target :: units_obj

    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence

    ! Input
    inquire(file=master_file,exist=file_existence)
    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
    end if

    call atom_tb_obj%write_xyz(file=trim(obj%pre_processing_dir) &
     // '/in_atom.xyz')
  end subroutine pre_process_txt2xyz

  subroutine process(obj)
    class(calculation),intent(in) :: obj

    ! Pre-processing: txt-to-xyz format conversion
    select case(obj%pre_processing)
    case('txt2xyz')
      call obj%pre_process_txt2xyz
    end select

    ! Processing: molecular dynamics, magnetic molecular dynamics,
    ! self-consistent field, spin dynamics
    select case(obj%processing)
    case('md')
      call obj%process_md()
    case('mmd')
      !call obj%process_mmd()
    case('scf')
      call obj%process_scf()
    case('sd')
      call obj%process_sd()
    end select

    ! Post-processing: band structure, density of states, txt-to-xyz format
    ! conversion
    select case(obj%post_processing)
    case('band')
      call obj%post_process_band()
    case('dos')
      call obj%post_process_dos()
    case('mae')
      call obj%post_process_mae()
    case('forces')
      call obj%post_process_forces()
    case('txt2xyz')
      call obj%post_process_txt2xyz()
    end select
  end subroutine process

  subroutine process_md(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target :: atom_tb_obj
    type(charge),target :: charge_obj
    type(element_tb),target :: element_tb_obj
    type(energy),target :: energy_obj
    type(forces),target :: forces_obj
    type(hamiltonian_tb),target :: hamiltonian_tb_obj
    type(lattice),target :: lattice_r, lattice_k
    type(mesh),target :: mesh_k
    type(mixing),target :: mixing_obj
    type(molecular_dynamics),target :: md_obj
    type(self_consistent_field),target :: scf_obj
    type(units),target :: units_obj

    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence
    integer :: unit
    logical :: isopen
    ! Read input
    inquire(file=master_file,exist=file_existence)
    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
      mesh_k = mesh(lattice_k)
      call mesh_k%read_txt(file=master_file)
      charge_obj = charge(atom_tb_obj)
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt(file=master_file)
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt(file=master_file)
      mixing_obj = mixing(charge_obj)
      call mixing_obj%read_txt(file=master_file)
      scf_obj = self_consistent_field(en=energy_obj,mx=mixing_obj)
      call scf_obj%read_txt(file=master_file)
      forces_obj = forces(scf=scf_obj)
      md_obj = molecular_dynamics(forces_obj)
      call md_obj%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
      mesh_k = mesh(lattice_k)
      call mesh_k%read_txt()
      charge_obj = charge(atom_tb_obj)
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt()
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt()
      mixing_obj = mixing(charge_obj)
      call mixing_obj%read_txt()
      scf_obj = self_consistent_field(en=energy_obj,mx=mixing_obj)
      call scf_obj%read_txt()
      forces_obj = forces(scf_obj)
      call forces_obj%read_txt()
      md_obj = molecular_dynamics(forces_obj)
      call md_obj%read_txt()
    end if

    !unit = output_unit
    ! Open log
    unit = 11
    inquire (unit=unit, opened = isopen)
    if (isopen) then
      write(error_unit,'(a)') 'calculation%process_md(): Unit 11 already open'
      error stop
    else
      open(unit=unit,file='out_log.txt',action='write')
    end if

    ! Run
    write(unit,'(a)') ''
    write(unit,'(a)') 'calculation%run(): Processing md'
    call hamiltonian_tb_obj%calculate_h_r()
    call hamiltonian_tb_obj%calculate_s_r()
    call md_obj%integrate(unit=unit)

    ! Close log
    close(unit)
  end subroutine process_md

  subroutine process_scf(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target :: atom_tb_obj
    type(charge),target :: charge_obj
    type(element_tb),target :: element_tb_obj
    type(energy),target :: energy_obj
    type(hamiltonian_tb),target :: hamiltonian_tb_obj
    type(lattice),target :: lattice_r, lattice_k
    type(mesh),target :: mesh_k
    type(mixing),target :: mixing_obj
    type(self_consistent_field),target :: scf_obj
    type(units),target :: units_obj

    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence, isopen
    character(len=4) :: old_mesh_type
    integer :: unit

    ! Read input
    inquire(file=master_file,exist=file_existence)
    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
      mesh_k = mesh(lattice_k)
      call mesh_k%read_txt(file=master_file)
      charge_obj = charge(atom_tb_obj)
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt(file=master_file)
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt(file=master_file)
      mixing_obj = mixing(charge_obj)
      call mixing_obj%read_txt(file=master_file)
      scf_obj = self_consistent_field(en=energy_obj,mx=mixing_obj)
      call scf_obj%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
      mesh_k = mesh(lattice_k)
      call mesh_k%read_txt()
      charge_obj = charge(atom_tb_obj)
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt()
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt()
      mixing_obj = mixing(charge_obj)
      call mixing_obj%read_txt()
      scf_obj = self_consistent_field(en=energy_obj,mx=mixing_obj)
      call scf_obj%read_txt()
    end if
    inquire(file='in_charge.txt',exist=file_existence)
    if(file_existence) then
      call charge_obj%read_txt()
    else
      call charge_obj%calculate_charge_in()
    end if

    ! Write units
    ! call units_obj%write_txt_formatted()
    ! Write mesh in list format
    old_mesh_type = mesh_k%type
    mesh_k%type = 'list'
    call mesh_k%dir2cart()
    call mesh_k%write_txt_formatted()
    call mesh_k%cart2dir()
    mesh_k%type = old_mesh_type
    ! Write out_hamiltonian_tb.txt
    ! call hamiltonian_tb_obj%write_txt_formatted(property= &
    ! [character(len=sl) :: 'nh','ns'])

    !unit = output_unit
    ! Open log
    unit = 11
    inquire (unit=unit, opened = isopen)
    if (isopen) then
      write(error_unit,'(a)') 'calculation%process_scf(): Unit 11 already open'
      error stop
    else
      open(unit=unit,file='out_log.txt',action='write')
    end if

    ! Write log
    call units_obj%write_txt_formatted(unit=unit)
    call element_tb_obj%write_txt_formatted(unit=unit)
    call lattice_r%write_txt_formatted(unit=unit)
    call atom_tb_obj%write_txt_formatted(unit=unit)
    call mesh_k%write_txt_formatted(unit=unit)
    call hamiltonian_tb_obj%write_txt_formatted(unit=unit)
    call energy_obj%write_txt_formatted(unit=unit)
    call mixing_obj%write_txt_formatted(unit=unit)
    call scf_obj%write_txt_formatted(unit=unit)
    call atom_tb_obj%write_txt_formatted(unit=unit,property= &
     [character(len=sl) :: 'nel','nn'])

    ! Run
    write(unit,'(a)') ''
    write(unit,'(a)') 'calculation%run(): Processing scf'
    call TBKOSTER_flush(unit)

    call hamiltonian_tb_obj%calculate_h_r()
    call hamiltonian_tb_obj%calculate_s_r()

    call scf_obj%run(unit)

    ! Close log
    inquire (unit=unit, opened = isopen)
    if (isopen) then
      close(unit)
    else
      write(error_unit,'(a)') 'calculation%process_scf(): Unit 11 already closed'
      error stop
    endif

    ! Write Fermi level en_f and total energy en_out
    call energy_obj%write_txt_formatted(property=[character(len=sl) :: 'en_f','en_out'])
  end subroutine process_scf

  subroutine process_sd(obj)
    class(calculation),intent(in) :: obj
    type(atom_tb),target :: atom_tb_obj
    type(charge),target :: charge_obj
    type(element_tb),target :: element_tb_obj
    type(energy),target :: energy_obj
    type(hamiltonian_tb),target :: hamiltonian_tb_obj
    type(lattice),target :: lattice_r, lattice_k
    type(mesh),target :: mesh_k
    type(mixing),target :: mixing_obj
    type(self_consistent_field),target :: scf_obj
    type(spin_dynamics),target :: sd_obj
    type(units),target :: units_obj

    character(len=*),parameter :: master_file = 'in_master.txt'
    logical :: file_existence
    integer :: unit
    logical :: isopen
    character(len=4) :: old_mesh_type

    ! Read input
    inquire(file=master_file,exist=file_existence)
    if(file_existence) then
      call units_obj%read_txt(file=master_file)
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt(file=master_file)
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt(file=master_file)
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt(file=master_file)
      mesh_k = mesh(lattice_k)
      call mesh_k%read_txt(file=master_file)
      charge_obj = charge(atom_tb_obj)
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt(file=master_file)
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt(file=master_file)
      mixing_obj = mixing(charge_obj)
      call mixing_obj%read_txt(file=master_file)
      scf_obj = self_consistent_field(en=energy_obj,mx=mixing_obj)
      call scf_obj%read_txt(file=master_file)
      sd_obj = spin_dynamics(scf_obj)
      call sd_obj%read_txt(file=master_file)
    else
      call units_obj%read_txt()
      element_tb_obj = element_tb(units_obj)
      call element_tb_obj%read_txt()
      lattice_r = lattice(units_obj)
      call lattice_r%read_txt()
      lattice_k = lattice_r%construct_reciprocal()
      atom_tb_obj = atom_tb(element_tb_obj,lattice_r,lattice_k)
      call atom_tb_obj%read_txt()
      mesh_k = mesh(lattice_k)
      call mesh_k%read_txt()
      charge_obj = charge(atom_tb_obj)
      hamiltonian_tb_obj = hamiltonian_tb(atom_tb_obj,charge_obj,mesh_k)
      call hamiltonian_tb_obj%read_txt()
      energy_obj = energy(hamiltonian_tb_obj)
      call energy_obj%read_txt()
      mixing_obj = mixing(charge_obj)
      call mixing_obj%read_txt()
      scf_obj = self_consistent_field(en=energy_obj,mx=mixing_obj)
      call scf_obj%read_txt()
      sd_obj = spin_dynamics(scf_obj)
      call sd_obj%read_txt()
    end if
    inquire(file='in_charge.txt',exist=file_existence)
    if(file_existence) then
      call charge_obj%read_txt()
    else
      call charge_obj%calculate_charge_in()
    end if

    ! Write units
    call units_obj%write_txt_formatted()

    ! Write mesh in list format
    old_mesh_type = mesh_k%type
    mesh_k%type = 'list'
    call mesh_k%dir2cart()
    call mesh_k%write_txt_formatted()
    call mesh_k%cart2dir()
    mesh_k%type = old_mesh_type

    !unit = output_unit
    ! Open log
    unit = 11
    inquire (unit=unit, opened = isopen)
    if (isopen) then
      write(error_unit,'(a)') 'calculation%process_sd(): Unit 11 already open'
      error stop
    else
      open(unit=unit,file='out_log.txt',action='write')
    end if
    ! Write log
    call units_obj%write_txt_formatted(unit=unit)
    call element_tb_obj%write_txt_formatted(unit=unit)
    call lattice_r%write_txt_formatted(unit=unit)
    call atom_tb_obj%write_txt_formatted(unit=unit)
    call mesh_k%write_txt_formatted(unit=unit)
    call hamiltonian_tb_obj%write_txt_formatted(unit=unit)
    call energy_obj%write_txt_formatted(unit=unit)
    call mixing_obj%write_txt_formatted(unit=unit)
    call scf_obj%write_txt_formatted(unit=unit)
    call sd_obj%write_txt_formatted(unit=unit)
    call atom_tb_obj%write_txt_formatted(unit=unit,property= &
     [character(len=sl) :: 'nel','nn'])
    ! Run
    write(unit,'(a)') ''
    write(unit,'(a)') 'calculation%run(): Processing sd'
    call hamiltonian_tb_obj%calculate_h_r()
    call hamiltonian_tb_obj%calculate_s_r()
    call sd_obj%integrate(unit)

    ! Close log
    inquire (unit=unit, opened = isopen)
    if (isopen) then
      close(unit)
    else
      write(error_unit,'(a)') 'calculation%process_sd(): Unit 11 already closed'
      error stop
    endif

  end subroutine process_sd

  !> Read object in text format from file (default: 'in_calculation.txt')
  subroutine read_txt(obj,file)
    class(calculation),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    ! Namelist variables
    character(len=7)  :: pre_processing
    character(len=sl) :: pre_processing_dir
    character(len=4)  :: processing
    character(len=7)  :: post_processing
    character(len=sl) :: post_processing_dir
    ! Namelist
    namelist /calculation/ pre_processing, pre_processing_dir, processing, &
     post_processing, post_processing_dir
    ! Local variables
    integer :: iostatus
    logical :: isopen

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_calculation.txt'
    end if

    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'calculation%read_txt(): Unit 10 already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if

    if(iostatus /= 0) then
      write(error_unit,*) 'calculation%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    call initialize(pre_processing,pre_processing_dir,processing, &
     post_processing,post_processing_dir)
    read(10,nml=calculation)
    ! Pre-processing
    pre_processing = lower(pre_processing)
    call check_pre_processing(trim(pre_processing))
    if(pre_processing /= 'none') then
      pre_processing_dir = './' // pre_processing
    end if
    ! Processing
    processing = lower(processing)
    call check_processing(trim(processing))
    ! Post-processing
    post_processing = lower(post_processing)
    call check_post_processing(trim(post_processing))
    if(post_processing /= 'none') then
      post_processing_dir = './' // post_processing
    end if

    obj%pre_processing = pre_processing
    obj%processing = processing
    obj%post_processing = post_processing

    rewind(10)
    read(10,nml=calculation)

    obj%pre_processing_dir = pre_processing_dir
    obj%post_processing_dir = post_processing_dir

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_calculation.txt')
  subroutine write_txt(obj,file,unit)
    class(calculation),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=len(obj%pre_processing)) :: pre_processing
    character(len=len(obj%pre_processing_dir)) :: pre_processing_dir
    character(len=len(obj%processing)) :: processing
    character(len=len(obj%post_processing)) :: post_processing
    character(len=len(obj%post_processing_dir)) :: post_processing_dir
    ! Namelist
    namelist /calculation/ pre_processing, pre_processing_dir, processing, &
     post_processing, post_processing_dir

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_calculation.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    pre_processing = obj%pre_processing
    pre_processing_dir = obj%pre_processing_dir
    processing = obj%processing
    post_processing = obj%post_processing
    post_processing_dir = obj%post_processing_dir

    write(unit_rt,nml=calculation)

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
  !> 'out_calculation.txt'), if tag (default: .true.) the namelist opening and
  !> closing tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(calculation),intent(in) :: obj
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
      file_rt = 'out_calculation.txt'
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
      write(unit_rt,'(a)') '&calculation'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('TBKOSTER_dir')
        write(unit_rt,'(a)') ' TBKOSTER_dir = ''' &
         // trim(obj%TBKOSTER_dir) // ''''
      case('post_processing')
        write(unit_rt,'(a)') ' post_processing = ''' &
         // trim(obj%post_processing) // ''''
      case('post_processing_dir')
        write(unit_rt,'(a)') ' post_processing_dir = ''' &
         // trim(obj%post_processing_dir) // ''''
      case('pre_processing')
        write(unit_rt,'(a)') ' pre_processing = ''' &
         // trim(obj%pre_processing) // ''''
      case('pre_processing_dir')
        write(unit_rt,'(a)') ' pre_processing_dir = ''' &
         // trim(obj%pre_processing_dir) // ''''
      case('processing')
        write(unit_rt,'(a)') ' processing = ''' // trim(obj%processing) // ''''
      end select
    end do

    if(tag_rt) then
      write(unit_rt,'(a)') ' /'
    end if

    if(.not. present(unit)) then
      call TBKOSTER_flush(unit_rt)
      close(unit_rt)
    else
      call TBKOSTER_flush(unit)
    end if

    !deallocate(file_rt,property_rt)
  end subroutine write_txt_formatted
end module calculation_mod
