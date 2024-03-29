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
!  self_consistent_field.f90
!  TBKOSTER
module self_consistent_field_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_mod
  use charge_mod
  use density_of_states_mod
  use band_structure_mod
  use energy_mod
  use hamiltonian_tb_mod
  use mesh_mod
  use mixing_mod
#if defined(OpenMP_Fortran_FOUND)
  use omp_lib
#endif
  use precision_mod, only: rp
  use string_mod, only: TBKOSTER_flush, int2str, log2str, lower, real2str, sl
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'ni_min', &
   'ni_max', &
   'delta_en', &
   'delta_q', &
   'verbose' &
   ]

  type,public :: self_consistent_field
    !> Units
    class(units),pointer :: u
    !> Atoms
    class(atom),pointer :: a
    !> Charges
    class(charge),pointer :: q
    !> Reciprocal space mesh
    class(mesh),pointer :: k
    !> Hamiltonian
    class(hamiltonian_tb),pointer :: h
    !> Energy
    class(energy),pointer :: en
    !> Mixing
    class(mixing),pointer :: mx
    !> Density of states
    class(density_of_states),pointer :: dos
    !> Band structure
    class(band_structure),pointer :: band   
    !> Logical forces (default: .false.)
    logical :: forces_logical
    !> Minimum number of iterations (default: 2)
    integer :: ni_min
    !> Maximum number of iterations (default: 50)
    integer :: ni_max
    !> Energy convergence criterion (default: 0.0)
    real(rp) :: delta_en
    !> Charge convergence criterion (default: 0.0)
    real(rp) :: delta_q
    !> Output verbosity flag (default: false)
    logical :: verbose
    !> Convergence flag
    logical :: converged

    ! SUBSYSTEM
    ! Define two subsystems which reference energy can be shifted to avoid large
    ! differences in Fermi energy of two subsystems (for example: molecule on
    ! surface) and help with local charge neutrality (LCN) convergence
    !
    ! Site index such that:
    ! if i<=nSubsystem i\in Sys 1
    ! if i>nSubsystem  i\in Sys 2
    integer  :: nSubsystem
    ! Potential such that
    ! if i\in Sys 1 Vshift=Vshift1
    ! if i\in Sys 2 Vshift=Vshift2
    real(rp) :: Vshift1,Vshift2

  contains
    ! Procedures
    procedure :: initialize
    procedure :: read_txt
    procedure :: run
    procedure :: update_m
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type self_consistent_field

  ! Constructor
  interface self_consistent_field
    procedure :: constructor
  end interface self_consistent_field

contains
  function constructor(en,mx,dos,band) result(obj)
    class(energy),target,intent(in),optional :: en
    class(mixing),target,intent(in),optional :: mx
    class(density_of_states),target,intent(in),optional :: dos
    class(band_structure),target,intent(in),optional :: band
    type(self_consistent_field) :: obj

    if(present(en)) then
      obj%u => en%u
      obj%a => en%a
      obj%q => en%q
      obj%k => en%k
      obj%h => en%h
      obj%en => en
    elseif(present(dos)) then
      obj%u => dos%u
      obj%a => dos%a
      obj%q => dos%en%q
      obj%k => dos%k
      obj%h => dos%h
      obj%en => dos%en
      obj%dos => dos
    elseif(present(band)) then
      obj%u => band%u
      obj%a => band%a
      obj%q => band%en%q
      obj%k => band%k
      obj%h => band%h
      obj%en => band%en
      obj%band => band
    else
      write(error_unit,*) 'self_consistent_field%constructor(): one of &
       & arguments en or dos must be specified'
      error stop
    end if

    if(present(mx)) then
      obj%mx => mx
    end if
  end function constructor

  subroutine initialize(obj,type)
    class(self_consistent_field),intent(inout) :: obj
    character(len=*),intent(in) :: type
    select case(lower(trim(type)))
    case('nscf')
      call initialize_nscf(obj%ni_min,obj%ni_max,obj%delta_en,obj%delta_q,obj%forces_logical)
    case('scf')
      call initialize_scf(obj%ni_min,obj%ni_max,obj%delta_en,obj%delta_q,obj%forces_logical)
    case('fscf')
      call initialize_fscf(obj%ni_min,obj%ni_max,obj%delta_en,obj%delta_q,obj%forces_logical)
    end select
  end subroutine initialize

  subroutine initialize_nscf(ni_min,ni_max,delta_en,delta_q,forces_logical)
    integer,intent(out) :: ni_min,ni_max
    real(rp),intent(out) :: delta_en,delta_q
    logical :: forces_logical

    ni_min = 1
    ni_max = 1
    delta_en = 0.0_rp
    delta_q  = 0.0_rp
    forces_logical = .false.
  end subroutine initialize_nscf

  subroutine initialize_scf(ni_min,ni_max,delta_en,delta_q,forces_logical)
    integer,intent(out) :: ni_min,ni_max
    real(rp),intent(out) :: delta_en,delta_q
    logical :: forces_logical

    ni_min = 2
    ni_max = 50
    delta_en = 0.0_rp
    delta_q  = 0.0_rp
    forces_logical = .false.
  end subroutine initialize_scf

  subroutine initialize_fscf(ni_min,ni_max,delta_en,delta_q,forces_logical)
    integer,intent(out) :: ni_min,ni_max
    real(rp),intent(out) :: delta_en,delta_q
    logical :: forces_logical

    ni_min = 1
    ni_max = 1
    delta_en = 0.0_rp
    delta_q  = 0.0_rp
    forces_logical = .true.
  end subroutine initialize_fscf

  !> Read object in text format from file (default: 'in_scf.txt')
  subroutine read_txt(obj,file)
    class(self_consistent_field),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    integer :: ni_min,ni_max
    real(rp) :: delta_en,delta_q
    logical :: verbose, forces_logical
    ! Namelist
    namelist /scf/ ni_min, ni_max, delta_en, delta_q, verbose, forces_logical

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_scf.txt'
    end if

    inquire(unit=10, opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'self_consistent_field%read_txt() : Unit 10 is &
      already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'self_consistent_field%read_txt(): file ', file_rt, &
       ' not found'
      error stop
    end if

    call initialize_scf(ni_min,ni_max,delta_en,delta_q,forces_logical)
    verbose = .false.
    read(10,nml=scf)

    obj%ni_min = ni_min
    obj%ni_max = ni_max
    obj%delta_en = delta_en * obj%u%convert_energy('to','hau')
    obj%delta_q  = delta_q
    obj%verbose = verbose
    obj%forces_logical = .false.
    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  subroutine run(obj,unit,post_processing)
    class(self_consistent_field),intent(inout) :: obj
    integer,intent(in),optional :: unit
    character(len=7),intent(in),optional  :: post_processing
    integer                     :: unit_rt
    ! Iteration index
    integer :: i
    ! Reciprocal space vector index
    integer :: ik
    ! Spin loop index
    integer :: isl
    ! Clock
    real(rp) :: time
    integer :: icount0, icount1, icount_rate, icount_max
    ! Eigenvalues
    !real(rp),dimension(obj%h%nh) :: w_k
    real(rp),dimension(:),allocatable :: w_k
    ! Eigenvectors
    !complex(rp),dimension(2,obj%h%nh,obj%h%nh) :: v_k
    complex(rp),dimension(:,:,:),allocatable :: v_k
    ! Output charge
    type(charge) :: q_out_k,om_k
    ! DOS
    !type(density_of_states) :: dos_k
    ! Convergence
    logical :: converged_en, converged_q

    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = output_unit
    end if

    !===========================================================================
    !                           WRITE INITIAL CHARGE
    !===========================================================================
    if(obj%verbose) then
      call obj%q%write_mulliken_charge_analysis(unit_rt,'in')
      !      call obj%q%write_atom_mag_on_the_fly(file='scf/out_atom_tb_0.txt',intent='in')
       call obj%a%write_txt_formatted(file='scf/out_atom_tb_0.txt')
    end if

    !===========================================================================
    !                              BEGIN SCF LOOP
    !===========================================================================
    obj%converged = .false.
    converged_en = .false.
    do i=1,obj%ni_max ! begin self consistent loop
      write(unit_rt,'(a)') '===== self_consistent_field%run(): iteration ' &
       // int2str(i) // ' ====='
      call TBKOSTER_flush(unit_rt)
      !=========================================================================
      !                        BUILD RENORMALIZATION
      !=========================================================================
      call obj%h%calculate_delta_h_eei()
      !call obj%h%calculate_delta_h_subsystem()
      call obj%h%calculate_delta_v_lcn()
      call obj%h%calculate_delta_v_pen()

      !=========================================================================
      !                     BEGIN K-POINT AND SPIN LOOPS
      !=========================================================================
      call system_clock(icount0,icount_rate,icount_max)
      if (.not.allocated(w_k)) allocate(w_k(obj%h%nh))
#if defined(OpenMP_Fortran_FOUND)
!$OMP PARALLEL DO &
!$OMP PRIVATE(ik,isl,w_k)
#endif
    do ik=1,obj%k%nx
      do isl=1,obj%a%nsl
        w_k = obj%h%build_w_k(ik, isl)
        call obj%en%save_en_k(ik, isl, w_k)
      end do ! end of isl
    end do ! end of ik
#if defined(OpenMP_Fortran_FOUND)
!$OMP END PARALLEL DO
#endif
    if (allocated(w_k)) deallocate(w_k)
    call system_clock(icount1,icount_rate,icount_max)
    time = real((icount1-icount0))/real(icount_rate)
    icount0 = icount1
    if(obj%verbose) then
      write(unit_rt,'(a)') 'self_consistent_field%run(): &
        &Time for k-loop 1 is ' // real2str(time) // ' s'
      call TBKOSTER_flush(unit_rt)
    end if
    !=========================================================================
    !                      END K-POINT AND SPIN LOOPS
    !=========================================================================

    !=========================================================================
    !                 SORT EIGENVALUES IN "INCREASING" ORDER
    !                 CALCULATE FERMI LEVEL AND BAND ENERGY
    !=========================================================================
    call obj%en%sort_en_k()
    call obj%en%calculate_en_f()
    call obj%en%calculate_en_band()
    if(obj%verbose) then
      call obj%en%write_txt_formatted(unit=unit_rt,property= &
       [character(len=sl) :: 'en_k_min','en_k_max','en_f'],tag=.false.)
       call TBKOSTER_flush(unit_rt)
    end if

    call obj%q%nullify_charge_out()
    call obj%q%nullify_orbital_moment()
    !=========================================================================
    !                    BEGIN K-POINT AND SPIN LOOPS
    !=========================================================================
    call system_clock(icount0,icount_rate,icount_max)
    if (.not.allocated(v_k)) allocate(v_k(2,obj%h%nh,obj%h%nh))
#if defined(OpenMP_Fortran_FOUND)
!$OMP PARALLEL DO &
!$OMP PRIVATE(ik,isl,v_k,q_out_k,om_k)
#endif
    do ik=1,obj%k%nx
      if(obj%ni_max>1) then
        q_out_k = charge(obj%a)
        om_k=charge(obj%a)
        if(obj%h%e_e_interaction == 'ujb') then
          q_out_k%rho_net_out_diagonal = .false.
        end if
      end if
      do isl=1,obj%a%nsl
        v_k = obj%h%build_v_k(ik, isl)
        if(obj%ni_max>1) then
          call q_out_k%calculate_charge_out(ik, isl, obj%k%nx, obj%k%w(ik), &
            obj%h%nh, obj%h%iaos2ih, obj%en%f_k, v_k)  
          call om_k%calculate_orbital_moment(ik,obj%k%nx,obj%h%iaos2ih,obj%k%w(ik),obj%en%f_k,v_k)
#if defined(OpenMP_Fortran_FOUND)
!$OMP CRITICAL
#endif
          call obj%q%add_charge_out_k(q_out_k)
          call obj%q%add_orbital_moment_k(om_k)
#if defined(OpenMP_Fortran_FOUND)
!$OMP END CRITICAL
#endif
        else if (obj%forces_logical) then

        else !if (obj%ni_max<=1)
#if defined(OpenMP_Fortran_FOUND)
!$OMP CRITICAL
#endif
        ! entering here when post-proc band and/or dos is activated
        if(present(post_processing)) then
          select case(post_processing)
            case('dos') 
              call obj%dos%add_dos_k(ik, isl)
              call obj%dos%add_dos_local_k(ik, isl, v_k)
            case('band') 
              if(TRIM(obj%band%proj)=='site') then
                call obj%band%save_proj_band_site(ik, isl, v_k)
              elseif(TRIM(obj%band%proj)=='spin') then
                call obj%band%save_proj_band_spin(ik, isl, v_k)
              elseif(TRIM(obj%band%proj)=='orbit') then
                call obj%band%save_proj_band_orbit(ik, isl, v_k)
              elseif(TRIM(obj%band%proj)=='spin,orbit') then
                call obj%band%save_proj_band_spin(ik, isl, v_k)
                call obj%band%save_proj_band_orbit(ik, isl, v_k) 
              endif
          end select
        end if
          
#if defined(OpenMP_Fortran_FOUND)
!$OMP END CRITICAL
#endif
        end if
      end do ! end of isl
    end do ! end of ik
#if defined(OpenMP_Fortran_FOUND)
!$OMP END PARALLEL DO
#endif
    if (allocated(v_k)) deallocate(v_k)
    call system_clock(icount1,icount_rate,icount_max)
    time = real((icount1-icount0))/real(icount_rate)
    icount0 = icount1
    if(obj%verbose) then
      write(unit_rt,'(a)') 'self_consistent_field%run(): &
        &Time for k-loop 2 is ' // real2str(time) // ' s'
      call TBKOSTER_flush(unit_rt)
    end if
    !=========================================================================
    !                      END K-POINT AND SPIN LOOPS
    !=========================================================================
    if(obj%ni_max>1) then
      !=======================================================================
      !           CALCULATE DOUBLE COUNTING TERMS AND TOTAL ENERGY
      !=======================================================================
      call obj%en%calculate_en_band()
      call obj%en%calculate_en_dc_eei()
      call obj%en%calculate_en_dc_lcn()
      call obj%en%calculate_en_dc_pen()
      call obj%en%calculate_en()
      if(obj%verbose) then
        call obj%en%write_txt_formatted(unit=unit_rt,property= &
         [character(len=sl) :: 'en_band','en_dc_eei','en_dc_lcn', &
         'en_dc_pen','en_out'],tag=.false.)
      end if

      !=======================================================================
      !                            TEST CONVERGENCE
      !=======================================================================
      ! Energy convergence
      if(i >= 2) then
        call obj%en%calculate_delta_en()
        if(obj%verbose) then
          call obj%en%write_txt_formatted(unit=unit_rt,property=&
            [character(len=sl) :: 'delta_en'],tag=.false.)
        end if
        converged_en = obj%en%is_converged(obj%delta_en)
      end if

      ! Charge convergence
      call obj%q%calculate_delta_q()
      if(obj%verbose) then
        call obj%q%write_txt_formatted(unit=unit_rt,property=&
           [character(len=sl) :: 'delta_q_mul','delta_rho_net'],tag=.false.)
      end if
      converged_q = obj%q%is_converged(obj%delta_q)

      ! Global convergence
      if(i >= obj%ni_min .and. converged_en .and. converged_q) then
        obj%converged = .true.
        exit
      end if

      !=======================================================================
      !                                  MIX
      !=======================================================================
      call obj%mx%mix(i)
      call obj%q%write_txt_formatted()
      if(obj%verbose) then
        call obj%q%write_mulliken_charge_analysis(unit_rt,'out')
        call obj%update_m
        !call obj%q%write_atom_mag_on_the_fly(file='scf/out_atom_tb_' // &int2str(i) // '.txt',intent='out')
        call obj%a%write_txt_formatted(file='scf/out_atom_tb_' // int2str(i) &
         // '.txt')
      end if
        call obj%q%set_q_in_to_q_out()
      end if
    end do ! end of self consistent loop

    if(obj%ni_max>1) then
      if(obj%converged) then
        write(unit_rt,'(a)') 'self_consistent_field%run(): SCF loop converged'
      else
        write(unit_rt,'(a)') 'self_consistent_field%run(): SCF loop did not &
         &converge'
      end if
    end if
    if(obj%verbose.and.obj%a%ns==4) then
      call obj%q%write_orbital_moment_analysis(unit_rt)
    endif
 !   write(*,*) obj%en%w_en_band_local
    call TBKOSTER_flush(unit_rt)
  end subroutine run

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_scf.txt')
  subroutine write_txt(obj,file,unit)
    class(self_consistent_field),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    integer :: ni_min,ni_max
    real(rp) :: delta_en,delta_q
    logical :: verbose
    ! Namelist
    namelist /scf/ ni_min, ni_max, delta_en, delta_q, verbose

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_scf.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    ni_min = obj%ni_min
    ni_max = obj%ni_max
    delta_en = obj%delta_en * obj%u%convert_energy('from','hau')
    delta_q = obj%delta_q
    verbose = obj%verbose

    write(unit_rt,nml=scf)
    call TBKOSTER_flush(unit_rt)
  
    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_scf.txt'), if tag (default: .true.) the namelist opening and closing
  !> tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(self_consistent_field),intent(in) :: obj
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
      file_rt = 'out_scf.txt'
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
      write(unit_rt,'(a)') '&scf'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('ni_min')
        write(unit_rt,'(a)') ' ni_min = ' // int2str(obj%ni_min)
      case('ni_max')
        write(unit_rt,'(a)') ' ni_max = ' // int2str(obj%ni_max)
      case('delta_en')
        write(unit_rt,'(a)') ' delta_en = ' // real2str(obj%delta_en &
         * obj%u%convert_energy('from','hau'))
      case('delta_q')
        write(unit_rt,'(a)') ' delta_q = ' // real2str(obj%delta_q)
      case('verbose')
        write(unit_rt,'(a)') ' verbose = ' // log2str(obj%verbose)
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

subroutine update_m(obj)
  use math_mod, only:  cart2sph
  class(self_consistent_field),intent(inout) :: obj
  ! LOCAL
  integer :: ia
  real(rp) :: m_s, m_p, m_d, m
  real(rp),dimension(3) :: m_sph, m_cart

  !write(output_unit,*) "====> Entering update_m"

  select case(obj%a%ns)
  case(1)
    do ia=1,obj%a%na
       obj%a%m(ia,:)=0_rp  ! update m
    end do
  case(2)
    do ia=1,obj%a%na
      m_s = obj%q%q_mul_out(ia,1,0)-obj%q%q_mul_out(ia,1,1)
      m_p = obj%q%q_mul_out(ia,2,0)-obj%q%q_mul_out(ia,2,1)
      m_d = obj%q%q_mul_out(ia,3,0)-obj%q%q_mul_out(ia,3,1)
      m = m_s + m_p + m_d
      obj%a%m(ia,1)=m  ! update m
    end do
  case(4)
    do ia=1,obj%a%na
      m_cart = sum(obj%q%q_mul_out(ia,:,1:3),1)
      m_sph = cart2sph(m_cart)
      obj%a%m(ia,:)=m_sph  ! update m
    end do
  end select

end subroutine update_m

end module self_consistent_field_mod
