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
!  spin_dynamics.f90
!  TBKOSTER
module spin_dynamics_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_mod
  use charge_mod
  use constant_mod, only: gamma_e,k_b,hbar,e_ha
  use element_mod
  use lattice_mod
  use math_mod, only: rad2deg, cart2sph, cross_product, nm2rho, normalize, &
  rho2nm, sph2cart, two_pi
  use precision_mod, only: rp
  use self_consistent_field_mod
  use string_mod, only: TBKOSTER_flush, int2str, log2str, lower, real2str, sl
  use units_mod
#if defined(INTEL_COMPILER)
  use ifport
#endif
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'integrator', &
   'engine', &
   't_i', &
   't_f', &
   'dt', &
   'fixed_time_step', &
   'quality_factor', &
   'alpha', &
   'temp', &
   'verbose', &
   'compute_effective_field_every', &
   'compute_effective_field_atomic_loop' &
   ]

  type,public :: spin_dynamics
    private
    !> Units
    class(units),pointer :: u
    !> Elements
    class(element),pointer :: e
    !> Atoms
    class(atom),pointer :: a
    !> Charges
    class(charge),pointer :: q
    !> Self-consistent field
    class(self_consistent_field),pointer :: scf

    !> Integrator ; options:
    !>  'euler' (default) : Euler explicit scheme
    !>  'st_1' : Suzuki-Trotter geometric decomposition of order 1
    character(len=5) :: integrator

    !> Engine ; options:
    !> 'approx' (default) : Approximate motion on a sphere for small variations
    !> 'exact' : Exact rotation on a sphere
    character(len=6) :: engine

    !> @defgroup Time Time-related variables
    !> @{

    !> Number of time points
    integer :: nt
    !> Initial time in fs (default: 0.0)
    real(rp) :: t_i
    !> Final time in fs (default: 0.0)
    real(rp) :: t_f
    !> Time step in fs (default: 0.0)
    real(rp) :: dt
    !> Fixed time step flag (default: true)
    logical :: fixed_time_step
    !> Quality factor for time step determination (default: 0.0)
    real(rp) :: quality_factor
    !> @}

    !> Gilbert damping constant (default: 0.0, typical values: 0.001-1.0)
    real(rp) :: alpha

    !> @defgroup Langevin_dynamics Langevin dynamics-related variables
    !> @{

    !> Spin temperature (default: 0.0)
    real(rp) :: temp
    !> Stochastic magnetic induction
    real(rp),dimension(:,:),allocatable :: b_stochastic
    !> @}

    !> Output verbosity flag (default: false)
    logical :: verbose

    !> @defgroup compute_effective_field Effective field computation-related variables
    !> @{
    !> Compute effective field every timesteps (default:1)
    integer :: compute_effective_field_every

    !> true if the effective field is computed inside the loop on each atoms
    !> false if the effective field is computed outside the loop on each atoms
    logical :: compute_effective_field_atomic_loop
    !> @}

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: advance_st
    procedure :: advance_st_approx
    procedure :: advance_st_exact
    procedure :: calculate_spin_temperature
    procedure :: estimate_timestep
    procedure :: initialize
    procedure :: integrate
    procedure :: integrate_euler
    procedure :: integrate_st_1
    procedure :: read_txt
    procedure :: set_b_stochastic
    procedure :: update_dt
    procedure :: write_txt
    procedure :: write_txt_formatted
    procedure :: write_vectors
  end type spin_dynamics

  interface spin_dynamics
    procedure :: constructor
  end interface spin_dynamics

contains
  function constructor(scf) result(obj)
    class(self_consistent_field),target,intent(in) :: scf
		type(spin_dynamics) :: obj

    obj%u => scf%u
    obj%e => scf%a%e
    obj%a => scf%a
    obj%q => scf%q
    obj%scf => scf

    call obj%initialize()
  end function constructor

  subroutine destructor(obj)
		type(spin_dynamics) :: obj

    if(allocated(obj%b_stochastic)) deallocate(obj%b_stochastic)
  end subroutine destructor

  !> Advance a single magnetization using the Suzuki-Trotter geometric
  !> decomposition on a unitary sphere,
  !> depending if the small angle approximation is done ='approx' or not ='exact'
  subroutine advance_st(obj,unit,ia,dt)
    class(spin_dynamics),intent(inout) :: obj
    integer,intent(in) :: unit
    integer,intent(in) :: ia
    real(rp),intent(in) :: dt

    select case(obj%engine)
    case('approx')
      call advance_st_approx(obj,unit,ia,dt)
    case('exact')
      call advance_st_exact(obj,unit,ia,dt)
    end select
  end subroutine advance_st

  !> Advance a single magnetization using the Suzuki-Trotter geometric
  !> approximate decomposition on a unitary sphere
  subroutine advance_st_approx(obj,unit,ia,dt)
		class(spin_dynamics),intent(inout) :: obj
    integer,intent(in) :: unit
    integer,intent(in) :: ia
    real(rp),intent(in) :: dt
    real(rp),dimension(3) :: b_o, b, m_o, m_n, v_sph
    real(rp) :: b_2, dt_2, norm_m, t_norm

    !integer :: ie,io,l
    !complex(rp),dimension(2,2) :: rho
    !real(rp) :: n
    !real(rp),dimension(3) :: m_cart

    ! Calculate the new magnetization m_n from the old magnetization m_o
    v_sph(:) = obj%a%m(ia,:)
    m_o = sph2cart(v_sph)

    b_o = -obj%a%b_pen(ia,:)*v_sph(1) ! twice the renormalization of b_pen by m (see paper)
    b = (b_o-(obj%alpha)*(cross_product(b_o,m_o)/v_sph(1)))/(1.0_rp+(obj%alpha)*(obj%alpha))

    b_2 = dot_product(b,b)
    dt_2 = dt*dt
    !write(output_unit,"(a,g10.4)") "bt_2*dt_2=", sqrt(b_2*dt_2)

    !write(unit,'(a)') 'spin_dynamics%advance_st(): m_o = ' &
    ! // real2str(obj%a%m(ia,1)) &
    ! // ', ' // real2str(obj%a%m(ia,2)*rad2deg) &
    ! // ', ' // real2str(obj%a%m(ia,3)*rad2deg)

    m_n = (m_o + dt*cross_product(b,m_o) + dt_2/2*(dot_product(b,m_o)*b &
       - b_2*m_o/2))/(1+b_2*dt_2/4)

    ! impose that m_o and m_n keep the same norm after the rotation
    t_norm = norm2(m_n)
    m_n(:) = m_n(:)*v_sph(1)/t_norm
    ! Update the magnetization and m_pen from cartesian to spherical coordinates
    obj%a%m(ia,:) = cart2sph(m_n)
    obj%a%m_pen(ia,:) = cart2sph(m_n)
    !write(unit,'(a)') 'spin_dynamics%advance_st(): m_n = ' &
    ! // real2str(obj%a%m(ia,1)) &
    ! // ', ' // real2str(obj%a%m(ia,2)*rad2deg) &
    ! // ', ' // real2str(obj%a%m(ia,3)*rad2deg)

    ! Update the charges
    !call obj%q%calculate_charge_in((/ia/))
    !obj%q%q_mul_in(ia,1:2,1:3) = 0.0_rp
    ! obj%q%q_mul_in(ia,3,1:3) = m_n
    ! ! Update the net charges
    ! ie = obj%a%ia2ie(ia)
    ! do io=1,obj%e%no(ie)
    !   l = obj%e%o2l(obj%e%o(ie,io))
    !   if(l==3) then
    !     rho(1,1) = obj%q%rho_net_in(ia,io,io,1)
    !     rho(2,2) = obj%q%rho_net_in(ia,io,io,2)
    !     rho(1,2) = cmplx(0.0_rp,0.0_rp,kind=rp)
    !     rho(2,1) = cmplx(0.0_rp,0.0_rp,kind=rp)
    !     call rho2nm(rho,n,m_cart)
    !     m_cart = sph2cart((/obj%a%m(ia,1)/5.0_rp,obj%a%m(ia,2),obj%a%m(ia,3)/))
    !     call nm2rho(n,m_cart,rho)
    !     obj%q%rho_net_in(ia,io,io,1) = rho(1,1)
    !     obj%q%rho_net_in(ia,io,io,2) = rho(2,2)
    !     obj%q%rho_net_in(ia,io,io,3) = rho(1,2)
    !     obj%q%rho_net_in(ia,io,io,4) = rho(2,1)
    !   end if
    ! end do
    ! call obj%q%write_mulliken_charge_analysis(unit,'in')
  end subroutine advance_st_approx

  !> Advance a single magnetization using the Suzuki-Trotter geometric
  !> exact decomposition on a unitary sphere
  subroutine advance_st_exact(obj,unit,ia,dt)
    class(spin_dynamics),intent(inout) :: obj
    integer,intent(in) :: unit
    integer,intent(in) :: ia
    real(rp),intent(in) :: dt
    real(rp),dimension(3) :: b_o, b, m_o, m_n, DS, v_sph
    real(rp) :: omega2, energy, xi, chi
    real(rp) :: norm_of_omega, norm_m, t_norm

    ! Calculate the new magnetization m_n from the old magnetization m_o
    v_sph(:) = obj%a%m(ia,:)
    m_o = sph2cart(v_sph)
        
    b_o = -obj%a%b_pen(ia,:)*v_sph(1) ! twice the renormalization of b_pen by m (see paper)
    b = (b_o-((obj%alpha)*(cross_product(b_o,m_o)/v_sph(1))))/(1.0_rp+(obj%alpha)*(obj%alpha))

    omega2 = dot_product(b,b)
    norm_of_omega = sqrt(omega2)

    energy = dot_product(b,m_o)
    DS = cross_product(b,m_o)
    xi = norm_of_omega*dt

    chi=energy/norm_of_omega

    !
    ! Thibaudeau, P. and Beaujouan, D. : Physica A: Statistical Mechanics and
    ! its Applications : doi:10.1016/j.physa.2011.11.030 : exact solution of a
    ! time evolution of spin in a constant field and damping
    ! The norm is automatically conserved upon rotations

    m_n(:)=(cos(xi)*m_o(:))+(sin(xi)*DS(:)/norm_of_omega)+&
    (chi*(1.0_rp-cos(xi))*b(:)/norm_of_omega)

    ! Update the magnetization and m_pen from cartesian to spherical coordinates
    obj%a%m(ia,:) = cart2sph(m_n)
    obj%a%m_pen(ia,:) = cart2sph(m_n)

  end subroutine advance_st_exact

  !> Calculate the spin temperature according to
  !> https://doi.org/10.1103/PhysRevE.61.3579
  function calculate_spin_temperature(obj) result(temp)
		class(spin_dynamics),intent(in) :: obj
    real(rp) :: energy, dots, temp, factor
    real(rp), dimension(3) :: s,torque, v_sph, b_pen
    integer :: ia

    energy = 0.0_rp
    dots = 0.0_rp

    do ia=1,obj%a%na
      v_sph = obj%a%m(ia,:)
      s = normalize(sph2cart(v_sph))
      b_pen = obj%a%b_pen(ia,:)
      energy = energy + dot_product(b_pen,s)
      torque(:) = cross_product(b_pen,s)
      dots = dots + dot_product(torque,torque)
    end do

    factor = abs(obj%t_f-obj%t_i)/obj%dt

    temp=dots/2.0_rp/energy/(k_b*obj%u%convert_energy('to','hau'))/factor

  end function calculate_spin_temperature

  !> Check the validity of the engine
  subroutine check_engine(engine)
    character(len=*),intent(in) :: engine

    if(engine /= 'approx' &
     .and. engine /= 'exact') then
      write(error_unit,*) 'spin_dynamics%check_engine(): &
       &spin_dynamics%engine must be one of: ''approx'', ''exact'''
      error stop
    end if
  end subroutine check_engine

  !> Check the validity of the integrator
  subroutine check_integrator(integrator)
    character(len=*),intent(in) :: integrator

    if(integrator /= 'euler' &
     .and. integrator /= 'st_1') then
      write(error_unit,*) 'spin_dynamics%check_integrator(): &
       &spin_dynamics%integrator must be one of: ''euler'', ''st_1'''
      error stop
    end if
  end subroutine check_integrator
  
  !> Estimate a better timestep according to the amplitude of effective pulsation
  subroutine estimate_timestep(obj,omega)
    class(spin_dynamics),intent(inout) :: obj
    real(rp), dimension(3), intent(in) :: omega
    real(rp) :: MaxNorm;
    MaxNorm = (obj%quality_factor)/norm2(omega)
    obj%dt=merge(obj%dt,MaxNorm,obj%dt<MaxNorm)
  end subroutine estimate_timestep

  subroutine initialize(obj)
    class(spin_dynamics),intent(inout) :: obj

    if(allocated(obj%b_stochastic)) deallocate(obj%b_stochastic)
    allocate(obj%b_stochastic(obj%a%na,3))

    obj%b_stochastic = 0.0_rp
  end subroutine initialize

  subroutine initialize_damping(alpha)
    real(rp),intent(out) :: alpha

    alpha = 0.0_rp
  end subroutine initialize_damping

  subroutine initialize_temperature(temp)
    real(rp),intent(out) :: temp

    temp = 0.0_rp
  end subroutine initialize_temperature

  subroutine initialize_time(integrator,engine,t_i,t_f,dt,fixed_time_step, &
   quality_factor, &
   compute_effective_field_every,compute_effective_field_atomic_loop)
    character(len=5),intent(out) :: integrator
    character(len=6),intent(out) :: engine
    real(rp),intent(out) :: t_i, t_f, dt
    logical,intent(out) :: fixed_time_step
    real(rp),intent(out) :: quality_factor
    integer,intent(out) :: compute_effective_field_every
    logical,intent(out) :: compute_effective_field_atomic_loop

    integrator = 'euler'
    engine = 'approx'
    t_i = 0.0_rp
    t_f = 0.0_rp
    dt = 0.0_rp
    fixed_time_step = .true.
    quality_factor = 0.0_rp
    compute_effective_field_every = 1
    compute_effective_field_atomic_loop = .true.
  end subroutine initialize_time

  !> Select the integration routine based on the integrator
  subroutine integrate(obj,unit)
    class(spin_dynamics),intent(inout) :: obj
    integer,intent(in),optional :: unit
    integer                     :: unit_rt

    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = output_unit
    end if

    ! Create the directory "sd" where the magnetizations will be written
    call execute_command_line('mkdir -p sd')

    ! Integrate the EOM
    select case(obj%integrator)
    case('euler')
      call obj%integrate_euler(unit_rt)
    case('st_1')
      call obj%integrate_st_1(unit_rt)
    end select
  end subroutine integrate

  !> Integrate the EOM using explicit Euler scheme
  subroutine integrate_euler(obj,unit)
    class(spin_dynamics),intent(inout) :: obj
    integer,intent(in) :: unit
    integer :: ia,it
    real(rp) :: t,temp

    ! Write the initial magnetizations of the atoms
    call obj%a%write_txt_formatted(file='sd/out_atom_tb_0.txt')
    !call obj%a%write_txt_formatted(file='sd/out_atom_tb_0.txt', &
    ! property=[character(len=sl) :: 'm'])
    if(obj%verbose) then
      call obj%a%write_txt_formatted(unit=unit,property= &
       [character(len=sl) :: 'm'],tag=.false.)
    end if

    ! Initialize time counter and time
    it = 1
    t = obj%t_i

    ! Time loop
    do while(t < obj%t_f)
      write(unit,'(a)') ''
      write(unit,'(a)') '===== spin_dynamics%integrate_euler(): iteration ' &
       // int2str(it) // ' ====='

      ! Update time step if necessary
      if(.not. obj%fixed_time_step) then
        call obj%update_dt()
      end if

      if(obj%verbose) then
        call obj%write_vectors(t)
      end if

      ! compute the stochastic field for all atoms from 1 by 1 step
      call obj%set_b_stochastic(1,1)

      ! update b_pen with the stochastic field
      obj%a%b_pen = obj%a%b_pen+obj%b_stochastic

      ! compute the effective field once according this flag
      if(.not. obj%compute_effective_field_atomic_loop) then
        call obj%scf%run(unit)
      end if

      ! Advance the na magnetizations by dt
      do ia=1,obj%a%na
        write(unit,'(a)') ''
        write(unit,'(a)') '===== spin_dynamics%integrate_euler(): atom ' &
         // int2str(ia) // ' ====='
        ! run a scf calculation to update the effective fields every known steps
        if (mod(it,obj%compute_effective_field_every)==0 .and. obj%compute_effective_field_atomic_loop) then
          call obj%scf%run(unit)
          if(.not. obj%scf%converged) then
            write(unit,'(a)') 'Convergence is not reached, SD stopped'
            exit
          end if
        end if
        call obj%advance_st(unit,ia,obj%dt)
      end do

      ! Write the instantaneous magnetizations of the atoms
      call obj%a%write_txt_formatted(file='sd/out_atom_tb_' // int2str(it) &
       // '.txt')
      !call obj%a%write_txt_formatted(file='sd/out_atom_tb_' // int2str(it) &
      ! // '.txt', property=[character(len=sl) :: 'm'])
      if(obj%verbose) then
        call obj%a%write_txt_formatted(unit=unit,property= &
         [character(len=sl) :: 'm'],tag=.false.)
         ! compute the spin themperature
        temp = calculate_spin_temperature(obj)
        write(unit,'(a)') ' timestep = '//real2str(t*obj%u%convert_time('from','hau'))&
        //' ['//trim(obj%u%time)//']'
        write(unit,'(a)') ' spin temperature = ' //real2str(temp)//' [K]'
      end if

      ! Write output to the lammps format
      call obj%a%write_lammps(t,file='out_lammps.lammpstrj')
      ! Update time counter and time
      it = it + 1
      t = t + obj%dt
    end do
  end subroutine integrate_euler

  !> Integrate the EOM using Suzuki-Trotter geometric decomposition scheme
  !> of first-order
  subroutine integrate_st_1(obj,unit)
		class(spin_dynamics),intent(inout) :: obj
    integer,intent(in) :: unit
    integer :: ia,it
    real(rp) :: t,temp

    ! Write the initial magnetizations of the atoms in a numbered file
    call obj%a%write_txt_formatted(file='sd/out_atom_tb_0.txt')
    !call obj%a%write_txt_formatted(file='sd/out_atom_tb_0.txt', &
    ! property=[character(len=sl) :: 'm'])

    ! Initialize time counter and time
    it = 1
    t = obj%t_i

    ! Time loop
    do while(t < obj%t_f)
      write(unit,'(a)') ''
      write(unit,'(a)') '===== spin_dynamics%integrate_st_1(): iteration ' &
       // int2str(it) // ' ====='
      ! Update time step
      if(.not. obj%fixed_time_step) then
        call obj%update_dt()
      end if

      if(obj%verbose) then
        call obj%write_vectors(t)
      end if

      ! compute the stochastic field for all atoms from 1 by 1 step
      call obj%set_b_stochastic(1,1)

      ! update b_pen with the stochastic field
      obj%a%b_pen = obj%a%b_pen+obj%b_stochastic

      ! compute the effective field once according this flag

      ! Advance the first na-1 magnetizations in ascending order by dt/2
      ! Update the effective field by a scf%run
      if(.not. obj%compute_effective_field_atomic_loop) then
        call obj%scf%run(unit)
      end if

      do ia=1,obj%a%na-1
        ! run a scf calculation to update the effective fields every known steps
        if (mod(it,obj%compute_effective_field_every)==0 .and. obj%compute_effective_field_atomic_loop) then
          call obj%scf%run(unit)
          if(.not. obj%scf%converged) then
            write(unit,'(a)') 'Convergence is not reached, SD stopped'
            exit
          end if
        end if
        call obj%advance_st(unit,ia,obj%dt/2)
      end do
      ! Advance the na-th magnetization by dt
      ! Update the effective field by a scf%run
      ! run a scf calculation to update the effective fields every known steps
      if(.not. obj%compute_effective_field_atomic_loop) then
        call obj%scf%run(unit)
      end if
      if (mod(it,obj%compute_effective_field_every)==0 .and. obj%compute_effective_field_atomic_loop) then
        call obj%scf%run(unit)
        if(.not. obj%scf%converged) then
          write(unit,'(a)') 'Convergence is not reached, SD stopped'
          exit
        end if
      end if
      call obj%advance_st(unit,obj%a%na,obj%dt)
      ! Advance the first na-1 magnetizations in descending order by dt/2
      ! Update the effective field by a scf%run
      if(.not. obj%compute_effective_field_atomic_loop) then
        call obj%scf%run(unit)
      end if
      do ia=obj%a%na-1,1,-1
        ! run a scf calculation to update the effective fields every known steps
        if (mod(it,obj%compute_effective_field_every)==0 .and. obj%compute_effective_field_atomic_loop) then
          call obj%scf%run(unit)
          if(.not. obj%scf%converged) then
            write(unit,'(a)') 'Convergence is not reached, SD stopped'
            exit
          end if
        end if
        call obj%advance_st(unit,ia,obj%dt/2)
      end do

      ! Write the instantaneous magnetizations of the atoms in a numbered file
      call obj%a%write_txt_formatted(file='sd/out_atom_tb_' // int2str(it) &
       // '.txt')
      !call obj%a%write_txt_formatted(file='sd/out_atom_' // int2str(it) &
      ! // '.txt', property=[character(len=sl) :: 'm'])

      if(obj%verbose) then
        call obj%a%write_txt_formatted(unit=unit,property= &
         [character(len=sl) :: 'm'],tag=.false.)
         write(unit,'(a)') ' timestep = '//real2str(t*obj%u%convert_time('from','hau'))&
         //' ['//trim(obj%u%time)//']'
         ! compute the spin themperature
        temp = calculate_spin_temperature(obj)
        write(unit,'(a)') ' spin temperature = ' // real2str(temp)
      end if

      ! Write output to the lammps format
      call obj%a%write_lammps(t,file='out_lammps.lammpstrj')
      ! Update time counter and time
      ! Update time counter and time
      it = it + 1
      t = t + obj%dt
    end do
  end subroutine integrate_st_1

  !> Read object in text format from file (default: 'in_sd.txt')
  subroutine read_txt(obj,file)
    class(spin_dynamics),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    character(len=5) :: integrator
    character(len=6) :: engine
    real(rp) :: t_i, t_f, dt
    logical :: fixed_time_step
    real(rp) :: quality_factor
    real(rp) :: alpha
    real(rp) :: temp
    logical :: verbose
    integer :: compute_effective_field_every
    logical :: compute_effective_field_atomic_loop
    ! Namelist
    namelist /sd/ integrator, engine, t_i, t_f, dt, fixed_time_step, &
    quality_factor, alpha, temp, verbose, compute_effective_field_every, &
    compute_effective_field_atomic_loop

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_sd.txt'
    end if

    inquire(unit=10, opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'sd%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'sd%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    call initialize_time(integrator,engine,t_i,t_f,dt,fixed_time_step,&
                         quality_factor,compute_effective_field_every,&
                         compute_effective_field_atomic_loop)
    call initialize_damping(alpha)
    call initialize_temperature(temp)
    verbose = .false.
    read(10,nml=sd)
    integrator = trim(lower(integrator))
    call check_integrator(integrator)
    engine = trim(lower(engine))
    call check_engine(engine)

    obj%integrator = integrator
    obj%engine = engine
    obj%t_i = t_i * obj%u%convert_time('to','hau')
    obj%t_f = t_f * obj%u%convert_time('to','hau')
    obj%dt = dt * obj%u%convert_time('to','hau')
    obj%fixed_time_step = fixed_time_step
    obj%quality_factor = quality_factor
    obj%alpha = alpha
    obj%temp = temp
    obj%verbose = verbose
    obj%compute_effective_field_every = compute_effective_field_every
    obj%compute_effective_field_atomic_loop = compute_effective_field_atomic_loop

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Set the stochastic magnetic field for Langevin dynamics
  subroutine set_b_stochastic(obj,init,step)
    class(spin_dynamics),intent(inout) :: obj
    integer, intent(in) :: init ! number of atom to start
    integer, intent(in) :: step ! stepping on atoms

    integer :: ia ! index on atoms
    integer :: id ! index on space dimensions
    integer :: count ! counter for time
    real(rp) :: coefficient ! noise amplitude

    real(rp) :: k_b_hau, hbar_hau ! constants in internal units

    k_b_hau = k_b*(obj%u%convert_energy('to','hau'))
    hbar_hau = hbar*(obj%u%convert_energy('to','hau'))

    call system_clock(count)
    call srand(count)

    coefficient = sqrt((obj%alpha)*two_pi*k_b_hau*(obj%temp)/hbar_hau/(obj%dt))

    do ia = init,obj%a%na,step
      do id = 1,3
        obj%b_stochastic(ia,id) = coefficient*(-1.0_rp+2.0_rp*rand())
      end do
    end do

  end subroutine set_b_stochastic

  !> Update the time step for an optimal stepping scheme
  subroutine update_dt(obj)
		class(spin_dynamics),intent(in) :: obj

  end subroutine update_dt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_sd.txt')
  subroutine write_txt(obj,file,unit)
    class(spin_dynamics) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=5) :: integrator
    character(len=6) :: engine
    real(rp) :: t_i, t_f, dt
    logical  :: fixed_time_step
    real(rp) :: quality_factor
    real(rp) :: alpha
    real(rp) :: temp
    logical  :: verbose
    integer  :: compute_effective_field_every
    logical  :: compute_effective_field_atomic_loop
    ! Namelist
    namelist /sd/ integrator, engine, t_i, t_f, dt, fixed_time_step, &
    quality_factor, alpha, temp, verbose, compute_effective_field_every, &
    compute_effective_field_atomic_loop

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_sd.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    integrator = obj%integrator
    engine = obj%engine
    t_i = obj%t_i * obj%u%convert_time('from','hau')
    t_f = obj%t_f * obj%u%convert_time('from','hau')
    dt  = obj%dt * obj%u%convert_time('from','hau')
    fixed_time_step = obj%fixed_time_step
    quality_factor  = obj%quality_factor
    alpha = obj%alpha
    temp = obj%temp
    verbose = obj%verbose
    compute_effective_field_every = obj%compute_effective_field_every
    compute_effective_field_atomic_loop = obj%compute_effective_field_atomic_loop

    write(unit_rt,nml=sd)

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
  !> 'out_sd.txt'), if tag (default: .true.) the namelist opening and closing
  !> tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(spin_dynamics) :: obj
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
      file_rt = 'out_sd.txt'
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
      write(unit_rt,'(a)') '&sd'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('integrator')
        write(unit_rt,'(a)') ' integrator = ''' // trim(obj%integrator) // ''''
      case('engine')
        write(unit_rt,'(a)') ' engine = ''' // trim(obj%engine) //''''
      case('t_i')
        write(unit_rt,'(a)') ' t_i = ' // real2str(obj%t_i &
         * obj%u%convert_time('from','hau'))
      case('t_f')
        write(unit_rt,'(a)') ' t_f = ' // real2str(obj%t_f &
         * obj%u%convert_time('from','hau'))
      case('dt')
        write(unit_rt,'(a)') ' dt = ' // real2str(obj%dt &
         * obj%u%convert_time('from','hau'))
      case('fixed_time_step')
        write(unit_rt,'(a)') ' fixed_time_step = ' &
         // log2str(obj%fixed_time_step)
      case('quality_factor')
        write(unit_rt,'(a)') ' quality_factor = ' &
         // real2str(obj%quality_factor)
       case('alpha')
         write(unit_rt,'(a)') ' alpha = ' // real2str(obj%alpha)
       case('temp')
         write(unit_rt,'(a)') ' temp = ' // real2str(obj%temp)
       case('b_stochastic')
         ! TO DO
       case('verbose')
         write(unit_rt,'(a)') ' verbose = ' // log2str(obj%verbose)
       case('compute_effective_field_every')
         write(unit_rt,'(a)') ' compute_effective_field_every = ' &
          // int2str(obj%compute_effective_field_every)
       case('compute_effective_field_atomic_loop')
         write(unit_rt,'(a)') ' compute_effective_field_atomic_loop = ' &
          // log2str(obj%compute_effective_field_atomic_loop)
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

  !> output time step, magnetizations and effective fields on selected unit (default=6)
  subroutine write_vectors(obj,timestep,unit)
    class(spin_dynamics), intent(in) :: obj
    real(rp),intent(in),optional :: timestep 
    integer,intent(in),optional :: unit
    integer :: unit_rt, ia
    real(rp),dimension(:,:),allocatable :: torque, m_cart
    real(rp),dimension(3) :: v_sph,b_cart,b_eff

    if (.not.present(unit)) then
      unit_rt=6
    else
      unit_rt=unit
    end if

    ! compute the torque b x m for all atoms
    allocate(m_cart(obj%a%na,3))
    allocate(torque(obj%a%na,3))
  
    do ia=1,obj%a%na
      v_sph(:) = obj%a%m(ia,:)
      b_cart(:) = -obj%a%b_pen(ia,:)
      m_cart(ia,:) = sph2cart(v_sph)
      b_eff(:) = (b_cart(:)-(obj%alpha)*cross_product(b_cart,m_cart)/v_sph(1))/(1.0_rp+(obj%alpha)*(obj%alpha))
      torque(ia,:) = cross_product(b_eff(:),m_cart(ia,:))
    end do

    ! torque output is hard coded in fs^{-1}
    write(unit_rt,*) timestep*obj%u%convert_time('from','hau'),&
    (m_cart(ia,:),ia=1,obj%a%na),&
    (torque(ia,:)*e_ha/hbar,ia=1,obj%a%na)
    
    deallocate(torque,m_cart)
  end subroutine write_vectors

end module spin_dynamics_mod
