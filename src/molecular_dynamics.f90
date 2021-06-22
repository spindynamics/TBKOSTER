!
! Copyright (C) 2020
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Pascal Thibaudeau <mailto:pascal.thibaudeau@cea.fr>,
! Ramon Cardias <mailto:ramon.cardias@cea.fr>.
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
!  molecular_dynamics.f90
!  DyNaMol
module molecular_dynamics_mod
    use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
    use atom_mod
    use charge_mod
    use constant_mod, only: gamma_e,k_b,hbar
    use element_mod
    use forces_mod
    use lattice_mod
    use math_mod, only: rad2deg, cross_product, nm2rho, normalize, &
    rho2nm, sph2cart, two_pi
    use precision_mod, only: rp
    use self_consistent_field_mod
    use string_mod, only: dynamol_flush, int2str, log2str, lower, real2str, sl
    use units_mod

    implicit none
    private

    !> Derived type properties for i/o methods
    character(len=sl),dimension(*),parameter :: property_list = &
    [character(len=sl) :: &
    'ensemble', &
    't_i', &
    't_f', &
    'dt', &
    'gamma', &
    'temperature', &
    'verbose' &
    ]

    type,public :: molecular_dynamics
        private
        !> Forces
        class(forces),pointer :: f

        !> ensemble ; options:
        !>  'NVE' (default) : Microcanonical aka NVE
        !>  'NVT' : aka Canonical
        !>  'muNV' : aka Grand Canonical
        !>  'NPH' : aka isoenthalpic-isobaric
        !>  'NPT' : aka isothermal-isobaric

        character(len=4) :: ensemble

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
        !> Thermostated temperature in K (default: 0.0)
        real(rp) :: temperature
        !> gamma damping constant for Langevin Dynamics (default: 0.0)
        !> unit of eV fs / nm^2
        real(rp) :: gamma
        !> @}

        !> @defgroup Molecular_dynamics Molecular- dynamics-related variables
        !> @{
        !> computed Kinetic temperature
        real(rp) :: temp
        !> computed stress tensor
        real(rp), dimension(6) :: stress
        !> @}

        !> Output verbosity flag (default: false)
        logical :: verbose

    contains
        ! Destructor
        final :: destructor
        ! Procedures
        procedure :: initialize
        procedure :: integrate
        procedure :: integrate_nve
        procedure :: integrate_nvt
        procedure :: integrate_muvt
        procedure :: integrate_nph
        procedure :: integrate_npt
        procedure :: integrate_cg
        procedure :: read_txt
        procedure :: write_txt
        procedure :: write_txt_formatted
    end type molecular_dynamics

    ! Constructor
    interface molecular_dynamics
        procedure :: constructor
    end interface molecular_dynamics

contains
    function constructor(f) result(obj)
        class(forces),target,intent(in) :: f
        type(molecular_dynamics) :: obj
  
        obj%f => f
  
        call obj%initialize()
    end function constructor 
    
    subroutine destructor(obj)
		type(molecular_dynamics) :: obj
    end subroutine destructor

    !> Check the validity of the ensemble
    subroutine check_ensemble(ensemble)
        character(len=*),intent(in) :: ensemble

        if(ensemble /= 'cg' &
            .and. ensemble /= 'nve' &
            .and. ensemble /= 'nvt' &
            .and. ensemble /= 'muvt' &
            .and. ensemble /= 'nph' &
            .and. ensemble /= 'npt') then
            write(error_unit,*) 'molecular_dynamics%check_ensemble(): &
            & molecular_dynamics%ensemble must be one of: ''nve'', ''nvt'', &
            ''muvt'', ''nph'', ''npt'', ''cg'''
            error stop
        end if
    end subroutine check_ensemble

    subroutine initialize(obj)
        class(molecular_dynamics),intent(inout) :: obj
    end subroutine initialize

    subroutine initialize_engine(ensemble,t_i,t_f,dt,gamma,temperature)
        character(len=4),intent(out) :: ensemble
        real(rp),intent(out) :: t_i, t_f, dt, gamma, temperature
     
        ensemble = trim('nve')
        t_i = 0.0_rp
        t_f = 0.0_rp
        dt = 0.0_rp
        temperature = 0.0_rp
        gamma = 0.0_rp
    end subroutine initialize_engine

    subroutine integrate(obj,unit)
        class(molecular_dynamics),intent(inout) :: obj
        integer,intent(in),optional :: unit
        integer                     :: unit_rt
    
        if(present(unit)) then
          unit_rt = unit
        else
          unit_rt = output_unit
        end if
    
        ! Integrate the EOM
        select case(obj%ensemble)
        case('nve')
          call obj%integrate_nve(unit_rt)
        case('nvt')
          call obj%integrate_nvt(unit_rt)
        case('muvt')
          call obj%integrate_muvt(unit_rt)
        case('nph')
          call obj%integrate_nph(unit_rt)
        case('npt')
          call obj%integrate_npt(unit_rt)
        case('cg')
          call obj%integrate_cg(unit_rt)
        end select
    end subroutine integrate

    !> Integration with Velocity Verlet scheme to produce NVE ensemble
    subroutine integrate_nve(obj,unit)
        class(molecular_dynamics),intent(inout) :: obj
        integer,intent(in) :: unit
        
        integer :: ia,ia1,it ! atom number and iteration step
        real(rp) :: t,temp ! time and kinetic temperature
        integer :: na ! number of atoms
        real(rp),dimension(:),allocatable :: mass ! atomic mass
        real(rp) :: dt ! time increment
        real(rp), dimension(3) :: old_forces ! forces before update 

        na = obj%f%a_tb%na
        dt = obj%dt ! already converted in hau unit

        ! Initialize time counter and time
        it = 1
        t = obj%t_i

        ! do a SCF calculation do get the forces for current positions
        do ia1=1,na
            obj%f%f_at(ia1,:)=0.0_rp
        end do
        call obj%f%scf%run(unit)
        call obj%f%calculate_forces()

        ! get the atomic masses once and put them in a corrresponding array
        allocate(mass(na))
        do ia=1,na
            mass(ia) = obj%f%e_tb%mass(obj%f%a_tb%ia2ie(ia)) * &
                   obj%f%u%convert_mass('to','hau')
        end do

        ! Time loop
        do while(t < obj%t_f)
            write(unit,'(a)') ''
            write(unit,'(a)') '===== molecular_dynamics%integrate_nve(): iteration ' &
             // int2str(it) // ' ====='

            if(obj%verbose) then
                write(unit,'(a,F15.8)') ' time = ',t * obj%f%u%convert_time('from','hau')
                call obj%f%scf%a%write_txt_formatted(file='',property= &
                 [character(len=sl) :: 'r','p'],tag=.false.,unit=unit,access='append')
            end if

            call obj%f%scf%a%write_txt_formatted(file='out_md.txt',property= &
                 [character(len=sl) :: 'r','p'],tag=.false.,access='append')

            do ia=1,na
                
                old_forces(:) = obj%f%f_at(ia,:)

                ! update the positions
                obj%f%a_tb%r(ia,:)=obj%f%a_tb%r(ia,:)+&
                                    (obj%f%a_tb%p(ia,:)*dt/mass(ia))+&
                                    (0.5_rp*old_forces(:)*dt*dt/mass(ia))
                
                ! do a SCF calculation do get the new forces
                call obj%f%scf%run(unit)
                call obj%f%calculate_forces()

                ! update the impulsions
                obj%f%a_tb%p(ia,:)=obj%f%a_tb%p(ia,:)+&
                                   0.5_rp*(old_forces(:)+obj%f%f_at(ia,:))*dt
            end do

            ! Update iteration step and time counter
            it = it + 1
            t = t + obj%dt
        end do
        deallocate(mass)    
    end subroutine integrate_nve

    !> Integration with the Langevin Dynamcis to produce an NVT ensemble
    subroutine integrate_nvt(obj,unit)
        class(molecular_dynamics),intent(inout) :: obj
        integer,intent(in) :: unit

        integer :: ia,ia1,it ! atom number and iteration step
        real(rp) :: t, ktemp ! time and kinetic temperature
        integer :: na ! number of atoms
        real(rp), dimension(:), allocatable :: mass ! atomic masses
        real(rp) :: dt ! time increment
        real(rp) :: gamma, temperature ! damping constant and temperature
        real(rp) :: b, D ! scaling factor and thermal noise amplitude
        real(rp), dimension(3) :: old_forces, old_positions ! quantities before update
        real(rp), dimension(3) :: force_rand ! random thermal force

        ! implementation follows from : 
        ! A simple and effective Verlet-type algorithm for simulating Langevin dynamics 
        ! https://www.tandfonline.com/doi/abs/10.1080/00268976.2012.760055

        na = obj%f%a_tb%na
        dt = obj%dt ! already converted in hau unit
        gamma = obj%gamma ! already converted in hau unit
        temperature = obj%temperature
        
        ! amplitude of the thermal noise
        D = sqrt(2.0_rp*k_b*obj%f%u%convert_energy('to','hau')*temperature*gamma/dt) 

        ! Initialize time counter and time
        it = 1
        t = obj%t_i

        ! do a SCF calculation do get the forces for current positions
        call obj%f%scf%run(unit)
        call obj%f%calculate_forces()

        ! get the atomic masses once and put them in a corrresponding array
        allocate(mass(na))
        do ia=1,na
            mass(ia) = obj%f%e_tb%mass(obj%f%a_tb%ia2ie(ia)) * &
                   obj%f%u%convert_mass('to','hau')
        end do 

        ! Time loop
        do while(t < obj%t_f)
            write(unit,'(a)') ''
            write(unit,'(a)') '===== molecular_dynamics%integrate_nvt(): iteration ' &
             // int2str(it) // ' ====='

            if(obj%verbose) then
                write(unit,'(a,F15.8)') ' time = ',t * obj%f%u%convert_time('from','hau')
                call obj%f%scf%a%write_txt_formatted(file='',property= &
                 [character(len=sl) :: 'r','p'],tag=.false.,unit=unit,access='append')
            end if

            call obj%f%scf%a%write_txt_formatted(file='out_md.txt',property= &
                 [character(len=sl) :: 'r','p'],tag=.false.,access='append')

            call obj%f%print_forces(unit=unit)        

            do ia=1,na
                old_forces(:)    = obj%f%f_at(ia,:)
                old_positions(:) = obj%f%a_tb%r(ia,:)
                !print *, "r(",ia,")=",old_positions(:) 
                !print *, "p(",ia,")=",obj%f%a_tb%p(ia,:)
                !print *, "f_old(",ia,")=",old_forces(:)
        
                b = 1.0_rp/(1.0_rp+(0.5*gamma*dt/mass(ia)))

                ! get a random thermal field
                force_rand(1) = D*(-1.0_rp+2.0_rp*rand())
                force_rand(2) = D*(-1.0_rp+2.0_rp*rand())
                force_rand(3) = D*(-1.0_rp+2.0_rp*rand())

                ! update the positions
                obj%f%a_tb%r(ia,:)=old_positions(:)+&
                                   (b*obj%f%a_tb%p(ia,:)*dt/mass(ia))+&
                                   (b*0.5_rp*old_forces(:)*dt*dt/mass(ia))+&
                                   0.5_rp*b*dt*dt*force_rand(:)/mass(ia)
                
                ! do a SCF calculation do get the new forces
                do ia1=1,na
                    obj%f%f_at(ia1,:)=0.0_rp
                end do
                call obj%f%a_tb%calculate_neighbours(obj%f%e_tb%r_c_max)
                call obj%f%scf%q%calculate_charge_in()
                call obj%f%scf%h%calculate_h_r()
                call obj%f%scf%h%calculate_s_r()
                call obj%f%scf%run(unit)
                call obj%f%calculate_forces()
                !print *, "f_new(",ia,")=",obj%f%f_at(ia,:)

                ! update the impulsions
                obj%f%a_tb%p(ia,:)=obj%f%a_tb%p(ia,:)+&
                                   0.5_rp*(obj%f%f_at(ia,:)+old_forces(:))*dt-&
                                   gamma*(obj%f%a_tb%r(ia,:)-old_positions(:))+&
                                   force_rand(:)*dt
            end do

            ! Update iteration step and time counter
            it = it + 1
            t = t + obj%dt
            !print *,"====="
        end do
        deallocate(mass)    
    end subroutine integrate_nvt

    subroutine integrate_muvt(obj,unit)
        class(molecular_dynamics),intent(inout) :: obj
        integer,intent(in) :: unit
    end subroutine integrate_muvt

    subroutine integrate_nph(obj,unit)
        class(molecular_dynamics),intent(inout) :: obj
        integer,intent(in) :: unit
    end subroutine integrate_nph

    subroutine integrate_npt(obj,unit)
        class(molecular_dynamics),intent(inout) :: obj
        integer,intent(in) :: unit
    end subroutine integrate_npt

    !> Integration with conjugate gradient to minimize forces
    subroutine integrate_cg(obj,unit)
        class(molecular_dynamics),intent(inout) :: obj
        integer,intent(in) :: unit
    end subroutine integrate_cg

    !> Read object in text format from file (default: 'in_md.txt')
    subroutine read_txt(obj,file)
        class(molecular_dynamics),intent(inout) :: obj
        character(len=*),intent(in),optional :: file
        character(len=:),allocatable :: file_rt
        integer :: iostatus
        logical :: isopen
        ! Namelist variables
        character(len=4) :: ensemble
        real(rp) :: t_i, t_f, dt, gamma, temperature
        logical :: verbose
        ! Namelist
        namelist /md/ ensemble, t_i, t_f, dt, gamma, temperature, verbose
        ! Local variables
        real(rp) :: coeff

        if (present(file)) then
            file_rt = trim(file)
        else
            file_rt = 'in_md.txt'
        end if

        inquire(unit=10, opened=isopen)
        if (isopen) then
            write(error_unit,'(a)') 'md%read_txt() : Unit 10 is already open'
            error stop
        else
            open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
        end if
        if (iostatus /= 0) then
            write(error_unit,*) 'md%read_txt(): file ', file_rt, ' not found'
            error stop
        end if

        call initialize_engine(ensemble, t_i, t_f, dt, gamma, temperature)
        verbose = .false.

        read(10,nml=md)
        ensemble = trim(lower(ensemble))
        call check_ensemble(ensemble)

        obj%ensemble = ensemble
        obj%verbose = verbose
        obj%t_i = t_i * obj%f%u%convert_time('to','hau')
        obj%t_f = t_f * obj%f%u%convert_time('to','hau')
        obj%dt = dt * obj%f%u%convert_time('to','hau')
        obj%temperature = temperature
        coeff = obj%f%u%convert_energy('to','hau')* &
                obj%f%u%convert_time('to','hau')/ &
                obj%f%u%convert_length('to','hau')/ &
                obj%f%u%convert_length('to','hau')
        obj%gamma = gamma * coeff

        close(unit=10)
    end subroutine read_txt

    subroutine write_txt(obj,file,unit)
        class(molecular_dynamics) :: obj
        character(len=*),intent(in),optional :: file
        character(len=:),allocatable         :: file_rt
        integer,intent(in),optional :: unit
        integer                     :: unit_rt
        ! Namelist variables
        character(len=4) :: ensemble
        real(rp) :: t_i, t_f, dt, gamma, temperature
        logical  :: verbose
        ! Namelist
        namelist /md/ ensemble, t_i, t_f, dt, gamma, temperature, verbose
        ! Local variables
        real(rp) :: coeff

        if (present(file)) then
            file_rt = file
        else
            file_rt = 'out_md.txt'
        end if
        if (present(unit)) then
            unit_rt = unit
        else
            unit_rt = 10
        end if

        if (.not. present(unit)) then
            open(unit=unit_rt,file=file_rt,action='write')
        end if

        ensemble = obj%ensemble
        t_i = obj%t_i * obj%f%u%convert_time('from','hau')
        t_f = obj%t_f * obj%f%u%convert_time('from','hau')
        dt  = obj%dt * obj%f%u%convert_time('from','hau')
        temperature = obj%temperature
        coeff = obj%f%u%convert_energy('from','hau')* &
                obj%f%u%convert_time('from','hau')/ &
                obj%f%u%convert_length('from','hau')/ &
                obj%f%u%convert_length('from','hau')
        gamma = obj%gamma * coeff
        verbose = obj%verbose

        write(unit_rt,nml=md)

        if(.not. present(unit)) then
            call dynamol_flush(unit_rt)
            close(unit_rt)
        else
            call dynamol_flush(unit)
        end if
    
    end subroutine write_txt

    subroutine write_txt_formatted(obj,file,property,tag,unit)
        class(molecular_dynamics) :: obj
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
        real(rp) :: coeff

        if (present(file)) then
            file_rt = file
        else
            file_rt = 'out_md.txt'
        end if
        if (present(property)) then
            property_rt = property
        else
            property_rt = property_list
        end if
        if (present(tag)) then
            tag_rt = tag
        else
            tag_rt = .true.
        end if
        if (present(unit)) then
            unit_rt = unit
        else
            unit_rt = 10
        end if

        if(.not. present(unit)) then
            open(unit=unit_rt,file=file_rt,action='write')
        end if
        if(tag_rt) then
            write(unit_rt,'(a)') '&md'
        end if

        do ip=1,size(property_rt)
            select case(lower(trim(property_rt(ip))))
            case('ensemble')
            write(unit_rt,'(a)') ' ensemble = ''' // trim(obj%ensemble) // ''''
            case('t_i')
            write(unit_rt,'(a)') ' t_i = ' // real2str(obj%t_i &
            * obj%f%u%convert_time('from','hau'))
            case('t_f')
            write(unit_rt,'(a)') ' t_f = ' // real2str(obj%t_f &
            * obj%f%u%convert_time('from','hau'))
            case('dt')
            write(unit_rt,'(a)') ' dt = ' // real2str(obj%dt &
            * obj%f%u%convert_time('from','hau'))
            case('gamma')
                coeff = obj%f%u%convert_energy('from','hau')* &
                obj%f%u%convert_time('from','hau')/ &
                obj%f%u%convert_length('from','hau')/ &
                obj%f%u%convert_length('from','hau')
            write(unit_rt,'(a)') ' gamma = ' // real2str(obj%gamma * coeff)
            case('temperature')
            write(unit_rt,'(a)') ' temperature = ' // real2str(obj%temperature)
            case('verbose')
            write(unit_rt,'(a)') ' verbose = ' // log2str(obj%verbose)
            end select
        end do

        if(tag_rt) then
            write(unit_rt,'(a)') ' /'
        end if
        if(.not. present(unit)) then
            call dynamol_flush(unit_rt)
            close(unit_rt)
        else
            call dynamol_flush(unit)
        end if
    
    end subroutine write_txt_formatted

end module molecular_dynamics_mod