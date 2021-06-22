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
!  forces.f90
!  DyNaMol
module forces_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_tb_mod
  use charge_mod
  use element_tb_mod
  use lattice_mod
  use math_mod, only: fermi_function
  use precision_mod, only: rp
  use string_mod, only: dynamol_flush, fixedreal2str, int2str, lower, real2str, sl, &
                        cmplx2str, log2str
  use self_consistent_field_mod
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'computed' &
   ]

  type,public :: forces
    !> Units
    class(units),pointer :: u
    ! Elements TB
    class(element_tb),pointer :: e_tb
    ! Atoms TB
    class(atom_tb),pointer :: a_tb
    !> Charges
    class(charge),pointer :: q
    !> Self-consistent field
    class(self_consistent_field),pointer :: scf
    !> Forces
    real(rp), dimension(:,:), allocatable :: f_at

    ! logical if forces are computed .true.=computed, .false.=not-computed
    logical :: computed
   
  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: build_forces_rho_per_point
    procedure :: calculate_forces_per_point
    procedure :: calculate_forces
    procedure :: initialize
    procedure :: print_forces
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type forces         

  ! Constructor
  interface forces        
    procedure :: constructor
  end interface forces          

contains
  function constructor(scf) result(obj)
    class(self_consistent_field),target,intent(in) :: scf
    type(forces) :: obj

    obj%u => scf%u
    obj%e_tb => scf%h%e_tb
    obj%a_tb => scf%h%a_tb
    obj%q => scf%q
    obj%scf => scf

    call obj%initialize()
  end function constructor

  subroutine destructor(obj)
    type(forces) :: obj

    if(allocated(obj%f_at)) deallocate(obj%f_at)
  end subroutine destructor

  function build_forces_rho_per_point(obj,ik,isl) result(rhodens)
    ! INPUT
    class(forces),intent(in) :: obj
    integer,intent(in) :: ik,isl
    ! OUTPUT
    complex(rp),dimension(2,obj%scf%h%nh,obj%scf%h%nh) :: rhodens
    ! LOCAL
    real(rp) :: w
    ! Energy filling f(E_n-E_f) and eigenvalues 
    real(rp),dimension(0:(obj%scf%k%nx)*(obj%scf%h%nh)*(obj%a_tb%nsl)) :: en_k 
    real(rp),dimension((obj%scf%k%nx)*(obj%scf%h%nh)*(obj%a_tb%nsl)) :: f_k
    complex(rp),dimension(2,obj%scf%h%nh,obj%scf%h%nh) :: v_k
    complex(rp) ::  temp ! temporary variable
    
    integer :: imat,jmat,ia1,ia2,ie1,ie2,io1,io2,imat_spin,jmat_spin
    integer :: kmat,kmat2,nn,ispin,jspin
    integer :: nh ! hamiltonian matrix dimension
    
    rhodens = cmplx(0.0_rp,0.0_rp,kind=rp)
    w = obj%scf%k%w(ik) ! get the weight of k from the mesh object
    f_k = obj%scf%en%f_k   ! get f_k from energy object
    en_k = obj%scf%en%en_k  ! get en_k from energy object
    nh = obj%scf%h%nh ! get the hamiltonian matrix dimension from the hamiltonian_tb object
    v_k = obj%scf%h%build_v_k(ik,isl)


    select case(obj%a_tb%ns)
    case(1,2) 
    do ia1=1,obj%a_tb%na
      ie1 = obj%a_tb%ia2ie(ia1)
      do io1=1,obj%e_tb%no(ie1)
        do ia2=1,obj%a_tb%na          
          ie2 = obj%a_tb%ia2ie(ia2) 
          do io2=1,obj%e_tb%no(ie2)
            imat = obj%scf%h%iaos2ih(ia1,io1,1)
            jmat = obj%scf%h%iaos2ih(ia2,io2,1) 
            do kmat=1,nh
              kmat2 = isl+(kmat-1)*obj%a_tb%ns    
              nn = ik+(kmat2-1)*(obj%scf%k%nx)
              temp = (obj%a_tb%g_s)*w*f_k(nn)*v_k(1,imat,kmat)*conjg(v_k(1,jmat,kmat))
              rhodens(1,imat,jmat) = rhodens(1,imat,jmat) + temp
              rhodens(2,imat,jmat) = rhodens(2,imat,jmat) + en_k(nn) * temp
            end do
          end do
        end do
      end do
    end do
    case(4)
    do ia1=1,obj%a_tb%na
      ie1 = obj%a_tb%ia2ie(ia1)
      do io1=1,obj%e_tb%no(ie1)
        do ia2=1,obj%a_tb%na
          ie2 = obj%a_tb%ia2ie(ia2)
          do io2=1,obj%e_tb%no(ie2)
            do ispin=1,2            
              do jspin=1,2
              jmat_spin = obj%scf%h%iaos2ih(ia2,io2,jspin)
              imat_spin = obj%scf%h%iaos2ih(ia1,io1,ispin)

                do kmat=1,nh
                  nn = ik+(kmat-1)*(obj%scf%k%nx)
                  temp = (obj%a_tb%g_s)*w*f_k(nn)*v_k(1,imat_spin,kmat)*conjg(v_k(1,jmat_spin,kmat))
                  rhodens(1,imat_spin,jmat_spin) = rhodens(1,imat_spin,jmat_spin) + temp
                  rhodens(2,imat_spin,jmat_spin) = rhodens(2,imat_spin,jmat_spin) + en_k(nn) * temp
                end do
              end do  
            end do
          end do
        end do
      end do
    end do
    end select

  end function build_forces_rho_per_point

  !> Compute the total atomic forces by summing on k-points and spin states
  subroutine calculate_forces(obj)
    class(forces),intent(inout) :: obj
    integer :: ik,isl

    ! loops on kpoints and spin index
    ! Atomic forces in obj%f_at are accumulated with these loops
    do ik=1,obj%scf%k%nx
      do isl=1,obj%a_tb%nsl
        call obj%calculate_forces_per_point(ik,isl)
      end do
    end do

  end subroutine calculate_forces

  !> Compute the atomic forces of a given k-point and spin state
  subroutine calculate_forces_per_point(obj,ik,isl)
    ! INPUT
    class(forces),intent(inout) :: obj ! object where the forces are
    integer,intent(in) :: ik,isl
    ! LOCAL
    real(rp),dimension(obj%a_tb%na,0:obj%a_tb%nn_max,obj%e_tb%no_max,obj%e_tb%no_max,3) :: d_B ! Derivative of the hopping matrix 
    real(rp),dimension(obj%a_tb%na,0:obj%a_tb%nn_max,obj%e_tb%no_max,obj%e_tb%no_max,3) :: d_S ! Derivative of the overlap matrix
    real(rp),dimension(obj%a_tb%na,obj%e_tb%no_max,3) :: d_en_intra ! Diagonal part of the Hamiltonian calculated using the NRL SK
    
    real(rp) :: ff ! dummy variable to save the forces before saving it in f_at
    real(rp) :: hdvov ! half the delta_v_ov
    real(rp) :: d_rho,r,r_0,r_l,f_cut,t1,t2,rr(3)

    complex(rp),dimension(obj%a_tb%na,obj%a_tb%nn_max,obj%a_tb%nsp) :: c_k
    complex(rp),dimension(2,obj%scf%h%nh,obj%scf%h%nh) :: rhodens
    complex(rp),dimension(obj%scf%h%nh,obj%scf%h%nh,3,obj%a_tb%na) :: d_B_k, d_S_k ! Projections in k-space
    complex(rp),dimension(obj%a_tb%na,3,obj%a_tb%ns) :: delta_v_ov

    integer :: ia1,ia2,io,ie1,ie2,ix,nn,ia2prev,ispin ! various indices
    integer :: nk ! the total number of k-points
    integer,dimension(obj%e_tb%no_max) :: imat, jmat,imat_spin,jmat_spin ! Matrix index

    d_B = 0.0_rp ; d_S = 0.0_rp
    d_B_k = cmplx(0.0_rp,0.0_rp,kind=rp) ; d_S_k = cmplx(0.0_rp,0.0_rp,kind=rp)
    
    d_en_intra = 0.0_rp ; ff = 0.0_rp
    c_k = cmplx(0.0_rp,0.0_rp,kind=rp)
    rhodens = cmplx(0.0_rp,0.0_rp,kind=rp)
    imat = 0.0_rp ; jmat = 0.0_rp
  
    d_B = obj%a_tb%build_d_b_r() ! calculate the derivative of the hopping matrix
    d_S = obj%a_tb%build_d_s_r() ! calculate the derivative of the overlap matrix
    d_en_intra = obj%a_tb%build_d_en_intra() ! calculate the diagonal elements of the Hamiltonian
    
    c_k = obj%a_tb%build_c_k(obj%scf%k%x(ik,:)) ! calculate the c_k in order to calculate 
                                                ! the projection of d_B and d_S into the k-space

    rhodens = obj%build_forces_rho_per_point(ik,isl) ! calculating both \rho and \rho^{\epsilon}
    nk = (obj%scf%k%nx)*(obj%scf%h%nh)*(obj%a_tb%nsl) ! compute the dimension from the content of known objects 

    ! putting together the d_B and its diagonal terms
    do ix=1,3
      do ia1=1,obj%a_tb%na
        ie1 = obj%a_tb%ia2ie(ia1)
        do io=1,obj%e_tb%no(ie1)
          d_B(ia1,0,io,io,ix) = d_en_intra(ia1,io,ix)
        end do
      end do
    end do

    ! projections in the k-space
    d_B_k = obj%scf%h%build_d_projection_k(d_B,c_k)  
    d_S_k = obj%scf%h%build_d_projection_k(d_S,c_k)
    
    delta_v_ov = obj%scf%h%delta_v_lcn + obj%scf%h%delta_v_pen ! self-consistence  

    ! Calculating the components of the forces f_at(1:na,1:3)
    select case(obj%a_tb%ns)
    case(1,2)
    do ix=1,3
      do ia1=1,obj%a_tb%na
        ia2prev=0
        ff=0.0_rp
        ie1 = obj%a_tb%ia2ie(ia1)
        imat(:) = obj%scf%h%iaos2ih(ia1,:,1)
        do nn=1,obj%a_tb%nn(ia1)
          ia2 = obj%a_tb%ian2ia(ia1,nn)
          ie2 = obj%a_tb%ia2ie(ia2)
          jmat(:) = obj%scf%h%iaos2ih(ia2,:,1)
          if (ia2/=ia2prev.and.ia2/=ia1) then
          ! half delta_v_ov
          hdvov = (delta_v_ov(ia2,1,1)+delta_v_ov(ia1,1,1))*0.5_rp
          ff = ff + &
          sum( &
          d_B_k(imat(:),jmat(:),ix,ia1)*conjg(rhodens(1,imat(:),jmat(:))) - &
          d_S_k(imat(:),jmat(:),ix,ia1)*conjg(rhodens(2,imat(:),jmat(:))) &
          ) + &
          sum( &
          d_S_k(imat(:),jmat(:),ix,ia1)*conjg(rhodens(1,imat(:),jmat(:))) &
          ) * hdvov      
          ! plus the CC
          ff = ff + &
          sum( &
          d_B_k(jmat(:),imat(:),ix,ia1)*conjg(rhodens(1,jmat(:),imat(:))) - &
          d_S_k(jmat(:),imat(:),ix,ia1)*conjg(rhodens(2,jmat(:),imat(:))) &
          ) + &
          sum( &
          d_S_k(jmat(:),imat(:),ix,ia1)*conjg(rhodens(1,jmat(:),imat(:))) &
          ) * hdvov
          ia2prev=ia2
          end if
        end do ! end of nn loop
        ! Compute the forces for a single k-point
        obj%f_at(ia1,ix) = obj%f_at(ia1,ix) - ff
      end do ! end of ia1 loop
    end do ! end of ix loop
    !  call obj%print_forces(6)
  
    do ix=1,3
      do ia1=1,obj%a_tb%na
        ie1 = obj%a_tb%ia2ie(ia1)
        imat(:) = obj%scf%h%iaos2ih(ia1,:,1)
        do nn=1,obj%a_tb%nn(ia1)
          ia2 = obj%a_tb%ian2ia(ia1,nn)
          ie2 = obj%a_tb%ia2ie(ia2)
          jmat(:) = obj%scf%h%iaos2ih(ia2,:,1)
          rr(:) = obj%a_tb%rn(ia1,nn,:)
          r = norm2(rr(:))
          r_0 = (obj%e_tb%r_0(ie2) + obj%e_tb%r_0(ie1))/2._rp
          r_l = (obj%e_tb%r_l(ie2) + obj%e_tb%r_l(ie1))/2._rp
          f_cut = fermi_function(r-r_0,1._rp/r_l)

          t1 = ((obj%e_tb%p(ie2,1)+obj%e_tb%p(ie1,1))/2)**2
          t2 = f_cut*exp(-t1*r)*(t1+(1._rp-f_cut)/r_l)

          d_rho = (rr(ix)/r)*t2

          do io=1,obj%e_tb%no(ie2)
            obj%f_at(ia1,ix) = obj%f_at(ia1,ix) - d_rho*(d_B(ia1,0,io,io,ix)*rhodens(1,imat(io),imat(io))+&
                                                         d_B(ia2,0,io,io,ix)*rhodens(1,jmat(io),jmat(io)))
          end do
        end do
      end do
    end do 

    case(4)
    do ix=1,3
      do ispin=1,2
        do ia1=1,obj%a_tb%na
          ia2prev=0
          ff=0.0_rp
          ie1 = obj%a_tb%ia2ie(ia1)
          do nn=1,obj%a_tb%nn(ia1)
            ia2 = obj%a_tb%ian2ia(ia1,nn)
            ie2 = obj%a_tb%ia2ie(ia2)
            imat_spin(:) = obj%scf%h%iaos2ih(ia1,:,ispin)
            jmat_spin(:) = obj%scf%h%iaos2ih(ia2,:,ispin)
            if (ia2/=ia2prev.and.ia2/=ia1) then
            ! half delta_v_ov
            hdvov = (delta_v_ov(ia2,1,1)+delta_v_ov(ia1,1,1))*0.5_rp
            ff = ff + &
            sum( &
            d_B_k(imat_spin(:),jmat_spin(:),ix,ia1)*conjg(rhodens(1,imat_spin(:),jmat_spin(:))) - &
            d_S_k(imat_spin(:),jmat_spin(:),ix,ia1)*conjg(rhodens(2,imat_spin(:),jmat_spin(:))) &
            ) + &
            sum( &
            d_S_k(imat_spin(:),jmat_spin(:),ix,ia1)*conjg(rhodens(1,imat_spin(:),jmat_spin(:))) &
            ) * hdvov
            ! plus the CC
            ff = ff + &
            sum( &
            d_B_k(jmat_spin(:),imat_spin(:),ix,ia1)*conjg(rhodens(1,jmat_spin(:),imat_spin(:))) - &
            d_S_k(jmat_spin(:),imat_spin(:),ix,ia1)*conjg(rhodens(2,jmat_spin(:),imat_spin(:))) &
            ) + &
            sum( &
            d_S_k(jmat_spin(:),imat_spin(:),ix,ia1)*conjg(rhodens(1,jmat_spin(:),imat_spin(:))) &
            ) * hdvov
            ia2prev=ia2  
            end if
          end do ! end of nn loop      
          ! Compute the forces for a single k-point
          obj%f_at(ia1,ix) = obj%f_at(ia1,ix) - ff
        end do ! end of ia1 loop
      end do ! end of spin loop
    end do ! end of ix loop
    
    do ix=1,3
      do ispin=1,2
        do ia1=1,obj%a_tb%na
          ie1 = obj%a_tb%ia2ie(ia1)
          do nn=1,obj%a_tb%nn(ia1)
            ia2 = obj%a_tb%ian2ia(ia1,nn)
            ie2 = obj%a_tb%ia2ie(ia2)
            imat_spin(:) = obj%scf%h%iaos2ih(ia1,:,ispin)
            jmat_spin(:) = obj%scf%h%iaos2ih(ia2,:,ispin)
            rr(:) = obj%a_tb%rn(ia1,nn,:)
            r = norm2(rr(:))
            r_0 = (obj%e_tb%r_0(ie2) + obj%e_tb%r_0(ie1))/2._rp
            r_l = (obj%e_tb%r_l(ie2) + obj%e_tb%r_l(ie1))/2._rp
            f_cut = fermi_function(r-r_0,1._rp/r_l)
  
            t1 = ((obj%e_tb%p(ie2,1)+obj%e_tb%p(ie1,1))/2)**2
            t2 = f_cut*exp(-t1*r)*(t1+(1._rp-f_cut)/r_l)
  
            d_rho = (rr(ix)/r)*t2
  
            do io=1,obj%e_tb%no(ie2)
              obj%f_at(ia1,ix) = obj%f_at(ia1,ix) - d_rho*(d_B(ia1,0,io,io,ix)*rhodens(1,imat_spin(io),imat_spin(io))+&
                                                           d_B(ia2,0,io,io,ix)*rhodens(1,jmat_spin(io),jmat_spin(io)))
            end do
          end do
        end do
      end do
    end do
    end select

  end subroutine calculate_forces_per_point

  subroutine initialize(obj)
    class(forces),intent(inout) :: obj
  
    if(allocated(obj%f_at)) deallocate(obj%f_at)
    allocate(obj%f_at(obj%a_tb%na,3))
    obj%f_at = 0.0_rp
    obj%computed = .false.
  end subroutine initialize

  !> Print the forces on screen
  !> Forces are computed internally in hau and then printed to the choosen units
  subroutine print_forces(obj,unit)
    class(forces),intent(inout) :: obj
    integer,intent(in),optional :: unit
    integer :: ia
    real(rp) :: factor

    factor = obj%u%convert_energy('from','hau')/obj%u%convert_length('from','hau')
    do ia=1,obj%a_tb%na
      write(unit,'(a)') ' f(' // int2str(ia) // ',:) = ' &
           // fixedreal2str(obj%f_at(ia,1)*factor) // ', ' &
           // fixedreal2str(obj%f_at(ia,2)*factor) // ', ' &
           // fixedreal2str(obj%f_at(ia,3)*factor)
      call dynamol_flush(unit=unit)
    end do 
  end subroutine print_forces
  
  !> Read object in text format from file (default: 'in_forces.txt')
  subroutine read_txt(obj,file)
    class(forces),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen

    ! Namelist variables
    logical :: computed
    ! Namelist
    namelist /forces/ computed

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_forces.txt'
    end if

    inquire(unit=10, opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'forces%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'forces%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    computed = .false.
    read(10,nml=forces)
    
    !call check_computed(computed)

    obj%computed = computed

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_forces.txt')
  subroutine write_txt(obj,file,unit)
    class(forces), intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    logical  :: computed ! true if forces are computed 
    ! Namelist
    namelist /forces/ computed

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_forces.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    computed = obj%computed

    write(unit_rt,nml=forces)
    call dynamol_flush(unit_rt)
    if(.not. present(unit)) then
      close(unit_rt)
    end if
    
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_forces.txt'), if tag (default: .true.) the namelist opening 
  !> and closing tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(forces),intent(in) :: obj
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
      file_rt = 'out_forces.txt'
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
      write(unit_rt,'(a)') '&forces'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('computed')
        write(unit_rt,'(a)') ' computed = ' // log2str(obj%computed)
      end select
    end do

    if(tag_rt) then
      write(unit_rt,'(a)') ' /'
    end if
    
    call dynamol_flush(unit_rt)
    
    if(.not. present(unit)) then
      close(unit_rt)
    end if
    
  end subroutine write_txt_formatted
    
end module forces_mod
