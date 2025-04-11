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
!  atom.f90
!  TBKOSTER
module atom_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use element_mod
  use lattice_mod
  use math_mod, only: deg2rad, epsilon, i_unit, pi, rad2deg, sph2cart
  use precision_mod, only: rp
  use string_mod, only: sl, int2str, lower, real2str, fixedreal2str, unique_str, TBKOSTER_flush
  use units_mod
  use constant_mod, only: a_0
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'ns', &
   'na', &
   'ntag', &
   'stag', &
   'tag', &
   'ia2ie', &
   'pbc', &
   'k_spiral', &
   'r_coord', &
   'm_coord', &
   'r', &
   'p', &
   'm', &
   'lambda_pen' &
   ]

  type,public :: atom
    !> Units
    class(units),pointer :: u
    !> Elements
    class(element),pointer :: e
    !> Real space lattice
    class(lattice),pointer :: l_r
    !> Reciprocal space lattice
    class(lattice),pointer :: l_k

    !> @defgroup Spin Spin-related variables
    !> @{

    !> Spin polarization ; options:
    !>  'unpolarized' : spin-unpolarized (default)
    !>  'collinear    : spin-polarized collinear
    !>  'noncollinear': spin-polarized noncollinear
    !character(len=12) :: spin_polarization
    !> Number of spins, options:
    !>  1: spin-unpolarized
    !>  2: spin-polarized collinear
    !>  4: spin-polarized noncollinear
    integer :: ns
    !> Number of spin loops (ns=1,4->1, ns=2->2)
    integer :: nsl
    !> Spin polarization (ns=1->1, ns=2,4->2)
    integer :: nsp
    !> Spin degeneracy (ns=1->2, ns=2,4->1)
    integer :: g_s
    !> 2spin/2spin-to-4spin index
    integer,dimension(2,2) :: iss2is
    !> @}

    !> @defgroup Atom Atom-related variables
    !> @{

    !> Number of atoms
    integer :: na
    !> Number of tags
    integer :: ntag
    !> Strides between tags
    integer,dimension(:),allocatable :: stag
    !> Tags
    character(len=sl),dimension(:),allocatable :: tag
    !> Atom-to-element index
    integer,dimension(:),allocatable :: ia2ie
    !> Number of electrons
    real(rp) :: nel
    !> @}

    !> Periodic boundary conditions
    integer,dimension(3) :: pbc
    !> Spin spiral vector
    real(rp),dimension(3) :: k_spiral

    !> @defgroup Dynamics Dynamics-related variables
    !> @{

    !> Position coordinates; options:
    !>  'cartesian'
    !>  'direct' (default)
    character(len=9) :: r_coord
    !> Positions, (na,3)
    real(rp),dimension(:,:),allocatable :: r
    !> Momentums, (na,3)
    real(rp),dimension(:,:),allocatable :: p
    !> Expression of the magnetization; options:
    !> 'cartesian'
    !> 'spherical' (default)
    character(len=9) :: m_coord
    !> Magnetizations in spherical coordinates (r,theta,phi) (radial, polar and
    !> azimuthal component), (na,3)
    real(rp),dimension(:,:),allocatable :: m, m_pen
    !> @}

    !> @defgroup Neighbours Neighbours-related variables
    !> @{

    !> Number of neighbours
    integer,dimension(:),allocatable :: nn                                  
    !> Maximum number of neighbours
    integer :: nn_max
    !> Atom/neighbour-to-atom index
    integer,dimension(:,:),allocatable :: ian2ia
    !> Atom/periodic-boundary-to-neighbour index
    integer,dimension(:,:,:,:,:),allocatable :: iapbc2in
    !> Neighbour position
    real(rp),dimension(:,:,:),allocatable :: rn        
    !> @}

    !> @defgroup Magnetic_penalization Magnetic penalization related variables
    !> @{

    !> Magnetic penalization.\n
    !> See \cite Gebauer1999 page 45-48.\n
    !> To impose a magnetic constraint on a system, an extra term is added to the
    !> total energy functional \f$ E_{\mathrm{tot}}[N,\mathbf{M}] \f$ in the form
    !> of a penalization energy functional.\n
    !> The penalization is chosen such that it is zero if the constraint is
    !> fulfilled and large and positive otherwise. In this way, the constraint is
    !> imposed approximately by energetically penalizing any deviation from a
    !> certain condition.\n
    !> When trying to impose a given magnetization \f$ \mathbf{M}_{\mathrm{pen}}
    !> \f$, the penalization can be chosen to be proportional to the mean square
    !> deviation of the actual magnetization \f$ \mathbf{M} \f$ from its target
    !> \f$ \mathbf{M}_{\mathrm{pen}} \f$.\n
    !> The new energy functional reads:
    !> \f{equation}{
    !>   E_{\lambda}[N,\mathbf{M}] = E_{\mathrm{tot}}[N,\mathbf{M}] + \lambda \int
    !>   \mathrm{d}\mathbf{r} (\mathbf{M}(\mathbf{r}) - \mathbf{M}_{\mathrm{pen}}
    !>   (\mathbf{r}))^2
    !> \f}
    !> where \f$ \lambda \f$ is a large positive multiplier.\n
    !> For a finite \f$ \lambda \f$, the constraint is not imposed exactly as the
    !> difference \f$ \mathbf{M}(\mathbf{r}) - \mathbf{M}_{\mathrm{pen}}
    !> (\mathbf{r}) \f$ remains finite too.\n
    !> In the limit \f$ \lambda \rightarrow \infty \f$, the difference \f$
    !> \mathbf{M}(\mathbf{r}) - \mathbf{M}_{\mathrm{pen}} (\mathbf{r}) \f$ tends
    !> to zero. This however comes at the cost of a slower convergence.\n
    !> In practical calculations, \f$ \lambda \f$ must therefore be chosen
    !> sufficiently large such that \f$ \mathbf{M} \f$ is close enough to \f$
    !> \mathbf{M}_{\mathrm{pen}} \f$, while still being as small as possible in
    !> order to allow fast convergence of the self-consistent cycle.
    !> Multiplier for the magnetic penalization
    real(rp),dimension(:),allocatable   :: lambda_pen
    !> Effective magnetic induction from the magnetic penalization, 'na,3)
    real(rp),dimension(:,:),allocatable :: b_pen
    !> @}

    ! On-site potential difference
    !real(rp),dimension(:),allocatable :: delta_v

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: build_c_k
    procedure :: calculate_neighbours
    procedure :: calculate_nel
    !procedure :: calculate_ns
    procedure :: calculate_spin
    procedure :: read_txt
    procedure :: write_lammps
    procedure :: write_txt
    procedure :: write_txt_formatted
    procedure :: write_xyz
  end type atom

  ! Constructor
  interface atom
    procedure :: constructor
  end interface atom

contains
  function constructor(e,l_r,l_k) result(obj)
    class(element),target,intent(in) :: e
    class(lattice),target,intent(in) :: l_r
    class(lattice),target,intent(in) :: l_k
    type(atom) :: obj

    obj%u => e%u
    obj%e => e
    obj%l_r => l_r
    obj%l_k => l_k

    obj%iss2is = reshape((/1,4,3,2/),(/2,2/))
  end function constructor

  subroutine destructor(obj)
    type(atom) :: obj

    if(allocated(obj%stag))       deallocate(obj%stag)
    if(allocated(obj%tag))        deallocate(obj%tag)
    if(allocated(obj%ia2ie))      deallocate(obj%ia2ie)
    if(allocated(obj%r))          deallocate(obj%r)
    if(allocated(obj%p))          deallocate(obj%p)
    if(allocated(obj%m))          deallocate(obj%m)
    if(allocated(obj%m_pen))      deallocate(obj%m_pen)
    if(allocated(obj%nn))         deallocate(obj%nn)
    if(allocated(obj%ian2ia))     deallocate(obj%ian2ia)
    if(allocated(obj%iapbc2in))     deallocate(obj%iapbc2in)
    if(allocated(obj%rn))         deallocate(obj%rn)
    if(allocated(obj%lambda_pen)) deallocate(obj%lambda_pen)
    if(allocated(obj%b_pen))      deallocate(obj%b_pen)
  end subroutine destructor

  function build_ia2ie(tag) result(ia2ie)
    character(len=*),dimension(:),intent(in) :: tag
    integer,dimension(size(tag)) :: ia2ie
    character(len=len(tag)),dimension(:),allocatable :: cdummy
    integer,dimension(:),allocatable :: idummy

    ! copy array of strings of exactly 2 characters to avoid gfortran warning
    character(len=2),dimension(:),pointer :: tug
    allocate(tug(size(tag)))
    tug(:)=tag(:)(1:2)
    call unique_str(tug,cdummy,idummy,ia2ie)
    deallocate(tug)
  end function build_ia2ie

  function build_c_k(obj,k) result(c_k)
    ! INPUT
    class(atom),intent(in) :: obj
    real(rp),dimension(3),intent(in) :: k
    ! OUTPUT
    complex(rp),dimension(obj%na,obj%nn_max,obj%nsp) :: c_k
    !	LOCAL
    real(rp) :: k_temp(3),RR_cart(3),RR_unit(3),prodscalR
    integer :: ia,in,is

    c_k = cmplx(0.0_rp,0.0_rp,kind=rp)

    select case(obj%ns)
    case(1,2)
      do ia=1,obj%na
        do in=1,obj%nn(ia)
          RR_cart = obj%rn(ia,in,:)
          RR_unit(1) = dot_product(RR_cart,obj%l_k%v(1,:))
          RR_unit(2) = dot_product(RR_cart,obj%l_k%v(2,:))
          RR_unit(3) = dot_product(RR_cart,obj%l_k%v(3,:))
          ProdscalR  = dot_product(k,RR_unit)
          do is=1,obj%ns
            c_k(ia,in,is) = exp(2*pi*i_unit*prodscalR)
          end do
        end do
      end do
    case(4)
      do ia=1,obj%na
        do in=1,obj%nn(ia)
          RR_cart = obj%rn(ia,in,:)
          RR_unit(1) = dot_product(RR_cart,obj%l_k%v(1,:))
          RR_unit(2) = dot_product(RR_cart,obj%l_k%v(2,:))
          RR_unit(3) = dot_product(RR_cart,obj%l_k%v(3,:))
          do is=1,2
            k_temp = k + (2*is-3)*0.5_rp*obj%k_spiral
            ProdscalR  = dot_product(k_temp,RR_unit)
            c_k(ia,in,is) = exp(2*pi*i_unit*prodscalR)
          end do
        end do
      end do
    end select
  end function build_c_k

  subroutine calculate_neighbours(obj,r_max,type)
    ! INPUT
    class(atom),intent(inout) :: obj
    real(rp),intent(inout) :: r_max ! cutoff distance to compute the neighbors 
    character(len=3):: type
    ! LOCAL
    integer,dimension(:),allocatable :: nn ! number of atoms of atom i
    integer,dimension(:,:),allocatable :: ian2ia ! list of atoms of a given atom
    integer,dimension(:,:,:,:,:),allocatable :: iapbc2in ! list of atoms of a given atom
    real(rp),dimension(:,:,:),allocatable :: rn ! distance vector between atoms
    integer :: ia1,ia2,in,ip1,ip2,ip3               
    integer,dimension(3) :: x            
    real(rp),dimension(3) :: r
    real(rp) :: r_max_2,r_2

  
    select case(trim(type))
    case('nrl')
    ! compute their numbers
    r_max_2 = r_max*r_max

  case('mod','wan')

    r_max_2=0.0
    do ia1=1,obj%na
      do ia2=1,obj%na
        do ip1=1,2*obj%pbc(1)+1
          do ip2=1,2*obj%pbc(2)+1
            do ip3=1,2*obj%pbc(3)+1
              x = (/ip1,ip2,ip3/) - obj%pbc - 1
              r = matmul(x,obj%l_r%v)
              r = r + obj%r(ia2,:) - obj%r(ia1,:)
              r_2 = dot_product(r,r)
              if(r_2>r_max_2) then
                r_max_2=r_2
              end if
            end do
          end do
        end do
      end do
    end do


  
 !   do ip1=1,2*obj%pbc(1)+1
 !     do ip2=1,2*obj%pbc(2)+1
 !       do ip3=1,2*obj%pbc(3)+1
 !         x = (/ip1,ip2,ip3/) - obj%pbc - 1
 !         r = matmul(x,obj%l_r%v)
 !         r_2 = dot_product(r,r)
 !         if(r_2>r_max_2) then
 !           r_max_2=r_2
 !         end if
 !       end do
 !     end do
 !   end do
    r_max_2=r_max_2+4*epsilon
    r_max=sqrt(r_max_2)
    end select

    allocate(nn(obj%na))


    nn = 0
    do ia1=1,obj%na
      do ia2=1,obj%na
        do ip1=1,2*obj%pbc(1)+1
          do ip2=1,2*obj%pbc(2)+1
            do ip3=1,2*obj%pbc(3)+1
              x = (/ip1,ip2,ip3/) - obj%pbc - 1
              r = matmul(x,obj%l_r%v)
              r = r + obj%r(ia2,:) - obj%r(ia1,:)
              r_2 = dot_product(r,r)
              if(r_2<r_max_2 .and. r_2>epsilon) then
                nn(ia1) = nn(ia1)+1
              end if
            end do
          end do
        end do
      end do
    end do

    obj%nn_max = maxval(nn)
    allocate(ian2ia(obj%na,obj%nn_max))
    allocate(iapbc2in(obj%na,obj%na,2*obj%pbc(1)+1,2*obj%pbc(2)+1,2*obj%pbc(3)+1))
    allocate(rn(obj%na,obj%nn_max,3))

    ! then compute their distances
    do ia1=1,obj%na
      in = 0
      do ia2=1,obj%na
        do ip1=1,2*obj%pbc(1)+1
          do ip2=1,2*obj%pbc(2)+1
            do ip3=1,2*obj%pbc(3)+1
              x = (/ip1,ip2,ip3/) - obj%pbc - 1
              r = matmul(x,obj%l_r%v)
              r = r + obj%r(ia2,:) - obj%r(ia1,:)
              r_2 = dot_product(r,r)
              if(r_2<r_max_2 .and. r_2>epsilon) then
                in = in+1
                ian2ia(ia1,in) = ia2
                iapbc2in(ia1,ia2,ip1,ip2,ip3)=in
                rn(ia1,in,:) = r
              end if
              if(r_2<=epsilon) then
                iapbc2in(ia1,ia2,ip1,ip2,ip3)=0
              end if
            end do
          end do
        end do
      end do
    end do


    call move_alloc(nn,obj%nn)
    call move_alloc(ian2ia,obj%ian2ia)
    call move_alloc(iapbc2in,obj%iapbc2in)
    call move_alloc(rn,obj%rn)

  end subroutine calculate_neighbours

  subroutine calculate_nel(obj)
    class(atom),intent(inout) :: obj
    integer :: ia,ie

    obj%nel = 0.0_rp
    do ia=1,obj%na
      ie = obj%ia2ie(ia)
      obj%nel = obj%nel + obj%e%q(ie)
    end do
  end subroutine calculate_nel

  ! subroutine calculate_ns(obj)
  !   class(atom),intent(inout) :: obj
  !
  !   select case(obj%spin_polarization)
  !   case('unpolarized')
  !     obj%ns = 1
  !   case('collinear')
  !     obj%ns = 2
  !   case('noncollinear')
  !     obj%ns = 4
  !   end select
  ! end subroutine calculate_ns

  subroutine calculate_spin(obj)
    class(atom),intent(inout) :: obj

    select case(obj%ns)
    case(1)
      obj%nsl = 1
      obj%nsp = 1
      obj%g_s = 2
    case(2)
      obj%nsl = 2
      obj%nsp = 2
      obj%g_s = 1
    case(4)
      obj%nsl = 1
      obj%nsp = 2
      obj%g_s = 1
    end select
  end subroutine calculate_spin

  !> Check the validity of the listing
  subroutine check_listing(listing)
    character(len=*),intent(in) :: listing

    if(listing /= 'by_atom' .and. listing /= 'by_tag') then
      write(error_unit,*) 'atom%check_listing(): listing must be one of: &
       &''by_atom'', ''by_tag'''
      error stop
    end if
  end subroutine check_listing

  !> Check the validity of the number of tags
  subroutine check_ntag(na,ntag)
    integer,intent(in) :: na, ntag

    if(ntag < 1 .or. ntag > na) then
      write(error_unit,*) 'atom%check_ntag(): ntag must verify &
       &1 <= atom%ntag <= atom%na'
      error stop
    end if
  end subroutine check_ntag

  !> Check the validity of the position coordinates
  subroutine check_r_coord(r_coord)
    character(len=*),intent(in) :: r_coord

    if(r_coord /= 'cartesian' .and. r_coord /= 'direct') then
      write(error_unit,*) 'atom%check_r_coord(): r_coord must be one of: &
       &''cartesian'', ''direct'''
      error stop
    end if
  end subroutine check_r_coord

  subroutine check_direction(direction)
    character(len=*), intent(in) :: direction

    if(direction /= 'from' &
     .and. direction /= 'to') then
      write(error_unit,*) 'units%check_direction(): direction must be one of: &
       &''from'', ''to'''
      error stop
    end if
  end subroutine check_direction

  !> Check the validity of the magnetization coordinates
  subroutine check_m_coord(m_coord)
    character(len=*),intent(in) :: m_coord

    if(m_coord /= 'cartesian' .and. m_coord /= 'spherical') then
      write(error_unit,*) 'atom%check_m_coord(): m_coord must be one of: &
       &''cartesian'', ''spherical'''
      error stop
    end if
  end subroutine check_m_coord

  !> Check the validity of the spin polarization
  ! subroutine check_spin_polarization(spin_polarization)
  !   character(len=*),intent(in) :: spin_polarization
  !
  !   if(.not. (spin_polarization == 'unpolarized' &
  !    .or. spin_polarization == 'collinear' &
  !    .or. spin_polarization == 'noncollinear')) then
  !     write(error_unit,*) 'atom%check_spin_polarization(): &
  !      &spin_polarization must be one of: ''unpolarized'', ''collinear'', &
  !      &''noncollinear'''
  !     error stop
  !   end if
  ! end subroutine check_spin_polarization

  !> Check the validity of the strides between tags
  subroutine check_stag(na,stag)
    integer,intent(in) :: na
    integer,dimension(:),intent(in) :: stag

    if(sum(stag) /= na) then
      write(error_unit,*) 'atom%check_stag(): stag must verify &
       &sum(atom%stag) == atom%na'
      error stop
    end if
  end subroutine check_stag

  subroutine convert_coordinates(obj,direction,m_coord)
    class(atom), intent(inout) :: obj
    character(len=*),intent(in) :: direction, m_coord
    character(len=9) :: coords_1, coords_2
    real(rp) :: m_x,m_y,m_z,m_r,m_theta,m_phi
    integer :: ia

    call check_direction(direction)
    select case(direction)
    case('from')
      coords_1 = m_coord
      coords_2 = obj%m_coord
    case('to')
      coords_1 = obj%m_coord
      coords_2 = m_coord
    end select

    select case(trim(coords_1))
    case('cartesian')
    select case(trim(coords_2))
    case('cartesian')
      ! do nothing
    case('spherical')
      do ia = 1,obj%na
        m_x=obj%m(ia,1)
        m_y=obj%m(ia,2)
        m_z=obj%m(ia,3)
        !        write(output_unit,*) "=> here mx my mz", m_x,m_y,m_z
        !        call TBKOSTER_flush(output_unit)
        m_r=sqrt(m_x*m_x+m_y*m_y+m_z*m_z)
        m_theta=acos(m_z/m_r)
        m_phi=atan(m_y/m_x)
        !        write(output_unit,*) "=> here mr mtheta mphi", m_r,m_theta,m_phi
        !        call TBKOSTER_flush(output_unit)
        obj%m(ia,1)=m_r
        obj%m(ia,2)=m_theta
        obj%m(ia,3)=m_phi
      end do
    end select
    case('spherical')
    select case(trim(coords_2))
    case('cartesian')
      do ia = 1,obj%na
        m_r=obj%m(ia,1)
        m_theta=obj%m(ia,2)
        m_phi=obj%m(ia,3)
        !        write(output_unit,*) "=> here mr mtheta mphi ", m_r,m_theta,m_phi
        !        call TBKOSTER_flush(output_unit)
        m_x=m_r*sin(m_theta)*cos(m_phi)
        m_y=m_r*sin(m_theta)*sin(m_phi)
        m_z=m_r*cos(m_theta)
        !        write(output_unit,*) "=> here mx my mz", m_x,m_y,m_z
        !        call TBKOSTER_flush(output_unit)
        obj%m(ia,1)=m_x
        obj%m(ia,2)=m_y
        obj%m(ia,3)=m_z
      end do

    case('spherical')
      ! do nothing
    end select
    end select

  end subroutine convert_coordinates

  subroutine initialize_dyn(na, r, p, m)
    integer,intent(in) :: na
    real(rp),dimension(:,:),allocatable,intent(out) :: r
    real(rp),dimension(:,:),allocatable,intent(out) :: p
    real(rp),dimension(:,:),allocatable,intent(out) :: m

    allocate(r(na,3))
    allocate(p(na,3))
    allocate(m(na,3))

    r = 0.0_rp
    p = 0.0_rp
    m = 0.0_rp
  end subroutine initialize_dyn

  subroutine initialize_ia2ie(ntag,stag,tag,ia2ie)
    integer,intent(in) :: ntag
    integer,dimension(ntag),intent(in) :: stag
    character(len=*),dimension(ntag),intent(in) :: tag
    integer,dimension(:),allocatable,intent(out) :: ia2ie
    integer,dimension(ntag) :: itag2ie
    integer :: ia1, ia2, itag

    allocate(ia2ie(sum(stag)))

    itag2ie = build_ia2ie(tag)
    ia1 = 1
    do itag=1,ntag
      ia2 = ia1+stag(itag)-1
      ia2ie(ia1:ia2) = spread(itag2ie(itag),1,stag(itag))
      ia1 = ia2+1
    end do
  end subroutine initialize_ia2ie

  subroutine initialize_pbc(pbc,k_spiral)
    integer,dimension(3),intent(out) :: pbc
    real(rp),dimension(3),intent(out) :: k_spiral

    pbc = 5
    k_spiral = 0.0_rp
  end subroutine initialize_pbc

  subroutine initialize_pen(na, m_pen,lambda_pen, b_pen)
    integer,intent(in) :: na
    real(rp),dimension(:),allocatable,intent(out) :: lambda_pen
    real(rp),dimension(:,:),allocatable,intent(out) :: b_pen
    real(rp),dimension(:,:),allocatable,intent(out) :: m_pen

    allocate(lambda_pen(na))
    allocate(b_pen(na,3))
    allocate(m_pen(na,3))

    lambda_pen = 0.0_rp
    b_pen = 0.0_rp
    m_pen= 0.0_rp
  end subroutine initialize_pen

  subroutine initialize_tag(na, ntag, stag, tag)
    integer,intent(in) :: na, ntag
    integer,dimension(:),allocatable,intent(out) :: stag
    character(len=sl),dimension(:),allocatable,intent(out) :: tag

    allocate(stag(ntag))
    allocate(tag(ntag))

    if(ntag == 1) then
      stag = na
    elseif(ntag == na) then
      stag = 1
    end if
  end subroutine initialize_tag

  !> Read object in text format from file (default: 'in_atom.txt')
  subroutine read_txt(obj,file)
    class(atom),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    ! Namelist variables
    integer :: ns, na, ntag
    integer,dimension(:),allocatable :: stag
    character(len=sl),dimension(:),allocatable :: tag
    integer,dimension(:),allocatable :: ia2ie
    integer,dimension(3) :: pbc
    real(rp),dimension(3) :: k_spiral
    character(len=9) :: r_coord
    real(rp),dimension(:,:),allocatable :: r
    real(rp),dimension(:,:),allocatable :: p
    character(len=7) :: m_listing
    character(len=9) :: m_coord
    real(rp),dimension(:,:),allocatable :: m,m_pen
    character(len=7) :: lambda_pen_listing
    real(rp),dimension(:)  ,allocatable :: lambda_pen
    real(rp),dimension(:,:),allocatable :: b_pen
    logical :: isopen
    ! Namelist
    namelist /atom/ ns, na, ntag, stag, tag, ia2ie, pbc, k_spiral, r_coord, &
    m_coord, r, p, m_listing, m, lambda_pen_listing, lambda_pen
    ! Local variables
    integer :: ia1, ia2, itag

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_atom.txt'
    end if

    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'atom%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'atom%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    call initialize_pbc(pbc,k_spiral)
    r_coord = 'direct'
    m_coord = 'spherical'
    m_listing = 'by_atom'
    lambda_pen_listing = 'by_atom'
    allocate(stag(0))
    read(10,nml=atom,iostat=iostatus)
    deallocate(stag)
    call check_ntag(na,ntag)
    call initialize_tag(na,ntag,stag,tag)
    call initialize_dyn(na,r,p,m)
    call initialize_pen(na,m_pen,lambda_pen,b_pen)
    allocate(ia2ie(0))
    rewind(10)
    read(10,nml=atom,iostat=iostatus)
    deallocate(ia2ie)
    call check_stag(na,stag)
    call initialize_ia2ie(ntag,stag,tag,ia2ie)
    rewind(10)
    read(10,nml=atom)
    r_coord = trim(lower(r_coord))
    call check_r_coord(r_coord)
    if(r_coord== 'cartesian') then
      r = r * obj%u%convert_length('to','hau')
    else !if(trim(r_coord) == 'direct') then
      r_coord = 'cartesian'
      r = obj%l_r%dir2cart(r)
    end if
    m_coord = trim(lower(m_coord))
    call check_m_coord(m_coord)
    m_listing = trim(lower(m_listing))
    call check_listing(m_listing)
    if(m_listing == 'by_tag') then
      ia2 = na
      do itag=ntag,1,-1
        ia1 = ia2-stag(itag)+1
        m(ia1:ia2,:) = spread(m(itag,:),1,stag(itag))
        ia2 = ia1-1
      end do
    end if
    lambda_pen_listing = trim(lower(lambda_pen_listing))
    call check_listing(lambda_pen_listing)
    if(lambda_pen_listing == 'by_tag') then
      ia2 = na
      do itag=ntag,1,-1
        ia1 = ia2-stag(itag)+1
        lambda_pen(ia1:ia2) = spread(lambda_pen(itag),1,stag(itag))
        ia2 = ia1-1
      end do
    end if

    obj%ns = ns
    obj%na = na
    obj%ntag = ntag
    call move_alloc(stag,obj%stag)
    call move_alloc(tag,obj%tag)
    call move_alloc(ia2ie,obj%ia2ie)
    obj%pbc = pbc
    obj%k_spiral = k_spiral
    obj%r_coord = r_coord

    call move_alloc(r,obj%r)
    call move_alloc(p,obj%p)
    call move_alloc(m,obj%m)
    call move_alloc(m_pen,obj%m_pen)
    if ((ns == 4).and.(m_coord == 'spherical')) then
      obj%m(:,2:3) = obj%m(:,2:3) * deg2rad
    end if

    obj%m_coord = 'spherical' ! internal representation is spherical
    if ((ns == 4).and.(m_coord == 'cartesian')) then
      call convert_coordinates(obj,'from','cartesian')
    end if
    obj%m_pen=obj%m
    lambda_pen = lambda_pen * obj%u%convert_energy('to','hau')
    call move_alloc(lambda_pen,obj%lambda_pen)
    call move_alloc(b_pen,obj%b_pen)

    call obj%calculate_nel()
    call obj%calculate_spin()

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> For timestep t, write lammps file to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_lammps.lammpstrj')
  subroutine write_lammps(obj,t,file,unit)
    class(atom), intent(in) :: obj

    real(rp), intent(in) :: t ! current timestep
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit

    integer :: unit_rt
    real(rp) :: lx,xy,xz,ly,yz,lz
    real(rp) :: xlo=0.0_rp,xhi,ylo=0.0_rp,yhi,zlo=0.0_rp,zhi
    real(rp) :: xlo_bound,xhi_bound,ylo_bound,yhi_bound,zlo_bound,zhi_bound
    real(rp), dimension(3) :: m, v_sph
    integer :: ia

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_lammps.lammpstrj'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write',position='append')
    end if

    write(unit_rt,'(a)') 'ITEM: TIMESTEP'
    write(unit_rt,'(a)') real2str(t*obj%u%convert_time('from','hau'))
    write(unit_rt,'(a)') 'ITEM: NUMBER OF ATOMS'
    write(unit_rt,'(a)') int2str(obj%na)
    write(unit_rt,'(a)') 'ITEM: BOX BOUNDS xy xz yz'
    lx = obj%l_r%v(1,1)
    xy = obj%l_r%v(1,2)
    yz = obj%l_r%v(1,3)
    ly = obj%l_r%v(2,2)
    yz = obj%l_r%v(2,3)
    lz = obj%l_r%v(3,3)
    do ia=1, obj%na
      xlo = merge (xlo, obj%r(ia,1), xlo < obj%r(ia,1))
      ylo = merge (ylo, obj%r(ia,2), ylo < obj%r(ia,2))
      zlo = merge (zlo, obj%r(ia,3), zlo < obj%r(ia,3))
    end do
    xhi = xlo + lx
    yhi = ylo + ly
    zhi = zlo + lz
    xlo_bound = xlo + min(0.0_rp, min(xy, min(xz,xy+xz)))
    xhi_bound = xhi + max(0.0_rp, max(xy, max(xz,xy+xz)))
    ylo_bound = ylo + min(0.0_rp, yz)
    yhi_bound = yhi + max(0.0_rp, yz)
    zlo_bound = zlo
    zhi_bound = zhi

    xlo_bound = xlo_bound*obj%u%convert_length('from','hau')
    xhi_bound = xhi_bound*obj%u%convert_length('from','hau')
    ylo_bound = ylo_bound*obj%u%convert_length('from','hau')
    yhi_bound = yhi_bound*obj%u%convert_length('from','hau')
    zlo_bound = zlo_bound*obj%u%convert_length('from','hau')
    zhi_bound = zhi_bound*obj%u%convert_length('from','hau')

    write(unit_rt,'(a)') real2str(xlo_bound)//' '//real2str(xhi_bound)&
                  //' '//real2str(xy)
    write(unit_rt,'(a)') real2str(ylo_bound)//' '//real2str(yhi_bound)&
                  //' '//real2str(xz)
    write(unit_rt,'(a)') real2str(zlo_bound)//' '//real2str(zhi_bound)&
                  //' '//real2str(yz)

    write(unit_rt,'(a)') 'ITEM: ATOMS type x y z vx vy vz'
    do ia=1, obj%na
      v_sph = obj%m(ia,:)
      m = sph2cart(v_sph)
      write(unit_rt, '(a)') int2str(obj%ia2ie(ia))//' '//&
      real2str(obj%r(ia,1)*obj%u%convert_length('from','hau'))//' '//&
      real2str(obj%r(ia,2)*obj%u%convert_length('from','hau'))//' '//&
      real2str(obj%r(ia,3)*obj%u%convert_length('from','hau'))//' '//&
      real2str(m(1))//' '//real2str(m(2))//' '//real2str(m(3))
    end do
    call TBKOSTER_flush(unit_rt)
    if(.not. present(unit)) then
      close(unit_rt)
    end if
  end subroutine write_lammps

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_atom.txt')
  subroutine write_txt(obj,file,unit)
    class(atom),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    integer :: ns, na, ntag
    integer,dimension(obj%ntag) :: stag
    character(len=sl),dimension(obj%ntag) :: tag
    integer,dimension(obj%na) :: ia2ie
    integer,dimension(3) :: pbc
    real(rp),dimension(3) :: k_spiral
    character(len=len(obj%r_coord)) :: r_coord
    real(rp),dimension(obj%na,3) :: r, p, m
    character(len=len(obj%m_coord)) :: m_coord
    real(rp),dimension(obj%na)   :: lambda_pen
    ! Namelist
    namelist /atom/ ns, na, ntag, stag, tag, ia2ie, pbc, k_spiral, r_coord, &
     m_coord, r, p, m, lambda_pen

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_atom.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    ns = obj%ns
    na = obj%na
    ntag = obj%ntag
    stag = obj%stag
    tag  = obj%tag
    ia2ie = obj%ia2ie
    pbc = obj%pbc
    r_coord = obj%r_coord
    if(obj%r_coord == 'cartesian') then
      r = obj%r * obj%u%convert_length('from','hau')
    else
      r = obj%r
    end if
    p = obj%p
    m_coord = obj%m_coord
    m = obj%m
    if((obj%ns == 4).and.(obj%m_coord == 'spherical')) then
      m(:,2:3) = m(:,2:3) * rad2deg
    end if

    lambda_pen = obj%lambda_pen * obj%u%convert_energy('from','hau')

    write(unit_rt,nml=atom)

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
  !> 'out_atom.txt'), if tag (default: .true.) the namelist opening and closing
  !> tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit,access)
    class(atom),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),intent(in),optional :: access
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    real(rp),dimension(obj%na,3) :: r
    real(rp),dimension(obj%na,obj%nn_max,3) :: rn
    real(rp),dimension(obj%na) :: lambda_pen
    ! Local variables
    integer :: ia, in, ip, itag

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_atom.txt'
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
      if (present(access)) then
        open(unit=unit_rt,file=file_rt,action='write',access=access)
      else
        open(unit=unit_rt,file=file_rt,action='write')
      end if
    end if
    if(tag_rt) then
      write(unit_rt,'(a)') '&atom'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('ns')
        write(unit_rt,'(a)') ' ns = ' // int2str(obj%ns)
      case('na')
        write(unit_rt,'(a)') ' na = ' // int2str(obj%na)
      case('ntag')
        write(unit_rt,'(a)') ' ntag = ' // int2str(obj%ntag)
      case('stag')
        do itag=1,obj%ntag
          write(unit_rt,'(a)') ' stag(' // int2str(itag) // ') = ' &
           // int2str(obj%stag(itag))
        end do
      case('tag')
        do itag=1,obj%ntag
          write(unit_rt,'(a)') ' tag(' // int2str(itag) // ') = ''' &
           // trim(obj%tag(itag)) // ''''
        end do
      case('ia2ie')
        do ia=1,obj%na
          write(unit_rt,'(a)') ' ia2ie(' // int2str(ia) // ') = ' &
           // int2str(obj%ia2ie(ia))
        end do
      case('nel')
        write(unit_rt,'(a)') ' nel = ' // real2str(obj%nel)
      case('pbc')
        write(unit_rt,'(a)') ' pbc = ' &
         // int2str(obj%pbc(1)) // ', ' &
         // int2str(obj%pbc(2)) // ', ' &
         // int2str(obj%pbc(3))
      case('k_spiral')
        write(unit_rt,'(a)') ' k_spiral = ' &
         // real2str(obj%k_spiral(1)) // ', ' &
         // real2str(obj%k_spiral(2)) // ', ' &
         // real2str(obj%k_spiral(3))
      case('r_coord')
        write(unit_rt,'(a)') ' r_coord = ''' // obj%r_coord // ''''
      case('m_coord')
        write(unit_rt,'(a)') ' m_coord = ''' // obj%m_coord // ''''
      case('r')
        do ia=1,obj%na
          if(obj%r_coord == 'cartesian') then
            r = obj%r * obj%u%convert_length('from','hau')
          else
            r = obj%r
          end if
          write(unit_rt,'(a)') ' r(' // int2str(ia) // ',:) = ' &
           // fixedreal2str(r(ia,1)) // ', ' &
           // fixedreal2str(r(ia,2)) // ', ' &
           // fixedreal2str(r(ia,3))
        end do
      case('p')
        do ia=1,obj%na
          write(unit_rt,'(a)') ' p(' // int2str(ia) // ',:) = ' &
           // real2str(obj%p(ia,1)) // ', ' &
           // real2str(obj%p(ia,2)) // ', ' &
           // real2str(obj%p(ia,3))
        end do
      case('m')
        do ia=1,obj%na
          write(unit_rt,'(a)') ' m(' // int2str(ia) // ',:) = ' &
            // real2str(obj%m(ia,1)) // ', ' &
            // real2str(obj%m(ia,2) * rad2deg) // ', ' &
            // real2str(obj%m(ia,3) * rad2deg)
        end do
      case('nn')
        do ia=1,obj%na
          write(unit_rt,'(a)') ' nn(' // int2str(ia) // ') = ' &
           // int2str(obj%nn(ia))
        end do
      case('nn_max')
        write(unit_rt,'(a)') ' nn_max = ' // int2str(obj%nn_max)
      case('ian2ia')
        do ia=1,obj%na
          do in=1,obj%nn(ia)
            write(unit_rt,'(a)') ' ian2ia(' // int2str(ia) // ',' &
             // int2str(in) // ') = ' // int2str(obj%ian2ia(ia,in))
          end do
        end do
      case('rn')
        rn = obj%rn * obj%u%convert_length('from','hau')
        do ia=1,obj%na
          do in=1,obj%nn(ia)
            write(unit_rt,'(a)') ' rn(' // int2str(ia) // ',' // int2str(in) &
             // ',:) = ' &
             // real2str(rn(ia,in,1)) // ', ' &
             // real2str(rn(ia,in,2)) // ', ' &
             // real2str(rn(ia,in,3))
          end do
        end do
      case('lambda_pen')
        lambda_pen = obj%lambda_pen * obj%u%convert_energy('from','hau')
        do ia=1,obj%na
          write(unit_rt,'(a)') ' lambda_pen(' // int2str(ia) // ') = ' &
           // real2str(lambda_pen(ia))
        end do
      case('b_pen')
        do ia=1,obj%na
          write(unit_rt,'(a)') ' b_pen(' // int2str(ia) // ',:) = ' &
           // real2str(obj%b_pen(ia,1)) // ', ' &
           // real2str(obj%b_pen(ia,2)) // ', ' &
           // real2str(obj%b_pen(ia,3))
        end do
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

  !> Write object in extended .xyz format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_atom.xyz')
  subroutine write_xyz(obj,file,unit)
    class(atom),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    integer :: ia
    real(rp),dimension(3) :: m_cart,local_m

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_atom.xyz'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    ! Number of atoms
    write(unit_rt,'(a)') int2str(obj%na)
    ! Serie of key/value pairs
    write(unit_rt,'(a)') 'Lattice="' &
     // real2str(obj%l_r%v(1,1)* obj%u%convert_length('from','hau')) // ' ' &
     // real2str(obj%l_r%v(1,2)* obj%u%convert_length('from','hau')) // ' ' &
     // real2str(obj%l_r%v(1,3)* obj%u%convert_length('from','hau')) // ' ' &
     // real2str(obj%l_r%v(2,1)* obj%u%convert_length('from','hau')) // ' ' &
     // real2str(obj%l_r%v(2,2)* obj%u%convert_length('from','hau')) // ' ' &
     // real2str(obj%l_r%v(2,3)* obj%u%convert_length('from','hau')) // ' ' &
     // real2str(obj%l_r%v(3,1)* obj%u%convert_length('from','hau')) // ' ' &
     // real2str(obj%l_r%v(3,2)* obj%u%convert_length('from','hau')) // ' ' &
     // real2str(obj%l_r%v(3,3)* obj%u%convert_length('from','hau')) // '" ' &
     // 'Properties=species:S:1:pos:R:3:dipoles:R:3'
    ! Atoms
    do ia=1,obj%na
      ! temporary array to avoid a compiler warning "An array temporary was created"
      local_m = obj%m(ia,:)
      m_cart = sph2cart(local_m)
      write(unit_rt,'(a)') obj%e%symbol(obj%ia2ie(ia)) // ' ' &
       // real2str(obj%r(ia,1)* obj%u%convert_length('from','hau')) // ' ' &
       // real2str(obj%r(ia,2)* obj%u%convert_length('from','hau')) // ' ' &
       // real2str(obj%r(ia,3)* obj%u%convert_length('from','hau')) // ' ' &
       // real2str(m_cart(1)) // ' ' &
       // real2str(m_cart(2)) // ' ' &
       // real2str(m_cart(3))
    end do

    if(.not. present(unit)) then
      call TBKOSTER_flush(unit_rt)
      close(unit_rt)
    else
      call TBKOSTER_flush(unit)
    end if
    !deallocate(file_rt)
  end subroutine write_xyz
end module atom_mod
