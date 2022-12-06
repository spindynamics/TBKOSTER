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
!  element.f90
!  DyNaMol
module element_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use math_mod, only: i_unit, sqrt_three
  use precision_mod, only: rp
  use string_mod, only: sl, int2str, lower, real2str, dynamol_flush
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'ne', &
   'symbol', &
   'no', &
   'o', &
   'q', &
   'q_s', &
   'q_p', &
   'q_d', &
   'u_lcn', &
   'u_lcn_d', &
   'i_stoner_d', &
   'b', &
   'j_dd', &
   'u_dd', &
   'xi_so_p', &
   'xi_so_d', &
   'mass' &
   ]

  type,public :: element
    !> @defgroup Element Element-related variables
    !> @{

    !> Number of elements
    integer :: ne
    !> Element symbols
    character(len=2),dimension(:),allocatable :: symbol
    !> Element names
    character(len=sl),dimension(:),allocatable :: name
    !> Element numbers
    integer,dimension(:),allocatable :: number
    !> @}

    !> @defgroup Orbital Orbital-related variables
    !> Orbitals are combinations of the azimuthal and magnetic quantum numbers
    !> \f$ (lm) \f$, often grouped under the symbol \f$ L \f$
    ! Format 1
    ! \f$ L=1: s           \; (l=0,m=   0) \f$\n
    ! \f$ L=2: p_x         \; (l=1,m=\pm1) \f$\n
    ! \f$ L=3: p_y         \; (l=1,m=\pm1) \f$\n
    ! \f$ L=4: p_z         \; (l=1,m=   0) \f$\n
    ! \f$ L=5: d_{xy}      \; (l=2,m=\pm2) \f$\n
    ! \f$ L=6: d_{yz}      \; (l=2,m=\pm1) \f$\n
    ! \f$ L=7: d_{xz}      \; (l=2,m=\pm1) \f$\n
    ! \f$ L=8: d_{x^2-y^2} \; (l=2,m=\pm2) \f$\n
    ! \f$ L=9: d_{z^2}     \; (l=2,m=   0) \f$\n
    ! Format 2
    !> \f{eqnarray*}{
    !>   L=1: &s           \; (l=0,m=   0)\\
    !>   L=2: &p_x         \; (l=1,m=\pm1)\\
    !>   L=3: &p_y         \; (l=1,m=\pm1)\\
    !>   L=4: &p_z         \; (l=1,m=   0)\\
    !>   L=5: &d_{xy}      \; (l=2,m=\pm2)\\
    !>   L=6: &d_{yz}      \; (l=2,m=\pm1)\\
    !>   L=7: &d_{xz}      \; (l=2,m=\pm1)\\
    !>   L=8: &d_{x^2-y^2} \; (l=2,m=\pm2)\\
    !>   L=9: &d_{z^2}     \; (l=2,m=   0)
    !> \f}
    !> @{

    !> Number of valence orbitals
    integer,dimension(:),allocatable :: no
    !> Maximum number of valence orbitals
    integer :: no_max
    !> Valence orbitals
    integer,dimension(:,:),allocatable :: o
    !> Orbital to azimuthal quantum number + 1
    integer,dimension(9) :: o2l
    !> Valence orbital symbols
    character(len=3),dimension(:),allocatable :: os
    !> @}

    !> @defgroup Charge Charge-related variables
    !> @{

    !> Valence charges \f$ q \f$ (total and orbital-projected)
    real(rp),dimension(:),allocatable :: q
    real(rp),dimension(:),allocatable :: q_s
    real(rp),dimension(:),allocatable :: q_p
    real(rp),dimension(:),allocatable :: q_d
    !> Local charge neutrality \f$ U \f$ (total and for d-orbitals)
    real(rp),dimension(:),allocatable :: u_lcn
    real(rp),dimension(:),allocatable :: u_lcn_d
    !> @}

    !> @defgroup Electron-electron_interaction Electron-electron interaction
    !> (only for d-orbitals)
    !> Stoner interaction or intra-atomic Coulomb interaction (Hartree-Fock
    !> method)\n
    !> See \cite Barreteau2016 page 20-21.
    !> @{

    !> Stoner parameters for d orbitals \f$ I^{\mathrm{Stoner}} \f$
    real(rp),dimension(:),allocatable :: i_stoner_d
    !> Racah parameters \f$ B \f$
    real(rp),dimension(:),allocatable :: b
    !> Exchange integrals \f$ J_{dd} \f$
    real(rp),dimension(:),allocatable :: j_dd
    !> Coulomb integrals \f$ U_{dd} \f$
    real(rp),dimension(:),allocatable :: u_dd
    !> On-site Coulomb potential \f$ V_{i\lambda,j\mu}^R \f$
    real(rp), dimension(:,:,:,:,:), allocatable :: v_osc
    !> @}

    !> @defgroup Spin-orbit Spin-orbit-related variables
    !> @{

    !> Spin-orbit constant \f$ \xi^{\mathrm{SO}} \f$ (for p-orbitals and
    !> d-orbitals)
    !> See \cite Barreteau2016 page 16.\n
    !> A spherical approximation to the spin-orbit (SO) potential:\n
    !> \f$ \hat{V}^{\mathrm{SO}} = \xi(r) \frac{\hat{\mathbf{L}}}{\hbar} \cdot
    !> \frac{\hat{\mathbf{S}}}{\hbar}\f$\n
    !> where the orbital moment operator is:\n
    !> \f$ \hat{\mathbf{L}} = \hbar \hat{\mathbf{l}} = \hat{\mathbf{r}} \times
    !> \hat{\mathbf{p}} \f$\n
    !> the spin moment operator is:\n
    !> \f$ \hat{\mathbf{S}} = \hbar \frac{\hat{\sigma}}{2} \f$\n
    !> and the function \f$ \xi(r) \f$ is in terms of the electrostatic
    !> potential \f$ V(r) \f$:\n
    !> \f$  \xi(r) = \frac{\hbar^2}{2m^2c^2} \frac{1}{r} \frac{\mathrm{d}V(r)}
    !> {\mathrm{d}r} \f$\n
    !> The matrix elements of \f$ \hat{V}^{\mathrm{SO}} \f$ in a basis of atomic
    !> spin orbitals \f$ |i\lambda\sigma\rangle \f$ take the form:\n
    !> \f$ \langle i\lambda\sigma | \hat{V}^{\mathrm{SO}} | j\mu\sigma' \rangle
    !> = \xi_{i\lambda\mu} \frac{1}{2} \langle \lambda\sigma | \hat{\mathbf{l}}
    !> \cdot \hat{\mathbf{\sigma}} | \mu\sigma' \rangle \delta_{i,j} \f$\n
    !> The SO potential only has site-diagonal elements since \f$ \xi_i(r) \f$
    !> is localized near \f$ r = 0 \f$. In addition, the angular moment operator
    !> does not couple orbitals of different nature (\f$ \lambda \neq \mu \f$)
    !> and is zero for s orbitals (the angular quantum number \f$ l = 0 \f$ by
    !> definition for these). Only \f$ \xi_p \f$ and \f$ \xi_d \f$ are thus
    !> relevant:\n
    !> \f$ \xi_{i\lambda} = \int_0^{\infty} R_{i\lambda}^2(r) r^2 \xi_i(r)
    !> \mathrm{d}r, \quad \lambda = p,d \f$
    real(rp),dimension(:),allocatable :: xi_so_p
    real(rp),dimension(:),allocatable :: xi_so_d
    complex(rp),dimension(:,:,:,:),allocatable :: l_dot_s

    !> atomic mass
    real(rp), dimension(:),allocatable :: mass

    !> @}

    !> Units
    class(units),pointer :: u

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: calculate_l_dot_s
    procedure :: calculate_v_osc
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type element

  ! Constructor
  interface element
    procedure :: constructor
  end interface element

contains
  function constructor(u) result(obj)
    class(units),target,intent(in) :: u
    type(element) :: obj

    obj%u => u

    obj%o2l = (/1, 2, 2, 2, 3, 3, 3, 3, 3/)
  end function constructor

  subroutine destructor(obj)
    type(element) :: obj

    if(allocated(obj%symbol))     deallocate(obj%symbol)
    if(allocated(obj%name))       deallocate(obj%name)
    if(allocated(obj%number))     deallocate(obj%number)
    if(allocated(obj%no))         deallocate(obj%no)
    if(allocated(obj%o))          deallocate(obj%o)
    if(allocated(obj%os))         deallocate(obj%os)
    if(allocated(obj%q))          deallocate(obj%q)
    if(allocated(obj%q_s))        deallocate(obj%q_s)
    if(allocated(obj%q_p))        deallocate(obj%q_p)
    if(allocated(obj%q_d))        deallocate(obj%q_d)
    if(allocated(obj%u_lcn))      deallocate(obj%u_lcn)
    if(allocated(obj%u_lcn_d))    deallocate(obj%u_lcn_d)
    if(allocated(obj%i_stoner_d)) deallocate(obj%i_stoner_d)
    if(allocated(obj%b))          deallocate(obj%b)
    if(allocated(obj%j_dd))       deallocate(obj%j_dd)
    if(allocated(obj%u_dd))       deallocate(obj%u_dd)
    if(allocated(obj%v_osc))      deallocate(obj%v_osc)
    if(allocated(obj%xi_so_p))    deallocate(obj%xi_so_p)
    if(allocated(obj%xi_so_d))    deallocate(obj%xi_so_d)
    if(allocated(obj%l_dot_s))    deallocate(obj%l_dot_s)
    if(allocated(obj%mass))       deallocate(obj%mass)
  end subroutine destructor

  !> Build the on-site Coulomb potential
  function build_v_osc(b,j_dd,u_dd) result(v_osc)
    real(rp),intent(in) :: b, j_dd, u_dd
    real(rp),dimension(5,5,5,5) :: v_osc
    integer :: i,j

    v_osc = 0.0_rp

    ! 1 <-> dxy
    ! 2 <-> dyz
    ! 3 <-> dzx
    ! 4 <-> dx^2-y^2
    ! 5 <-> d(3z^2-r^2)

    !
    ! one orbital terms
    !
    do i=1,5
      v_osc(i,i,i,i) = u_dd+2*j_dd
    end do

    !
    ! two orbitals terms
    !
    ! Uijij:  Coulomb
    v_osc(1,2,1,2) = u_dd-b
    v_osc(1,3,1,3) = u_dd-b
    v_osc(2,3,2,3) = u_dd-b
    v_osc(2,4,2,4) = u_dd-b
    v_osc(3,4,3,4) = u_dd-b
    v_osc(1,4,1,4) = u_dd+5*b
    v_osc(1,5,1,5) = u_dd-3*b
    v_osc(4,5,4,5) = u_dd-3*b
    v_osc(2,5,2,5) = u_dd+3*b
    v_osc(3,5,3,5) = u_dd+3*b
    do i=2,5
      do j=1,i-1
        v_osc(i,j,i,j) = v_osc(j,i,j,i)
      end do
    end do
    ! Uijji:  exchange integral
    call symmetrize_v_osc(1,2,2,1,j_dd+b/2,  v_osc)
    call symmetrize_v_osc(1,3,3,1,j_dd+b/2,  v_osc)
    call symmetrize_v_osc(2,3,3,2,j_dd+b/2,  v_osc)
    call symmetrize_v_osc(2,4,4,2,j_dd+b/2,  v_osc)
    call symmetrize_v_osc(3,4,4,3,j_dd+b/2,  v_osc)
    call symmetrize_v_osc(1,4,4,1,j_dd-5*b/2,v_osc)
    call symmetrize_v_osc(1,5,5,1,j_dd+3*b/2,v_osc)
    call symmetrize_v_osc(4,5,5,4,j_dd+3*b/2,v_osc)
    call symmetrize_v_osc(2,5,5,2,j_dd-3*b/2,v_osc)
    call symmetrize_v_osc(3,5,5,3,j_dd-3*b/2,v_osc)

    !
    ! three orbitals terms
    !
    call symmetrize_v_osc(4,2,2,5,  -b*sqrt(3.0_rp),v_osc)
    call symmetrize_v_osc(4,3,3,5,   b*sqrt(3.0_rp),v_osc)
    call symmetrize_v_osc(4,2,5,2, 2*b*sqrt(3.0_rp),v_osc)
    call symmetrize_v_osc(4,3,5,3,-2*b*sqrt(3.0_rp),v_osc)

    !
    ! four orbitals terms
    !
    call symmetrize_v_osc(1,2,5,3,-2*b*sqrt(3.0_rp),v_osc)
    call symmetrize_v_osc(1,2,3,5,   b*sqrt(3.0_rp),v_osc)
    call symmetrize_v_osc(1,4,3,2,-  b*3,           v_osc)
    call symmetrize_v_osc(1,3,2,5,   b*sqrt(3.0_rp),v_osc)
    call symmetrize_v_osc(1,4,2,3,   b*3,           v_osc)
  end function build_v_osc

  subroutine calculate_l_dot_s(obj)
    class(element),intent(inout) :: obj
    integer :: ie,io1,io2,is1,is2,ils1,ils2,ils3,ils4

    if(allocated(obj%l_dot_s)) deallocate(obj%l_dot_s)
    allocate(obj%l_dot_s(obj%ne,3,obj%no_max*2,obj%no_max*2))
    obj%l_dot_s = cmplx(0.0_rp,0.0_rp,kind=rp)

    do ie=1,obj%ne
      is1=1
      is2=1
      do io1=1,obj%no(ie)
        do io2=1,obj%no(ie)
          ils1=is1+(io1-1)*2
          ils2=is2+(io2-1)*2

          if(obj%o(ie,io1)==2 .and. obj%o(ie,io2)==3) then     !p bands
            obj%l_dot_s(ie,2,ils1,ils2) = -i_unit
          elseif(obj%o(ie,io1)==3 .and. obj%o(ie,io2)==2) then !p bands
            obj%l_dot_s(ie,2,ils1,ils2) = i_unit
          elseif(obj%o(ie,io1)==5 .and. obj%o(ie,io2)==8) then !d bands
            obj%l_dot_s(ie,3,ils1,ils2) = 2*i_unit
          elseif(obj%o(ie,io1)==6 .and. obj%o(ie,io2)==7) then !d bands
            obj%l_dot_s(ie,3,ils1,ils2) = i_unit
          elseif(obj%o(ie,io1)==8 .and. obj%o(ie,io2)==5) then !d bands
            obj%l_dot_s(ie,3,ils1,ils2) = -2*i_unit
          elseif(obj%o(ie,io1)==7 .and. obj%o(ie,io2)==6) then !d bands
            obj%l_dot_s(ie,3,ils1,ils2) = -i_unit
          end if
        end do
      end do

      is1=1
      is2=2
      do io1=1,obj%no(ie)
        do io2=1,obj%no(ie)
          ils1=is1+(io1-1)*2
          ils2=is2+(io2-1)*2

          if(obj%o(ie,io1)==2 .and. obj%o(ie,io2)==4) then
            obj%l_dot_s(ie,2,ils1,ils2) = cmplx(1.0_rp,0.0_rp,kind=rp)
          elseif(obj%o(ie,io1)==3 .and. obj%o(ie,io2)==4) then
            obj%l_dot_s(ie,2,ils1,ils2) = -i_unit
          elseif(obj%o(ie,io1)==4 .and. obj%o(ie,io2)==2) then
            obj%l_dot_s(ie,2,ils1,ils2) = cmplx(-1.0_rp,0.0_rp,kind=rp)
          elseif(obj%o(ie,io1)==4 .and. obj%o(ie,io2)==3) then
            obj%l_dot_s(ie,2,ils1,ils2) = i_unit
          elseif(obj%o(ie,io1)==5 .and. obj%o(ie,io2)==6) then
            obj%l_dot_s(ie,3,ils1,ils2) = cmplx(1.0_rp,0.0_rp,kind=rp)
          elseif(obj%o(ie,io1)==5 .and. obj%o(ie,io2)==7) then
            obj%l_dot_s(ie,3,ils1,ils2) = -i_unit
          elseif(obj%o(ie,io1)==6 .and. obj%o(ie,io2)==8) then
            obj%l_dot_s(ie,3,ils1,ils2) = -i_unit
          elseif(obj%o(ie,io1)==6 .and. obj%o(ie,io2)==9) then
            obj%l_dot_s(ie,3,ils1,ils2) = -sqrt_three*i_unit
          elseif(obj%o(ie,io1)==7 .and. obj%o(ie,io2)==8) then
            obj%l_dot_s(ie,3,ils1,ils2) = cmplx(-1.0_rp,0.0_rp,kind=rp)
          elseif(obj%o(ie,io1)==7 .and. obj%o(ie,io2)==9) then
            obj%l_dot_s(ie,3,ils1,ils2) = sqrt_three
          elseif(obj%o(ie,io1)==6 .and. obj%o(ie,io2)==5) then
            obj%l_dot_s(ie,3,ils1,ils2) = cmplx(-1.0_rp,0.0_rp,kind=rp)
          elseif(obj%o(ie,io1)==7 .and. obj%o(ie,io2)==5) then
            obj%l_dot_s(ie,3,ils1,ils2) = i_unit
          elseif(obj%o(ie,io1)==8 .and. obj%o(ie,io2)==6) then
            obj%l_dot_s(ie,3,ils1,ils2) = i_unit
          elseif(obj%o(ie,io1)==9 .and. obj%o(ie,io2)==6) then
            obj%l_dot_s(ie,3,ils1,ils2) = sqrt_three*i_unit
          elseif(obj%o(ie,io1)==8 .and. obj%o(ie,io2)==7) then
            obj%l_dot_s(ie,3,ils1,ils2) = cmplx(1.0_rp,0.0_rp,kind=rp)
          elseif(obj%o(ie,io1)==9 .and. obj%o(ie,io2)==7) then
            obj%l_dot_s(ie,3,ils1,ils2) = -sqrt_three
          end if
        end do
      end do

      is1=2
      is2=2
      do io1=1,obj%no(ie)
        do io2=1,obj%no(ie)
          ils1=is1+(io1-1)*2
          ils2=is2+(io2-1)*2
          ils3=1+(io1-1)*2
          ils4=1+(io2-1)*2
          obj%l_dot_s(ie,2,ils1,ils2) = conjg(obj%l_dot_s(ie,2,ils3,ils4))
          obj%l_dot_s(ie,3,ils1,ils2) = conjg(obj%l_dot_s(ie,3,ils3,ils4))
        end do
      end do

      is1=2
      is2=1
      do io1=1,obj%no(ie)
        do io2=1,obj%no(ie)
          ils1=is1+(io1-1)*2
          ils2=is2+(io2-1)*2
          ils3=1+(io1-1)*2
          ils4=2+(io2-1)*2
          obj%l_dot_s(ie,2,ils1,ils2) = -conjg(obj%l_dot_s(ie,2,ils3,ils4))
          obj%l_dot_s(ie,3,ils1,ils2) = -conjg(obj%l_dot_s(ie,3,ils3,ils4))
        end do
      end do
      obj%l_dot_s(ie,2,:,:) = obj%xi_so_p(ie)/2*obj%l_dot_s(ie,2,:,:)
      obj%l_dot_s(ie,3,:,:) = obj%xi_so_d(ie)/2*obj%l_dot_s(ie,3,:,:)
    end do ! end loop on ie
  end subroutine calculate_l_dot_s

  subroutine calculate_v_osc(obj)
    class(element),intent(inout) :: obj
    integer :: ie

    if(allocated(obj%v_osc)) deallocate(obj%v_osc)
    allocate(obj%v_osc(obj%ne,5,5,5,5))
    obj%v_osc = 0.0_rp

    do ie=1,obj%ne
      obj%v_osc(ie,:,:,:,:) = build_v_osc(obj%b(ie),obj%j_dd(ie),obj%u_dd(ie))
    end do
  end subroutine calculate_v_osc

  ! Initialize electron-electron interaction
  subroutine initialize_eei(ne, i_stoner_d, b, j_dd, u_dd)
    integer,intent(in) :: ne
    real(rp),dimension(:),allocatable,intent(out) :: i_stoner_d, b, j_dd, u_dd

    allocate(i_stoner_d(ne))
    allocate(b(ne))
    allocate(j_dd(ne))
    allocate(u_dd(ne))

    i_stoner_d = 0.0_rp
    b    = 0.0_rp
    j_dd = 0.0_rp
    u_dd = 0.0_rp
  end subroutine initialize_eei

  ! Initialize local charge neutrality
  subroutine initialize_lcn(ne, u_lcn, u_lcn_d)
    integer,intent(in) :: ne
    real(rp),dimension(:),allocatable,intent(out) :: u_lcn, u_lcn_d

    allocate(u_lcn(ne))
    allocate(u_lcn_d(ne))

    u_lcn   = 0.0_rp
    u_lcn_d = 0.0_rp
  end subroutine initialize_lcn

  ! Initialize orbitals
  subroutine initialize_o(ne, symbol, no_max, o, os)
    integer,intent(in) :: ne
    character(len=*),dimension(:),intent(in) :: symbol
    character(len=1) :: first_c
    integer,intent(in) :: no_max
    integer,dimension(:,:),allocatable,intent(out) :: o
    character(len=3),dimension(:),allocatable,intent(out) :: os
    integer :: ie

    allocate(o(ne,no_max))
    allocate(os(ne))
    o = 0.0_rp
    do ie=1,ne
      select case(lower(trim(symbol(ie))))
      ! 1st line
      case('h')
        o(ie,1) = 1
      case('he')
        o(ie,1) = 1
      ! 2nd line
      case('li')
        o(ie,1) = 1
      case('be')
        o(ie,1) = 1
      case('b')
        o(ie,1:4) = (/1,2,3,4/)
      case('c')
        o(ie,1:4) = (/1,2,3,4/)
      case('n')
        o(ie,1:4) = (/1,2,3,4/)
      case('o')
        o(ie,1:4) = (/1,2,3,4/)
      case('f')
        o(ie,1:4) = (/1,2,3,4/)
      case('ne')
        o(ie,1:4) = (/1,2,3,4/)
      ! 3rd line
      case('na')
        o(ie,1) = 1
      case('mg')
        o(ie,1) = 1
      case('al')
        o(ie,1:4) = (/1,2,3,4/)
      case('si')
        o(ie,1:4) = (/1,2,3,4/)
      case('p')
        o(ie,1:4) = (/1,2,3,4/)
      case('s')
        o(ie,1:4) = (/1,2,3,4/)
      case('cl')
        o(ie,1:4) = (/1,2,3,4/)
      case('ar')
        o(ie,1:4) = (/1,2,3,4/)
      ! 4th line
      case('k')
        o(ie,1) = 1
      case('ca')
        o(ie,1) = 1
      case('sc')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('ti')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('v')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('cr')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('mn')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('fe')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('co')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('ni')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('cu')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('zn')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('ga')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('ge')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('as')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('se')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('br')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('kr')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      ! 5th line
      case('rb')
        o(ie,1) = 1
      case('sr')
        o(ie,1) = 1
      case('y')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('zr')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('nb')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('mo')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('tc')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('ru')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('rh')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('pd')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('ag')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('cd')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('in')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('sn')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('sb')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('te')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('i')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('xe')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      ! 6th line
      case('cs')
        o(ie,1) = 1
      case('ba')
        o(ie,1) = 1
      case('la')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      ! [Lanthanides]
      case('hf')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('ta')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('w')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('re')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('os')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('ir')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('pt')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('au')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('hg')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('tl')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('pb')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('bi')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('po')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('at')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case('rn')
        o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
      case default
        first_c=lower(trim(symbol(ie)))
        if(first_c=='j') then
       ! TB model
          o(ie,1:9) = (/1,2,3,4,5,6,7,8,9/)
        else
        write(error_unit,*) 'element%parse_symbol(): element symbol ', &
         symbol(ie), ' not known'
        error stop
        endif
      end select
      os(ie) = parse_orbital_to_string(o(ie,:))
    end do
  end subroutine initialize_o

  ! Initialize charge-related variables
  subroutine initialize_q(ne, symbol, name, number, no, no_max, q, q_s, q_p, &
   q_d, mass)
    integer,intent(in) :: ne
    character(len=*),dimension(:),intent(in) :: symbol
    character(len=1) :: first_c
    character(len=sl),dimension(:),allocatable,intent(out) :: name
    integer,dimension(:),allocatable,intent(out) :: number
    integer,dimension(:),allocatable,intent(out) :: no
    integer,intent(out) :: no_max
    real(rp),dimension(:),allocatable,intent(out) :: q, q_s, q_p, q_d
    real(rp),dimension(:),allocatable,intent(out) :: mass
    integer :: ie

    allocate(name(ne))
    allocate(number(ne))
    allocate(no(ne))
    allocate(q(ne))
    allocate(q_s(ne))
    allocate(q_p(ne))
    allocate(q_d(ne))
    allocate(mass(ne))
    do ie=1,ne
      ! Periodic table
      select case(lower(trim(symbol(ie))))
      ! 1st line
      case('h')
        name(ie) = 'hydrogen'
        number(ie) = 1
        no(ie) = 1
        q(ie)   = 1.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 1.0079_rp
      case('he')
        name(ie) = 'helium'
        number(ie) = 2
        no(ie) = 1
        q(ie)   = 2.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 4.0026_rp
      ! 2nd line
      case('li')
        name(ie) = 'lithium'
        number(ie) = 3
        no(ie) = 1
        q(ie)   = 1.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 6.941_rp
      case('be')
        name(ie) = 'beryllium'
        number(ie) = 4
        no(ie) = 1
        q(ie)   = 2.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 9.0122_rp
      case('b')
        name(ie) = 'boron'
        number(ie) = 5
        no(ie) = 4
        q(ie)   = 3.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 1.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 10.811_rp
      case('c')
        name(ie) = 'carbon'
        number(ie) = 6
        no(ie) = 4
        q(ie)   = 4.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 2.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 12.0107_rp
      case('n')
        name(ie) = 'nitrogen'
        number(ie) = 7
        no(ie) = 4
        q(ie)   = 5.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 3.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 14.0067_rp
      case('o')
        name(ie) = 'oxygen'
        number(ie) = 8
        no(ie) = 4
        q(ie)   = 6.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 4.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 15.9994_rp
      case('f')
        name(ie) = 'fluorine'
        number(ie) = 9
        no(ie) = 4
        q(ie)   = 7.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 5.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 18.9984_rp
      case('ne')
        name(ie) = 'neon'
        number(ie) = 10
        no(ie) = 4
        q(ie)   = 8.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 6.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 20.1797_rp
      ! 3rd line
      case('na')
        name(ie) = 'sodium'
        number(ie) = 11
        no(ie) = 1
        q(ie)   = 1.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 22.9897_rp
      case('mg')
        name(ie) = 'magnesium'
        number(ie) = 12
        no(ie) = 1
        q(ie)   = 2.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 24.305_rp
      case('al')
        name(ie) = 'aluminium'
        number(ie) = 13
        no(ie) = 4
        q(ie)   = 3.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 1.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 26.9815_rp
      case('si')
        name(ie) = 'silicon'
        number(ie) = 14
        no(ie) = 4
        q(ie)   = 4.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 2.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 28.0855_rp
      case('p')
        name(ie) = 'phosphorus'
        number(ie) = 15
        no(ie) = 4
        q(ie)   = 5.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 3.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 30.9738_rp
      case('s')
        name(ie) = 'sulfur'
        number(ie) = 16
        no(ie) = 4
        q(ie)   = 6.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 4.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 32.065_rp
      case('cl')
        name(ie) = 'chlorine'
        number(ie) = 17
        no(ie) = 4
        q(ie)   = 7.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 5.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 35.453_rp
      case('ar')
        name(ie) = 'argon'
        number(ie) = 18
        no(ie) = 4
        q(ie)   = 8.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 6.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 39.948
      ! 4th line
      case('k')
        name(ie) = 'potassium'
        number(ie) = 19
        no(ie) = 1
        q(ie)   = 1.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 39.0983_rp
      case('ca')
        name(ie) = 'calcium'
        number(ie) = 20
        no(ie) = 1
        q(ie)   = 2.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 40.078_rp
      case('sc')
        name(ie) = 'scandium'
        number(ie) = 21
        no(ie) = 9
        q(ie)   = 3.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 1.0_rp
        mass(ie)= 44.9559_rp
      case('ti')
        name(ie) = 'titanium'
        number(ie) = 22
        no(ie) = 9
        q(ie)   = 4.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 2.0_rp
        mass(ie)= 47.867_rp
      case('v')
        name(ie) = 'vanadium'
        number(ie) = 23
        no(ie) = 9
        q(ie)   = 5.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 3.0_rp
        mass(ie)= 50.9415_rp
      case('cr')
        name(ie) = 'chromium'
        number(ie) = 24
        no(ie) = 9
        q(ie)   = 6.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 5.0_rp
        mass(ie)= 51.9961_rp
      case('mn')
        name(ie) = 'manganese'
        number(ie) = 25
        no(ie) = 9
        q(ie)   = 7.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 5.0_rp
        mass(ie)= 54.938_rp
      case('fe')
        name(ie) = 'iron'
        number(ie) = 26
        no(ie) = 9
        q(ie)   = 8.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 6.0_rp
        mass(ie)= 55.845_rp
      case('co')
        name(ie) = 'cobalt'
        number(ie) = 27
        no(ie) = 9
        q(ie)   = 9.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 7.0_rp
        mass(ie)= 58.9332_rp
      case('ni')
        name(ie) = 'nickel'
        number(ie) = 28
        no(ie) = 9
        q(ie)   = 10.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 8.0_rp
        mass(ie)= 58.6934_rp
      case('cu')
        name(ie) = 'copper'
        number(ie) = 29
        no(ie) = 9
        q(ie)   = 11.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 63.546_rp
      case('zn')
        name(ie) = 'zinc'
        number(ie) = 30
        no(ie) = 9
        q(ie)   = 12.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 65.39_rp
      case('ga')
        name(ie) = 'gallium'
        number(ie) = 31
        no(ie) = 9
        q(ie)   = 13.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 1.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 69.723_rp
      case('ge')
        name(ie) = 'germanium'
        number(ie) = 32
        no(ie) = 9
        q(ie)   = 14.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 2.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 72.64_rp
      case('as')
        name(ie) = 'arsenic'
        number(ie) = 33
        no(ie) = 9
        q(ie)   = 15.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 3.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 74.9216_rp
      case('se')
        name(ie) = 'selenium'
        number(ie) = 34
        no(ie) = 9
        q(ie)   = 16.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 4.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 78.96_rp
      case('br')
        name(ie) = 'bromine'
        number(ie) = 35
        no(ie) = 9
        q(ie)   = 17.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 5.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 79.904_rp
      case('kr')
        name(ie) = 'krypton'
        number(ie) = 36
        no(ie) = 9
        q(ie)   = 18.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 6.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 83.8_rp
      ! 5th line
      case('rb')
        name(ie) = 'rubidium'
        number(ie) = 37
        no(ie) = 1
        q(ie)   = 1.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 85.4678_rp
      case('sr')
        name(ie) = 'strontium'
        number(ie) = 38
        no(ie) = 1
        q(ie)   = 2.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 87.62_rp
      case('y')
        name(ie) = 'yttrium'
        number(ie) = 39
        no(ie) = 9
        q(ie)   = 3.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 1.0_rp
        mass(ie)= 88.9059_rp
      case('zr')
        name(ie) = 'zirconium'
        number(ie) = 40
        no(ie) = 9
        q(ie)   = 4.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 2.0_rp
        mass(ie)= 91.224_rp
      case('nb')
        name(ie) = 'niobium'
        number(ie) = 41
        no(ie) = 9
        q(ie)   = 5.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 4.0_rp
        mass(ie)= 92.9064_rp
      case('mo')
        name(ie) = 'molybdenum'
        number(ie) = 42
        no(ie) = 9
        q(ie)   = 6.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 5.0_rp
        mass(ie)= 95.94_rp
      case('tc')
        name(ie) = 'technetium'
        number(ie) = 43
        no(ie) = 9
        q(ie)   = 7.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 5.0_rp
        mass(ie)= 98_rp
      case('ru')
        name(ie) = 'ruthenium'
        number(ie) = 44
        no(ie) = 9
        q(ie)   = 8.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 7.0_rp
        mass(ie)= 101.07_rp
      case('rh')
        name(ie) = 'rhodium'
        number(ie) = 45
        no(ie) = 9
        q(ie)   = 9.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 8.0_rp
        mass(ie)= 102.9055_rp
      case('pd')
        name(ie) = 'palladium'
        number(ie) = 46
        no(ie) = 9
        q(ie)   = 10.0_rp
        q_s(ie) = 0.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 106.42_rp
      case('ag')
        name(ie) = 'silver'
        number(ie) = 47
        no(ie) = 9
        q(ie)   = 11.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 107.8682_rp
      case('cd')
        name(ie) = 'cadmium'
        number(ie) = 48
        no(ie) = 9
        q(ie)   = 12.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 112.411_rp
      case('in')
        name(ie) = 'indium'
        number(ie) = 49
        no(ie) = 9
        q(ie)   = 13.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 1.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 114.818_rp
      case('sn')
        name(ie) = 'tin'
        number(ie) = 50
        no(ie) = 9
        q(ie)   = 14.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 2.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 118.71_rp
      case('sb')
        name(ie) = 'antimony'
        number(ie) = 51
        no(ie) = 9
        q(ie)   = 15.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 3.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 121.76_rp
      case('te')
        name(ie) = 'tellurium'
        number(ie) = 52
        no(ie) = 9
        q(ie)   = 16.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 4.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 127.6_rp
      case('i')
        name(ie) = 'iodine'
        number(ie) = 53
        no(ie) = 9
        q(ie)   = 17.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 5.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 126.9045_rp
      case('xe')
        name(ie) = 'xenon'
        number(ie) = 54
        no(ie) = 9
        q(ie)   = 18.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 6.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 131.293_rp
      ! 6th line
      case('cs')
        name(ie) = 'caesium'
        number(ie) = 55
        no(ie) = 1
        q(ie)   = 1.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 132.9055_rp
      case('ba')
        name(ie) = 'barium'
        number(ie) = 56
        no(ie) = 1
        q(ie)   = 2.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 137.327_rp
      case('la')
        name(ie) = 'lanthanum'
        number(ie) = 57
        no(ie) = 9
        q(ie)   = 3.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 1.0_rp
        mass(ie)= 138.9055_rp
      ! [Lanthanides]
      case('hf')
        name(ie) = 'hafnium'
        number(ie) = 72
        no(ie) = 9
        q(ie)   = 4.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 2.0_rp
        mass(ie)= 178.49_rp
      case('ta')
        name(ie) = 'tantalum'
        number(ie) = 73
        no(ie) = 9
        q(ie)   = 5.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 3.0_rp
        mass(ie)= 180.9479_rp
      case('w')
        name(ie) = 'tungsten'
        number(ie) = 74
        no(ie) = 9
        q(ie)   = 6.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 4.0_rp
        mass(ie)= 183.84_rp
      case('re')
        name(ie) = 'rhenium'
        number(ie) = 75
        no(ie) = 9
        q(ie)   = 7.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 5.0_rp
        mass(ie)= 186.207_rp
      case('os')
        name(ie) = 'osmium'
        number(ie) = 76
        no(ie) = 9
        q(ie)   = 8.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 6.0_rp
        mass(ie)= 190.23_rp
      case('ir')
        name(ie) = 'iridium'
        number(ie) = 77
        no(ie) = 9
        q(ie)   = 9.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 7.0_rp
        mass(ie)= 192.217_rp
      case('pt')
        name(ie) = 'platinum'
        number(ie) = 78
        no(ie) = 9
        q(ie)   = 10.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 9.0_rp
        mass(ie)= 195.078_rp
      case('au')
        name(ie) = 'gold'
        number(ie) = 79
        no(ie) = 9
        q(ie)   = 11.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 196.9665_rp
      case('hg')
        name(ie) = 'mercury'
        number(ie) = 80
        no(ie) = 9
        q(ie)   = 12.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 200.59_rp
      case('tl')
        name(ie) = 'thallium'
        number(ie) = 81
        no(ie) = 9
        q(ie)   = 13.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 1.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 204.3833_rp
      case('pb')
        name(ie) = 'lead'
        number(ie) = 82
        no(ie) = 9
        q(ie)   = 14.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 2.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 207.2_rp
      case('bi')
        name(ie) = 'bismuth'
        number(ie) = 83
        no(ie) = 9
        q(ie)   = 15.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 3.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 208.9804_rp
      case('po')
        name(ie) = 'polonium'
        number(ie) = 84
        no(ie) = 9
        q(ie)   = 16.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 4.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 209_rp
      case('at')
        name(ie) = 'astatine'
        number(ie) = 85
        no(ie) = 9
        q(ie)   = 17.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 5.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 210_rp
      case('rn')
        name(ie) = 'radon'
        number(ie) = 86
        no(ie) = 9
        q(ie)   = 18.0_rp
        q_s(ie) = 2.0_rp
        q_p(ie) = 6.0_rp
        q_d(ie) = 10.0_rp
        mass(ie)= 222_rp
      case default
        first_c=lower(trim(symbol(ie)))
        if(first_c=='j') then
     ! TB model
        name(ie) = lower(trim(symbol(ie)))
        number(ie) = 0
        no(ie) = 9
        q(ie)   = 1.0_rp
        q_s(ie) = 1.0_rp
        q_p(ie) = 0.0_rp
        q_d(ie) = 0.0_rp
        mass(ie)= 0.0_rp
      else
        write(error_unit,*) 'element%parse_symbol(): element symbol ', &
         symbol(ie), 'not known'
         error stop
      endif
      end select
    end do
    no_max = maxval(no)
  end subroutine initialize_q

  ! Initialize spin-related variables
  subroutine initialize_so(ne, xi_so_p, xi_so_d)
    integer,intent(in) :: ne
    real(rp),dimension(:),allocatable,intent(out) :: xi_so_p, xi_so_d

    allocate(xi_so_p(ne))
    allocate(xi_so_d(ne))

    xi_so_p    = 0.0_rp
    xi_so_d    = 0.0_rp
  end subroutine initialize_so

  ! Parse orbital into spd string format
  function parse_orbital_to_string(o) result(os)
    integer,dimension(:) :: o
    character(len=3) :: os
    integer :: io

    os = '---'
    do io=1,size(o)
      if(o(io)==1) then
        os(1:1) = 's'
      elseif(o(io)>=2 .and. o(io)<=4) then
        os(2:2) = 'p'
      elseif(o(io)>=5 .and. o(io)<=9) then
        os(3:3) = 'd'
      end if
    end do
  end function parse_orbital_to_string

  !> Read object in text format from file (default: 'in_element.txt')
  subroutine read_txt(obj,file)
    class(element),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    integer :: ne
    character(len=2),dimension(:),allocatable :: symbol
    character(len=sl),dimension(:),allocatable :: name
    integer,dimension(:),allocatable :: number
    integer,dimension(:),allocatable :: no
    integer :: no_max
    integer,dimension(:,:),allocatable :: o
    character(len=3),dimension(:),allocatable :: os
    real(rp),dimension(:),allocatable :: q, q_s, q_p, q_d
    real(rp),dimension(:),allocatable :: u_lcn, u_lcn_d
    real(rp),dimension(:),allocatable :: i_stoner_d, b, j_dd, u_dd
    real(rp),dimension(:),allocatable :: xi_so_p, xi_so_d
    real(rp),dimension(:),allocatable :: mass
    ! Namelist
    namelist /element/ ne, symbol, no, o, q, q_s, q_p, q_d, u_lcn, u_lcn_d, &
     i_stoner_d, b, j_dd, u_dd, xi_so_p, xi_so_d, mass
    ! Local variables
    integer :: ie

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_element.txt'
    end if

    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'element%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if

    if(iostatus /= 0) then
      write(error_unit,*) 'element%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    allocate(symbol(0))
    read(10,nml=element,iostat=iostatus)
    deallocate(symbol)
    allocate(symbol(ne))
    allocate(no(0))
    rewind(10)
    read(10,nml=element,iostat=iostatus)
    deallocate(no)
    call initialize_q(ne,symbol,name,number,no,no_max,q,q_s,q_p,q_d,mass)
    call initialize_o(ne,symbol,no_max,o,os)
    call initialize_lcn(ne,u_lcn,u_lcn_d)
    call initialize_eei(ne,i_stoner_d,b,j_dd,u_dd)
    call initialize_so(ne,xi_so_p,xi_so_d)
    rewind(10)
    read(10,nml=element)
    no_max = maxval(no)

    ! Charge initialization recipe
    do ie=1,ne
      os(ie) = parse_orbital_to_string(o(ie,:))
    end do
    do ie=1,ne
      select case(os(ie))
      case('spd')
        q_d(ie) = q_d(ie)
        q_s(ie) = (q(ie)-q_d(ie))/2.0_rp
        q_p(ie) = (q(ie)-q_d(ie))/2.0_rp
      case('--d')
        q_d(ie) = q(ie)
        q_s(ie) = 0.0_rp
        q_p(ie) = 0.0_rp
      case('s--')
        q_d(ie) = 0.0_rp
        q_s(ie) = q(ie)
        q_p(ie) = 0.0_rp
      case('-p-')
        q_d(ie) = 0.0_rp
        q_s(ie) = 0.0_rp
        q_p(ie) = q(ie)
      case('sp-')
        q_d(ie) = 0.0_rp
        q_s(ie) = q(ie)/2.0_rp
        q_p(ie) = q(ie)/2.0_rp
      case('s-d')
        q_d(ie) = q_d(ie)
        q_s(ie) = q(ie)-q_d(ie)
        q_p(ie) = 0.0_rp
      case('-pd')
        q_d(ie) = q_d(ie)
        q_s(ie) = 0.0_rp
        q_p(ie) = q(ie)-q_d(ie)
      end select
    end do

    obj%ne = ne
    call move_alloc(symbol,obj%symbol)
    call move_alloc(name,obj%name)
    call move_alloc(number,obj%number)
    call move_alloc(mass,obj%mass)
    call move_alloc(no,obj%no)
    obj%no_max = no_max
    call move_alloc(o,obj%o)
    call move_alloc(os,obj%os)
    call move_alloc(q,obj%q)
    call move_alloc(q_s,obj%q_s)
    call move_alloc(q_p,obj%q_p)
    call move_alloc(q_d,obj%q_d)
    u_lcn = u_lcn * obj%u%convert_energy('to','hau')
    call move_alloc(u_lcn,obj%u_lcn)
    u_lcn_d = u_lcn_d * obj%u%convert_energy('to','hau')
    call move_alloc(u_lcn_d,obj%u_lcn_d)
    i_stoner_d = i_stoner_d * obj%u%convert_energy('to','hau')
    call move_alloc(i_stoner_d,obj%i_stoner_d)
    b = b * obj%u%convert_energy('to','hau')
    call move_alloc(b,obj%b)
    j_dd = j_dd * obj%u%convert_energy('to','hau')
    call move_alloc(j_dd,obj%j_dd)
    u_dd = u_dd * obj%u%convert_energy('to','hau')
    call move_alloc(u_dd,obj%u_dd)
    xi_so_p = xi_so_p * obj%u%convert_energy('to','hau')
    call move_alloc(xi_so_p,obj%xi_so_p)
    xi_so_d = xi_so_d * obj%u%convert_energy('to','hau')
    call move_alloc(xi_so_d,obj%xi_so_d)

    call obj%calculate_v_osc()
    call obj%calculate_l_dot_s()

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Symmetrize the on-site Coulomb potential
  subroutine symmetrize_v_osc(i,j,k,l,val,v_osc)
    integer,intent(in)  :: i,j,k,l
    real(rp),intent(in) :: val
    real(rp),dimension(5,5,5,5),intent(out) :: v_osc

    v_osc(i,j,k,l) = val
    v_osc(k,j,i,l) = val
    v_osc(i,l,k,j) = val
    v_osc(k,l,i,j) = val
    v_osc(j,i,l,k) = val
    v_osc(l,i,j,k) = val
    v_osc(j,k,l,i) = val
    v_osc(l,k,j,i) = val
  end subroutine symmetrize_v_osc

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_element.txt')
  subroutine write_txt(obj,file,unit)
    class(element),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    integer :: ne
    character(len=2),dimension(obj%ne) :: symbol
    integer,dimension(obj%ne) :: no
    integer,dimension(obj%ne,obj%no_max) :: o
    real(rp),dimension(obj%ne) :: q, q_s, q_p, q_d
    real(rp),dimension(obj%ne) :: u_lcn, u_lcn_d
    real(rp),dimension(obj%ne) :: i_stoner_d, b, j_dd, u_dd
    real(rp),dimension(obj%ne) :: xi_so_p, xi_so_d
    real(rp),dimension(obj%ne) :: mass
    ! Namelist
    namelist /element/ ne, symbol, no, o, q, q_s, q_p, q_d, u_lcn, u_lcn_d, &
     i_stoner_d, b, j_dd, u_dd, xi_so_p, xi_so_d, mass

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_element.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    ne = obj%ne
    symbol = obj%symbol
    no = obj%no
    o = obj%o
    q   = obj%q
    q_s = obj%q_s
    q_p = obj%q_p
    q_d = obj%q_d
    u_lcn   = obj%u_lcn * obj%u%convert_energy('from','hau')
    u_lcn_d = obj%u_lcn_d * obj%u%convert_energy('from','hau')
    i_stoner_d = obj%i_stoner_d * obj%u%convert_energy('from','hau')
    b    = obj%b * obj%u%convert_energy('from','hau')
    j_dd = obj%j_dd * obj%u%convert_energy('from','hau')
    u_dd = obj%u_dd * obj%u%convert_energy('from','hau')
    xi_so_p = obj%xi_so_p * obj%u%convert_energy('from','hau')
    xi_so_d = obj%xi_so_d * obj%u%convert_energy('from','hau')
    mass = obj%mass

    write(unit_rt,nml=element)

    if(.not. present(unit)) then
      call dynamol_flush(unit_rt)
      close(unit_rt)
    else
      call dynamol_flush(unit)
      close(unit)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_element.txt'), if tag (default: .true.) the namelist opening and
  !> closing tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(element),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    real(rp),dimension(obj%ne) :: u_lcn, u_lcn_d
    real(rp),dimension(obj%ne) :: i_stoner_d, b, j_dd, u_dd
    real(rp),dimension(obj%ne) :: xi_so_p, xi_so_d
    real(rp),dimension(obj%ne) :: mass
    ! Local variables
    integer :: ie, io, ip

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_element.txt'
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
      write(unit_rt,'(a)') '&element'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('ne')
        write(unit_rt,'(a)') ' ne = ' // int2str(obj%ne)
      case('symbol')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' symbol(' // int2str(ie) // ') = ''' &
           // trim(obj%symbol(ie)) // ''''
        end do
      case('name')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' name(' // int2str(ie) // ') = ''' &
           // trim(obj%name(ie)) // ''''
        end do
      case('number')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' number(' // int2str(ie) // ') = ' &
           // int2str(obj%number(ie))
        end do
      case('no')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' no(' // int2str(ie) // ') = ' &
           // int2str(obj%no(ie))
        end do
      case('o')
        do ie=1,obj%ne
          write(unit_rt,'(' // int2str(obj%no(ie)) // 'a)') ' o(' // int2str(ie) &
           // ',1:' // int2str(obj%no(ie)) // ') = ' // int2str(obj%o(ie,1)), &
           (', ' // int2str(obj%o(ie,io)), io=2, obj%no(ie))
        end do
      case('q')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' q(' // int2str(ie) // ') = ' &
           // real2str(obj%q(ie))
        end do
      case('q_s')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' q_s(' // int2str(ie) // ') = ' &
           // real2str(obj%q_s(ie))
        end do
      case('q_p')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' q_p(' // int2str(ie) // ') = ' &
           // real2str(obj%q_p(ie))
        end do
      case('q_d')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' q_d(' // int2str(ie) // ') = ' &
           // real2str(obj%q_d(ie))
        end do
      case('u_lcn')
        u_lcn = obj%u_lcn * obj%u%convert_energy('from','hau')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' u_lcn(' // int2str(ie) // ') = ' &
           // real2str(u_lcn(ie))
        end do
      case('u_lcn_d')
        u_lcn_d = obj%u_lcn_d * obj%u%convert_energy('from','hau')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' u_lcn_d(' // int2str(ie) // ') = ' &
           // real2str(u_lcn_d(ie))
        end do
      case('i_stoner_d')
        i_stoner_d = obj%i_stoner_d * obj%u%convert_energy('from','hau')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' i_stoner_d(' // int2str(ie) // ') = ' &
           // real2str(i_stoner_d(ie))
        end do
      case('b')
        b = obj%b * obj%u%convert_energy('from','hau')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' b(' // int2str(ie) // ') = ' &
           // real2str(b(ie))
        end do
      case('j_dd')
        j_dd = obj%j_dd * obj%u%convert_energy('from','hau')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' j_dd(' // int2str(ie) // ') = ' &
           // real2str(j_dd(ie))
        end do
      case('u_dd')
        u_dd = obj%u_dd * obj%u%convert_energy('from','hau')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' u_dd(' // int2str(ie) // ') = ' &
           // real2str(u_dd(ie))
        end do
      case('xi_so_p')
        xi_so_p = obj%xi_so_p * obj%u%convert_energy('from','hau')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' xi_so_p(' // int2str(ie) // ') = ' &
           // real2str(xi_so_p(ie))
        end do
      case('xi_so_d')
        xi_so_d = obj%xi_so_d * obj%u%convert_energy('from','hau')
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' xi_so_d(' // int2str(ie) // ') = ' &
           // real2str(xi_so_d(ie))
        end do
      case('mass')
        mass = obj%mass
        do ie=1,obj%ne
          write(unit_rt,'(a)') ' mass(' // int2str(ie) // ') = ' &
           // real2str(mass(ie))
        end do
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
    !deallocate(file_rt,property_rt)
  end subroutine write_txt_formatted
end module element_mod
