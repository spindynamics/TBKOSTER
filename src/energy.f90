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
!  energy.f90
!  TBKOSTER
module energy_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_mod
  use charge_mod
  use element_mod
  use hamiltonian_tb_mod
  use math_mod
  use mesh_mod
  use precision_mod, only: rp
  use string_mod, only: TBKOSTER_flush, int2str, log2str, lower, real2str, sl
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'smearing', &
   'degauss', &
   'fixed_fermi_level', &
   'en_f_ffl', &
   'fixed_spin_moment', &
   'm_fsm', &
   'en_min', &
   'en_max' &
   ]

  type,public :: energy
    ! Units
    class(units),pointer :: u
    ! Elements
    class(element),pointer :: e
    ! Atoms
    class(atom),pointer :: a
    ! Charges
    class(charge),pointer :: q
    ! Reciprocal space mesh
    class(mesh),pointer :: k
    ! Hamiltonian
    class(hamiltonian_tb),pointer :: h
    !> @defgroup Smearing Smearing related variables
    !> @{

    !> Smearing type ; options:
    !>  'fd': Fermi-Dirac
    !>  'g' : Gaussian
    !>  'mp': Methfessel-Paxton
    !>  'mv': Marzari-Vanderbilt
    character(len=2) :: smearing
    !> Gaussian spreading
    real(rp) :: degauss
    !> @}

    !> @defgroup Fixed_Fermi_level Fixed Fermi level related variables
    !> @{

    ! Fixed Fermi level flag
    logical :: fixed_fermi_level
    ! Fixed Fermi level
    real(rp) :: en_f_ffl
    !> @}

    !> @defgroup Fixed-spin-moment Fixed-spin-moment related variables
    !> @{

    !> Fixed-spin-moment.\n
    !> See \cite Schwarz1984 page 2-3.\n
    !> In systems such as transition metals where the orbital contributions to the
    !> magnetisation can be neglected, the latter is mostly determined by the spin
    !> moment \f$ M \f$.\n
    !> The number of valence electrons \f$ N \f$ and the spin
    !> moment per unit cell \f$ M \f$ are given by:
    ! \f$ N = N^\uparrow + N^\downarrow = \int_{-\infty}^{E_{\mathrm{F}}^\uparrow}
    ! n^\uparrow(E) dE + \int_{-\infty}^{E_{\mathrm{F}}^\downarrow}
    ! n^\downarrow(E) dE \f$\n
    ! \f$ M = N^\uparrow - N^\downarrow = \int_{-\infty}^{E_{\mathrm{F}}^\uparrow}
    ! n^\uparrow(E) dE - \int_{-\infty}^{E_{\mathrm{F}}^\downarrow}
    ! n^\downarrow(E) dE \f$\n
    !> \f{eqnarray}{
    !>   N = N^\uparrow + N^\downarrow = \int_{-\infty}^{E_{\mathrm{F}}^\uparrow}
    !>    n^\uparrow(E) dE + \int_{-\infty}^{E_{\mathrm{F}}^\downarrow}
    !>    n^\downarrow(E) dE\\
    !>   M = N^\uparrow - N^\downarrow = \int_{-\infty}^{E_{\mathrm{F}}^\uparrow}
    !>    n^\uparrow(E) dE - \int_{-\infty}^{E_{\mathrm{F}}^\downarrow}
    !>    n^\downarrow(E) dE
    !> \f}
    !> where \f$ N \f$ represents electron numbers, \f$ n \f$ densities of states
    !> (DOS), \f$ E_{\mathrm{F}} \f$ Fermi energies, \f$ \uparrow \f$ stands for
    !> majority spin electrons and \f$ \downarrow \f$ for minority spin electrons.
    !> \n
    !> To obtain the ground state of a system, one traditionally requires \f$
    !> E_{\mathrm{F}}^\uparrow = E_{\mathrm{F}}^\downarrow \f$. \f$ N \f$ alone
    !> defines \f$ E_{\mathrm{F}} \f$ in equation (1) and \f$ M \f$ is thus
    !> determined by equation (2).\n
    !> Alternatively, relaxing the \f$ E_{\mathrm{F}}^\uparrow =
    !> E_{\mathrm{F}}^\downarrow \f$ constraint allows one to manually set
    !> \f$ M \f$: this is the so-called fixed-spin-moment method.\n
    !> This time, \f$ N^\uparrow \f$ and \f$ N^\downarrow \f$ are obtained from
    !> the first parts of equations (1) and (2) as functions of \f$ N \f$ and
    !> \f$ M \f$. The integrals in the second part of equation (1) then determine
    !> \f$ E_{\mathrm{F}}^\uparrow \f$ and \f$ E_{\mathrm{F}}^\downarrow \f$.\n
    !> Such a method allows one to perform self-consistent spin-polarized
    !> calculations for a given \f$ M \f$ without being restricted to the ground
    !> state. The obtained lowest total energy state will generally have \f$
    !> E_{\mathrm{F}}^\uparrow \neq E_{\mathrm{F}}^\downarrow \f$.
    ! Fixed-spin-moment flag
    logical :: fixed_spin_moment
    !> Fixed-spin-moment magnetization \f$ M^{\mathrm{FSM}} \f$
    real(rp) :: m_fsm
    !> Fixed-spin-moment Fermi level \f$ E_{\mathrm{F}}^{\mathrm{FSM}} \f$
    real(rp),dimension(2) :: en_f_fsm
    !> @}

    ! Energy lower bound (for band structure visualisation)
    real(rp) :: en_min
    ! Energy upper bound (for band structure visualisation)
    real(rp) :: en_max

    !> @defgroup Eigenvalues Eigenvalues related variables
    !> @{

    !> Number of eigenvalues
    integer :: nen_k
    !> Eigenvalues
    real(rp),dimension(:)  ,allocatable :: en_k
    real(rp),dimension(:,:,:),allocatable :: en_k_2
    !> Eigenvalues per spin
    real(rp),dimension(:,:),allocatable :: en_k_fsm
    !> Sorted eigenvalue index
    integer, dimension(:)  ,allocatable :: indx
    !> Sorted eigenvalue index per spin
    integer, dimension(:,:),allocatable :: indx_fsm
    !> Energy filling f(E_n-E_f)
    real(rp),dimension(:)  ,allocatable :: f_k
    !> Energy filling per spin f(E_n-E_f)
    real(rp),dimension(:,:),allocatable :: f_k_fsm
    !> @}

    !> Fermi level
    real(rp) :: en_f

    !> @defgroup Energy Energy related variables
    !> @{

    !> Band energy
    real(rp) :: en_band
    real(rp) :: en_band_f
    !> Local band energy
    real(rp),dimension(:,:,:),allocatable :: en_band_local
    real(rp),dimension(:,:,:),allocatable :: en_band_local_f
    !> Entropy times temperature (for FD smearing)
    real(rp) :: s_t
    !> Double counting energies: EEI, LCN, magnetic penalization
    real(rp) :: en_dc_eei
    real(rp) :: en_dc_lcn
    real(rp) :: en_dc_pen
    !> Total energy
    real(rp) :: en_in, en_out
    !> Total energy difference
    real(rp) :: delta_en
    !> @}

    !> @defgroup Local_band_energy_weight Local band energy weight-related
    !> variables
    !> @{
  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: calculate_delta_en
    procedure :: calculate_en
    procedure :: calculate_en_band
    procedure :: calculate_en_band_local_k
    procedure :: calculate_en_dc_eei
    procedure :: calculate_en_dc_lcn
    procedure :: calculate_en_dc_pen
    procedure :: calculate_en_f
    procedure :: calculate_en_f_std
    procedure :: calculate_en_f_fsm
    procedure :: initialize
    procedure :: is_converged
    procedure :: print_indx
    procedure :: read_txt
    procedure :: save_en_k
    procedure :: sort_en_k
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type energy

  ! Constructor
  interface energy
    procedure :: constructor
  end interface energy

contains
  function constructor(h) result(obj)
    class(hamiltonian_tb),target,intent(in) :: h
    type(energy) :: obj

    obj%u => h%u
    obj%e => h%e_tb
    obj%a => h%a_tb
    obj%q => h%q
    obj%k => h%k
    obj%h => h
  end function constructor

  subroutine destructor(obj)
    type(energy) :: obj

    if(allocated(obj%indx))            deallocate(obj%indx)
    if(allocated(obj%indx_fsm))        deallocate(obj%indx_fsm)
    if(allocated(obj%en_k))            deallocate(obj%en_k)
    if(allocated(obj%en_k_2))          deallocate(obj%en_k_2)
    if(allocated(obj%en_k_fsm))        deallocate(obj%en_k_fsm)
    if(allocated(obj%f_k))             deallocate(obj%f_k)
    if(allocated(obj%f_k_fsm))         deallocate(obj%f_k_fsm)
    if(allocated(obj%en_band_local))   deallocate(obj%en_band_local)
    if(allocated(obj%en_band_local_f)) deallocate(obj%en_band_local_f)
  end subroutine destructor

  subroutine calculate_delta_en(obj)
    class(energy),intent(inout) :: obj

    obj%delta_en = abs(obj%en_out-obj%en_in)
  end subroutine calculate_delta_en

  subroutine calculate_en(obj)
    class(energy),intent(inout) :: obj

    obj%en_in  = obj%en_out
    obj%en_out = obj%en_band - obj%en_dc_eei - obj%en_dc_lcn - obj%en_dc_pen
  end subroutine calculate_en

  subroutine calculate_en_band(obj)
    ! INPUT
    class(energy),intent(inout) :: obj
    !	LOCAL
    real(rp) :: ff,ff_i
    integer :: nn,ik,kspin,kmat,jj

    obj%en_band   = 0.0_rp
    obj%en_band_f = 0.0_rp

    if(obj%fixed_spin_moment) then
      do kspin=1,2
        do ik=1,obj%k%nx
          do kmat=1,obj%h%nh

            jj = ik+(kmat-1)*obj%k%nx
            ff = obj%f_k_fsm(jj,kspin)
            obj%en_band = obj%en_band + ff*obj%k%w(ik)*obj%en_k_fsm(jj,kspin)
            ff_i = integrated_delta_function(-(obj%en_k_fsm(jj,kspin) &
             - obj%en_f_fsm(kspin))/obj%degauss, obj%smearing)
            obj%s_t = obj%s_t-obj%k%w(ik)*ff_i

          end do
        end do
      end do
    else
      do nn=1,obj%nen_k
        ik=mod(nn,obj%k%nx)
        if(ik==0) ik=obj%k%nx

        ff = obj%f_k(nn)
        obj%en_band   = obj%en_band + ff*obj%en_k(nn)*obj%a%g_s*obj%k%w(ik)
        obj%en_band_f = obj%en_band_f + ff*(obj%en_k(nn)-obj%en_f) &
         *obj%a%g_s*obj%k%w(ik)
        ff_i = integrated_delta_function(-(obj%en_k(nn)-obj%en_f) &
         /obj%degauss, obj%smearing)
        obj%s_t = obj%s_t-obj%k%w(ik)*ff_i*obj%a%g_s
      end do
    end if
  end subroutine calculate_en_band

  subroutine calculate_en_band_local_k(obj,ik,isl,v_k)
    ! INPUT
    class(energy),intent(inout) :: obj
    integer :: ik,isl
    complex(rp),intent(in),dimension(2,obj%h%nh,obj%h%nh) :: v_k
    ! LOCAL
    integer :: ia,ie,io,l,ispin,jspin,imat_ispin,imat_jspin,imat,jmat,jmat2,nn,q
    real(rp) :: ffe


    select case(obj%a%ns)
    case(1,2)
      do ia=1,obj%a%na
        ie = obj%a%ia2ie(ia)
        do io=1,obj%e%no(ie)
          imat = obj%h%iaos2ih(ia,io,1)
          do jmat=1,obj%h%nh
            jmat2 = isl+(jmat-1)*obj%a%ns
            nn = ik+(jmat2-1)*obj%k%nx
            ffe = obj%f_k(nn)*obj%en_k(nn)
            obj%en_band_local(ia,io,isl) = obj%en_band_local(ia,io,isl) &
             + ffe*obj%k%w(ik)*obj%a%g_s*real(v_k(1,imat,jmat) &
             *conjg(v_k(2,imat,jmat)))
            ffe = ffe - obj%en_f*obj%f_k(nn)
            obj%en_band_local_f(ia,io,isl) = obj%en_band_local_f(ia,io,isl) &
             + ffe*obj%k%w(ik)*obj%a%g_s*real(v_k(1,imat,jmat) &
             *conjg(v_k(2,imat,jmat)))
          end do
        end do
      end do

    case(4)
      do ia=1,obj%a%na
        ie = obj%a%ia2ie(ia)
        do io=1,obj%e%no(ie)
          do ispin=1,2

            imat_ispin = obj%h%iaos2ih(ia,io,ispin)
            
            do jspin=1,2

              imat_jspin = obj%h%iaos2ih(ia,io,jspin)

              q = obj%a%iss2is(ispin,jspin)

              do jmat=1,obj%h%nh
                nn = ik+(jmat-1)*obj%k%nx

                ffe = obj%f_k(nn)*obj%en_k(nn)

                obj%en_band_local(ia,io,q) = obj%en_band_local(ia,io,q) &
                 + ffe*obj%k%w(ik) &
                 *real(conjg(v_k(1,imat_jspin,jmat))*v_k(2,imat_ispin,jmat) &
                 + v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)))/2

                ffe = ffe-obj%f_k(nn)*obj%en_f

                obj%en_band_local_f(ia,io,q) = obj%en_band_local_f(ia,io,q) &
                 + ffe*obj%k%w(ik) &
                 *real(conjg(v_k(1,imat_jspin,jmat))*v_k(2,imat_ispin,jmat)&
                 + v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)))/2
              end do

            end do
          end do
        end do
      end do
    end select
  end subroutine calculate_en_band_local_k

  subroutine calculate_en_dc_eei(obj)
    ! INPUT
    class(energy),intent(inout) :: obj
    ! LOCAL
    integer :: ia,ie,io1,io2,ispin,jspin

    obj%en_dc_eei = 0.0_rp

    select case(obj%a%ns)
    case(1,2)
      do ia=1,obj%a%na
        ie = obj%a%ia2ie(ia)
        do io1=1,obj%e%no(ie)
          do io2=1,obj%e%no(ie)
            do ispin=1,obj%a%ns
              obj%en_dc_eei = obj%en_dc_eei + 0.5_rp &
               *real(obj%h%delta_h_eei(ia,io1,io2,ispin) &
               *obj%q%rho_net_out(ia,io2,io1,ispin))
            end do
          end do
        end do
      end do
    case(4)
      do ia=1,obj%a%na
        ie = obj%a%ia2ie(ia)
        do io1=1,obj%e%no(ie)
          do io2=1,obj%e%no(ie)
            do ispin=1,2
              do jspin=1,2
                obj%en_dc_eei = obj%en_dc_eei + 0.5_rp &
                 *real(obj%h%delta_h_eei(ia,io1,io2,obj%a%iss2is(ispin,jspin)) &
                 *obj%q%rho_net_out(ia,io2,io1,obj%a%iss2is(jspin,ispin)))
              end do
            end do
          end do
        end do
      end do
    end select
  end subroutine calculate_en_dc_eei

  subroutine calculate_en_dc_lcn(obj)
    ! INPUT
    class(energy),intent(inout) :: obj
    ! LOCAL
    real(rp) :: n_s, n_p, n_d, n_tot, dn
    integer :: ia,ie

    obj%en_dc_lcn = 0.0_rp

    select case(obj%a%ns)
    case(1,4)
      do ia=1,obj%a%na
        ie = obj%a%ia2ie(ia)

        n_s = obj%q%q_mul_out(ia,1,0)
        n_p = obj%q%q_mul_out(ia,2,0)
        n_d = obj%q%q_mul_out(ia,3,0)
        n_tot = n_s + n_p + n_d
        dn = n_tot - obj%e%q(ie)

        obj%en_dc_lcn = obj%en_dc_lcn &
         + obj%e%u_lcn(ie)/2*(n_tot**2-obj%e%q(ie)**2) &
         + obj%e%u_lcn_d(ie)/2*(n_d**2-(obj%e%q_d(ie))**2)
      end do
    case(2)
      do ia=1,obj%a%na
        ie = obj%a%ia2ie(ia)

        n_s = obj%q%q_mul_out(ia,1,0) + obj%q%q_mul_out(ia,1,1)
        n_p = obj%q%q_mul_out(ia,2,0) + obj%q%q_mul_out(ia,2,1)
        n_d = obj%q%q_mul_out(ia,3,0) + obj%q%q_mul_out(ia,3,1)
        n_tot = n_s + n_p + n_d
        dn = n_tot - obj%e%q(ie)

        obj%en_dc_lcn = obj%en_dc_lcn &
         + obj%e%u_lcn(ie)/2*(n_tot**2-obj%e%q(ie)**2) &
         + obj%e%u_lcn_d(ie)/2*(n_d**2-(obj%e%q_d(ie))**2)
      end do
    end select
  end subroutine calculate_en_dc_lcn

  subroutine calculate_en_dc_pen(obj)
    ! INPUT
    class(energy),intent(inout) :: obj
    ! LOCAL
    integer :: ia,l
    real(rp) :: m_z
    real(rp),dimension(3) :: m_cart

    obj%en_dc_pen = 0.0_rp

    select case(obj%a%ns)
    case(2)
      select case(obj%h%m_penalization)
      case('none')
        ! Do nothing
      case('r')
        do ia=1,obj%a%na
          m_z = 0.0_rp
          do l=1,3
            m_z = m_z + obj%q%q_mul_out(ia,l,0) - obj%q%q_mul_out(ia,l,1)
          end do

          obj%en_dc_pen = obj%en_dc_pen + obj%a%lambda_pen(ia) &
           *(m_z**2 - obj%a%m_pen(ia,1)**2)
        end do
      end select
    case(4)
      select case(obj%h%m_penalization)
      case('none')
        ! Do nothing
      case('r,theta,phi')
        do ia=1,obj%a%na
          m_cart = sum(obj%q%q_mul_out(ia,:,1:3),1)

          obj%en_dc_pen = obj%en_dc_pen + obj%a%lambda_pen(ia) &
           *(dot_product(m_cart,m_cart) - obj%a%m_pen(ia,1)**2)
        end do
      end select
    end select
  end subroutine calculate_en_dc_pen

  subroutine calculate_en_f(obj)
    ! INPUT
    class(energy),intent(inout) :: obj

    if(obj%fixed_spin_moment) then
      call obj%calculate_en_f_fsm()
    else
      call obj%calculate_en_f_std()
    end if
  end subroutine calculate_en_f

  ! Calculate the Fermi level with broadening technique
  subroutine calculate_en_f_std(obj)
    ! INPUT
    class(energy),intent(inout) :: obj
    ! LOCAL
    real(rp) :: en_0,en_1,en_mid,en_f_init,ff,ffd,dos_en_0
    real(rp) :: dos_int,dos_int_pred
    integer  :: ii,jj,ik,nn
    integer  :: n_homo,ndeg_homo,itest,n_homo_min,n_homo_max,ne_homo

    if(obj%degauss>0.0_rp) then  ! if beta>0 means fermi broadening... otherwise discrete levels
      dos_int = 0.0_rp      ! integrated number of electrons
      !--------------------------------------------------------
      !         first evaluation of the fermi level
      !         the levels are filled in increasing order
      !         the total number of electrons is equal to nel
      !---------------------------------------------------------
      jj = 1

      do
        ik = mod(obj%indx(jj),obj%k%nx) ! ik is the k point index of the corresponding eigenvalue
        if(ik==0) ik = obj%k%nx
        dos_int_pred = dos_int
        dos_int = dos_int + obj%a%g_s*obj%k%w(ik)
        if(dos_int_pred<obj%a%nel .and. dos_int>=obj%a%nel) exit
        jj = jj+1
      end do

      en_f_init = obj%en_k(obj%indx(jj))


      !------------------------------------------------------
      !        more precise evaluation of the fermi level
      !        with a broadening technique.
      !------------------------------------------------------
      select case(obj%smearing)
      case('fd')

        en_0 = en_f_init
        dos_int = 0.0_rp
        do while(abs(dos_int-obj%a%nel)>epsilon)
          dos_en_0 = 0.0_rp
          dos_int = 0.0_rp
          do jj=1,obj%nen_k
            ik = mod(obj%indx(jj),obj%k%nx)
            if(ik==0) ik=obj%k%nx

            !****|******************************************************************|
            ff  = theta_function(-(obj%en_k(obj%indx(jj))-en_0)/obj%degauss,obj%smearing)
            ffd = delta_function(-(obj%en_k(obj%indx(jj))-en_0)/obj%degauss,obj%smearing)/obj%degauss
            !****|******************************************************************|
            dos_int = dos_int + obj%k%w(ik)*ff*obj%a%g_s
            dos_en_0 = dos_en_0 + obj%k%w(ik)*ffd*obj%a%g_s
          end do
          en_0 = en_0 + (obj%a%nel-dos_int)/dos_en_0  ! we use a Newton-like algorithm to determine the Fermi Level: N(en_f)-nel=0
        end do

        !
        ! en_f: precise fermi level
        !
        obj%en_f = en_0

      case('g','mp','mv')

        en_0 = en_f_init - 5*obj%degauss
        en_1 = en_f_init + 5*obj%degauss
        en_mid = (en_0+en_1)/2.0_rp
        dos_int = 0.0_rp
        do while(abs(dos_int-obj%a%nel)>epsilon)
          dos_en_0 = 0.0_rp
          dos_int = 0.0_rp

          do jj=1,obj%nen_k
            ik = mod(obj%indx(jj),obj%k%nx)
            if(ik==0) ik = obj%k%nx

            !****|******************************************************************|
            ff  = theta_function(-(obj%en_k(obj%indx(jj))-en_mid)/obj%degauss,obj%smearing)
            ffd = delta_function(-(obj%en_k(obj%indx(jj))-en_mid)/obj%degauss,obj%smearing)/obj%degauss
            !****|******************************************************************|
            dos_int = dos_int + obj%k%w(ik)*ff*obj%a%g_s
            dos_en_0 = dos_en_0 + obj%k%w(ik)*ffd*obj%a%g_s
          end do

          if((dos_int-obj%a%nel)>0) then
            en_1 = en_mid
            en_mid = (en_0+en_1)/2.0_rp
          else
            en_0 = en_mid
            en_mid = (en_0+en_1)/2.0_rp
          end if
        end do

        !
        ! en_f: precise fermi level
        !
        obj%en_f = en_mid
      end select
      !
      if(obj%fixed_fermi_level) then
        obj%en_f = obj%en_f_ffl
      end if

      do nn=1,obj%nen_k
        ! calculate the band filling with fermi broadening
        obj%f_k(nn) = theta_function(-(obj%en_k(nn)-obj%en_f) &
         /obj%degauss,obj%smearing)
      end do

    elseif(obj%degauss==0.0_rp .and. obj%k%nx==1) then ! case without fermi broadening -> (homo/lumo) case

      n_homo = int(obj%a%nel/obj%a%g_s)
      obj%en_f = obj%en_k(obj%indx(n_homo))
      ndeg_homo = 1
      ii = 0
      itest = 0
      do while(itest==0)
        ii=ii+1
        if(abs(obj%en_k(obj%indx(n_homo)) - obj%en_k(obj%indx(n_homo-ii))) &
         < epsilon) then
          ndeg_homo = ndeg_homo+1
          itest = 0
        else
          itest = 1
        end if
      end do
      n_homo_min = n_homo - ii + 1
      ii = 0
      itest = 0
      do while(itest==0)
        ii = ii+1
        if(abs(obj%en_k(obj%indx(n_homo)) - obj%en_k(obj%indx(n_homo+ii))) &
         < epsilon) then
          ndeg_homo = ndeg_homo+1
          itest = 0
        else
          itest = 1
        end if
      end do
      n_homo_max = n_homo + ii - 1

      do nn=1,obj%nen_k
        if(nn < n_homo_min) then
          obj%f_k(obj%indx(nn)) = 1.0_rp
        elseif(nn > n_homo_max) then
          obj%f_k(obj%indx(nn)) = 0.0_rp
        elseif(nn>=n_homo_min .and. nn<=n_homo_max) then
          ne_homo = nint(obj%a%nel/obj%a%g_s) - n_homo_min + 1
          obj%f_k(obj%indx(nn)) = real(ne_homo)/ndeg_homo
        end if
      end do
    end if
  end subroutine calculate_en_f_std

  ! Calculate the Fermi level with broadening technique for FSM method
  subroutine calculate_en_f_fsm(obj)
    ! INPUT
    class(energy),intent(inout) :: obj
    ! LOCAL
    real(rp),dimension(2) :: nel_fsm
    real(rp) :: en_0,en_1,en_mid,en_f_init(2),ff,ffd,dos_en_0,dos_en_f_tot
    real(rp) :: dos_int,dos_int_pred
    integer  :: ii,ispin,jmat,jj,ik,nn,jmat2
    integer  :: n_homo,ndeg_homo,itest,n_homo_min,n_homo_max,ne_homo

    nel_fsm(1) = (obj%a%nel+obj%m_fsm)/2.0_rp
    nel_fsm(2) = (obj%a%nel-obj%m_fsm)/2.0_rp
    !--------------------------------------------------------
    !         first evaluation of the fermi level
    !         the levels are filled in increasing order
    !         up to the nel level (total number of electrons)
    !---------------------------------------------------------
    if(obj%degauss>0.0_rp) then

      do ispin=1,2

        dos_int=0.0_rp

        jj=1
        do
          jmat = obj%indx_fsm(jj,ispin)/obj%k%nx+1
          ik = mod(obj%indx_fsm(jj,ispin),obj%k%nx)
          if(ik==0) then
            ik = obj%k%nx
            jmat = jmat-1
          end if
          dos_int_pred = dos_int
          dos_int = dos_int + obj%k%w(ik)
          if(dos_int_pred<nel_fsm(ispin) .and. dos_int>=nel_fsm(ispin)) exit
          jj = jj+1
        end do

        en_f_init(ispin) = obj%en_k_fsm(obj%indx_fsm(jj,ispin),ispin)
      end do

      !------------------------------------------------------
      !        more precise evaluation of the fermi level
      !        with a broadening technique.
      !------------------------------------------------------

      do ispin=1,2

        select case(obj%smearing)
        case('fd')
          en_0 = en_f_init(ispin)
          dos_int = 0.0_rp

          do while(abs(dos_int-nel_fsm(ispin))>epsilon)
            dos_int = 0.0_rp
            dos_en_0 = 0.0_rp
            do jj=1,obj%nen_k
              jmat = obj%indx_fsm(jj,ispin)/obj%k%nx+1
              ik = mod(obj%indx_fsm(jj,ispin),obj%k%nx)
              if(ik==0) then
                ik = obj%k%nx
                jmat = jmat-1
              end if

              !****|******************************************************************|
              ff  = theta_function(-(obj%en_k_fsm(obj%indx_fsm(jj,ispin),ispin)-en_0)/obj%degauss,obj%smearing)
              ffd = delta_function(-(obj%en_k_fsm(obj%indx_fsm(jj,ispin),ispin)-en_0)/obj%degauss,obj%smearing)/obj%degauss
              !****|******************************************************************|
              dos_int = dos_int + obj%k%w(ik)*ff
              dos_en_0 = dos_en_0 + obj%k%w(ik)*ffd

            end do
            en_0 = en_0 + (nel_fsm(ispin)-dos_int)/dos_en_0

          end do  !end do while

          dos_en_f_tot = dos_en_0
          obj%en_f_fsm(ispin) = en_0    !en_f: precise fermi level

        case('g','mp','mv')
          en_0 = en_f_init(ispin)-5*obj%degauss
          en_1 = en_f_init(ispin)+5*obj%degauss
          en_mid = (en_0+en_1)/2._rp
          dos_int = 0.0_rp
          do while(abs(dos_int-nel_fsm(ispin))>epsilon)
            dos_int = 0.0_rp
            dos_en_0 = 0.0_rp
            do jj=1,obj%nen_k
              jmat = obj%indx_fsm(jj,ispin)/obj%k%nx+1
              ik = mod(obj%indx_fsm(jj,ispin),obj%k%nx)
              if(ik==0) then
                ik = obj%k%nx
                jmat = jmat-1
              end if

              !****|******************************************************************|
              ff  = theta_function(-(obj%en_k_fsm(obj%indx_fsm(jj,ispin),ispin)-en_mid)/obj%degauss,obj%smearing)
              ffd = delta_function(-(obj%en_k_fsm(obj%indx_fsm(jj,ispin),ispin)-en_mid)/obj%degauss,obj%smearing)/obj%degauss
              !****|******************************************************************|
              dos_int = dos_int + obj%k%w(ik)*ff
              dos_en_0 = dos_en_0 + obj%k%w(ik)*ffd

            end do
            if((dos_int-nel_fsm(ispin))>0) then
              en_1 = en_mid
              en_mid = (en_0+en_1)/2._rp
            else
              en_0 = en_mid
              en_mid = (en_0+en_1)/2._rp
            end if
          end do
          dos_en_f_tot = dos_en_0
          obj%en_f_fsm(ispin) = en_mid

        end select
          !        write(*,*) 'en_f     =', obj%en_f_fsm(ispin)
          !        write(*,*) 'dos_en_f =', dos_en_f_tot
          !        write(*,*) 'dos_int=', dos_int
        do ik=1,obj%k%nx
          do jmat=1,obj%nen_k/obj%k%nx
            jmat2 = ispin + (jmat-1)*2
            jj = ik + (jmat-1)*obj%k%nx
            nn = ik + (jmat2-1)*obj%k%nx
            obj%f_k(nn) = theta_function(-(obj%en_k_fsm(jj,ispin) &
             - obj%en_f_fsm(ispin))/obj%degauss,obj%smearing)
            obj%f_k_fsm(jj,ispin) = theta_function(-(obj%en_k_fsm(jj,ispin) &
             - obj%en_f_fsm(ispin))/obj%degauss,obj%smearing)
          end do
        end do

      end do  ! end do ispin

    elseif(obj%degauss==0.0_rp .and. obj%k%nx==1) then

      do ispin=1,2

        if(ispin==1) write(*,*) 'spin up'
        if(ispin==2) write(*,*) 'spin down'

        n_homo = int(nel_fsm(ispin))
        obj%en_f_fsm(ispin) = obj%en_k_fsm(obj%indx_fsm(n_homo,ispin),ispin)

        ndeg_homo = 1
        ii = 0
        itest = 0
        do while(itest==0)
          ii=ii+1
          if(abs(obj%en_k_fsm(obj%indx_fsm(n_homo,ispin),ispin) &
           - obj%en_k_fsm(obj%indx_fsm(n_homo-ii,ispin),ispin)) < epsilon) then
            ndeg_homo = ndeg_homo + 1
            itest = 0
          else
            itest = 1
          end if
        end do
        n_homo_min = n_homo - ii + 1
        ii = 0
        itest = 0
        do while(itest==0)
          ii=ii+1
          if(abs(obj%en_k_fsm(obj%indx_fsm(n_homo,ispin),ispin) &
           - obj%en_k_fsm(obj%indx_fsm(n_homo+ii,ispin),ispin)) < epsilon) then
            ndeg_homo = ndeg_homo + 1
            itest = 0
          else
            itest = 1
          end if
        end do
        n_homo_max = n_homo + ii - 1

        do ik=1,obj%k%nx
          do jmat=1,obj%nen_k/obj%k%nx
            jmat2 = ispin + (jmat-1)*2
            jj = ik + (jmat-1)*obj%k%nx
            nn = ik + (jmat2-1)*obj%k%nx

            obj%f_k(nn) = theta_function(-(obj%en_k_fsm(jj,ispin) &
             - obj%en_f_fsm(ispin))/obj%degauss,obj%smearing)
            obj%f_k_fsm(jj,ispin) = theta_function(-(obj%en_k_fsm(jj,ispin) &
             - obj%en_f_fsm(ispin))/obj%degauss,obj%smearing)

            if(jj<n_homo_min) then
              obj%f_k(nn) = 1.0_rp
              obj%f_k_fsm(jj,ispin) = 1.0_rp
            elseif(jj>n_homo_max) then
              obj%f_k(nn) = 0.0_rp
              obj%f_k_fsm(jj,ispin) = 0.0_rp
            elseif(nn>=n_homo_min .and. nn<=n_homo_max) then
              ne_homo = nint(nel_fsm(ispin)) - n_homo_min + 1
              obj%f_k(nn) = real(ne_homo)/ndeg_homo
              obj%f_k_fsm(jj,ispin) = real(ne_homo)/ndeg_homo
            end if
          end do
        end do

        write(*,*) 'E homo                =', obj%en_k_fsm(n_homo,ispin)
        write(*,*) 'E lumo                =', obj%en_k_fsm(n_homo_max+1,ispin)
        write(*,*) 'number of homo states =', ndeg_homo
        write(*,*) 'filling of homo states=', real(ne_homo)/ndeg_homo

      end do
    end if

    obj%en_f = (obj%en_f_fsm(1)+obj%en_f_fsm(2))/2.0_rp
  end subroutine calculate_en_f_fsm

  subroutine check_fixed_spin_moment(ns,fixed_fermi_level,fixed_spin_moment)
    integer,intent(in) :: ns
    logical,intent(in) :: fixed_fermi_level, fixed_spin_moment

    if((ns==1 .or. ns==4) .and. fixed_spin_moment) then
      write(error_unit,*) 'energy%check_fixed_spin_moment(): Cannot do &
       &fixed-spin-moment for spin-unpolarized or spin-polarized noncollinear &
       &calculations, only for spin-polarized collinear calculations'
      error stop
    end if

    if(fixed_fermi_level .and. fixed_spin_moment) then
      write(error_unit,*) 'energy%check_fixed_spin_moment(): Cannot do &
       &fixed-spin-moment for fixed Fermi level calculations'
      error stop
    end if
  end subroutine check_fixed_spin_moment

  subroutine check_smearing(smearing)
    character(len=*),intent(in) :: smearing

    if(smearing /= 'fd' &
     .and. smearing /= 'g' &
     .and. smearing /= 'mp' &
     .and. smearing /= 'mv') then
      write(error_unit,*) 'energy%check_smearing(): energy%smearing must be &
       &one of ''fd'', ''g'', ''mp'', ''mv'''
      error stop
    end if
  end subroutine check_smearing

  subroutine initialize(obj)
    class(energy),intent(inout) :: obj

    if(obj%fixed_spin_moment) then
      obj%nen_k = obj%k%nx*obj%h%nh
    else
      obj%nen_k = obj%k%nx*obj%h%nh*obj%a%nsl
    end if

    ! Eigenvalues and energy fillings
    if(allocated(obj%en_k)) deallocate(obj%en_k)
    allocate(obj%en_k(0:obj%k%nx*obj%h%nh*obj%a%nsl))
    if(allocated(obj%en_k_2)) deallocate(obj%en_k_2)
    allocate(obj%en_k_2(obj%h%nh,obj%k%nx,obj%a%nsl))
    if(allocated(obj%f_k)) deallocate(obj%f_k)
    allocate(obj%f_k(obj%k%nx*obj%h%nh*obj%a%nsl))
    if(.not. obj%fixed_spin_moment) then
      if(allocated(obj%indx)) deallocate(obj%indx)
      allocate(obj%indx(0:obj%k%nx*obj%h%nh*obj%a%nsl))
    else
      if(allocated(obj%en_k_fsm)) deallocate(obj%en_k_fsm)
      allocate(obj%en_k_fsm(0:obj%k%nx*obj%h%nh,2))
      if(allocated(obj%f_k_fsm)) deallocate(obj%f_k_fsm)
      allocate(obj%f_k_fsm(obj%k%nx*obj%h%nh,2))
      if(allocated(obj%indx_fsm)) deallocate(obj%indx_fsm)
      allocate(obj%indx_fsm(0:obj%k%nx*obj%h%nh,2))
    end if

    ! Local band energy
    if(allocated(obj%en_band_local)) deallocate(obj%en_band_local)
    allocate(obj%en_band_local(obj%a%na,obj%e%no_max,obj%a%ns))
    if(allocated(obj%en_band_local_f)) deallocate(obj%en_band_local_f)
    allocate(obj%en_band_local_f(obj%a%na,obj%e%no_max,obj%a%ns))

    obj%en_k = 0.0_rp
    obj%f_k = 0.0_rp
    if(.not. obj%fixed_spin_moment) then
      obj%indx = 0
    else
      obj%en_k_fsm = 0.0_rp
      obj%f_k_fsm = 0.0_rp
      obj%indx_fsm = 0
    end if
    obj%en_f = 0.0_rp
    obj%en_band = 0.0_rp
    obj%en_band_f = 0.0_rp
    obj%en_band_local = 0.0_rp
    obj%en_band_local_f = 0.0_rp
    obj%s_t = 0.0_rp
    obj%en_dc_eei = 0.0_rp
    obj%en_dc_lcn = 0.0_rp
    obj%en_dc_pen = 0.0_rp
    obj%en_in = 0.0_rp
    obj%en_out = 0.0_rp
    obj%delta_en = 0.0_rp
  end subroutine initialize

  subroutine initialize_en_bounds(en_min,en_max)
    real(rp),intent(out) :: en_min,en_max

    en_min = 0.0_rp
    en_max = 0.0_rp
  end subroutine initialize_en_bounds

  subroutine initialize_ffl(fixed_fermi_level,en_f_ffl)
    logical,intent(out) :: fixed_fermi_level
    real(rp),intent(out) :: en_f_ffl

    fixed_fermi_level = .false.
    en_f_ffl = 0.0_rp
  end subroutine initialize_ffl

  subroutine initialize_fsm(fixed_spin_moment,m_fsm)
    logical,intent(out) :: fixed_spin_moment
    real(rp),intent(out) :: m_fsm

    fixed_spin_moment = .false.
    m_fsm = 0.0_rp
  end subroutine initialize_fsm

  subroutine initialize_smearing(smearing,degauss)
    character(len=2), intent(out) :: smearing
    real(rp),intent(out) :: degauss

    smearing = 'mp'
    degauss = 0.0_rp
  end subroutine initialize_smearing

  !> Test if the total energy is converged
  function is_converged(obj,delta_en) result(l)
    class(energy),intent(in) :: obj
    real(rp),intent(in) :: delta_en
    logical :: l

    l = obj%delta_en<delta_en
  end function is_converged

  subroutine print_indx(obj)
    class(energy), intent(in) :: obj
    integer :: i
    logical, save :: to_print=.false.

    if (to_print) then
      to_print = .false.
      do i=0,obj%k%nx*obj%h%nh*obj%a%nsl
        write(*,'(a,i3,a,E22.16,a,i3,a,E22.16)') &
        ' en(',i,')=', obj%en_k(i),&
        ' en(',obj%indx(i),')=', obj%en_k(obj%indx(i))
      end do
    end if
  end subroutine print_indx

  !> Read object in text format from file (default: 'in_energy.txt')
  subroutine read_txt(obj,file)
    class(energy),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    character(len=2) :: smearing
    real(rp) :: degauss
    logical  :: fixed_fermi_level
    real(rp) :: en_f_ffl
    logical  :: fixed_spin_moment
    real(rp) :: m_fsm
    real(rp) :: en_min,en_max
    ! Namelist
    namelist /energy/ smearing, degauss, fixed_fermi_level, en_f_ffl, &
     fixed_spin_moment, m_fsm, en_min, en_max

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_energy.txt'
    end if

    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'energy%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'energy%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    call initialize_smearing(smearing,degauss)
    call initialize_ffl(fixed_fermi_level,en_f_ffl)
    call initialize_fsm(fixed_spin_moment,m_fsm)
    call initialize_en_bounds(en_min,en_max)
    read(10,nml=energy,iostat=iostatus)
    smearing = lower(smearing)
    call check_smearing(trim(smearing))
    call check_fixed_spin_moment(obj%a%ns,fixed_fermi_level,fixed_spin_moment)

    obj%smearing = smearing
    obj%degauss = degauss * obj%u%convert_energy('to','hau')
    obj%fixed_fermi_level = fixed_fermi_level
    obj%en_f_ffl = en_f_ffl * obj%u%convert_energy('to','hau')
    obj%fixed_spin_moment = fixed_spin_moment
    obj%m_fsm = m_fsm
    obj%en_min = en_min * obj%u%convert_energy('to','hau')
    obj%en_max = en_max * obj%u%convert_energy('to','hau')

    call obj%initialize()
    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  subroutine save_en_k(obj,ik,isl,w_k)
    ! INPUT
    class(energy),intent(inout) :: obj
    integer,intent(in) :: ik,isl
    real(rp),dimension(obj%h%nh),intent(in) :: w_k
    ! LOCAL
    integer :: ih,jmat2,jj

    if(obj%fixed_spin_moment) then
      do ih=1,obj%h%nh
        jmat2 = isl+(ih-1)*2
        jj = ik+(jmat2-1)*obj%k%nx
        obj%en_k(jj) = w_k(ih)
      end do

      do ih=1,obj%h%nh
        jj = ik+(ih-1)*obj%k%nx
        obj%en_k_fsm(jj,isl) = w_k(ih)
      end do
    else
      select case(obj%a%ns)
      case(1,2)
        do ih=1,obj%h%nh
          jmat2 = isl+(ih-1)*obj%a%ns
          jj = ik+(jmat2-1)*obj%k%nx
          obj%en_k(jj) = w_k(ih)
        end do
      case(4)
        do ih=1,obj%h%nh
          jj = ik+(ih-1)*obj%k%nx
          obj%en_k(jj) = w_k(ih)
        end do
      end select
      obj%en_k_2(:,ik,isl) = w_k
    end if
  end subroutine save_en_k

  !> Sorting eigenvalues in ascendig order
  subroutine sort_en_k(obj)
    class(energy),intent(inout) :: obj
    integer :: is

    type(group), dimension(:), allocatable :: temp
    integer :: i
    allocate(temp(obj%nen_k))

    if(obj%fixed_spin_moment) then
      do is=1,2
        do i=1,obj%nen_k
          temp(i)%value=obj%en_k_fsm(i,is)
          temp(i)%order=i
        end do
        !call QSort(obj%nen_k,temp,obj%indx_fsm(:,is))
        call indexx(obj%nen_k,obj%en_k_fsm(:,is),obj%indx_fsm(:,is))
      end do
    else
      do i=1,obj%nen_k
        temp(i)%value=obj%en_k(i)
        temp(i)%order=i
      end do
      !call QSort(obj%nen_k,temp,obj%indx)
      call indexx(obj%nen_k,obj%en_k,obj%indx)
      call obj%print_indx()
    endif
    deallocate(temp)
  end subroutine sort_en_k

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_energy.txt')
  subroutine write_txt(obj,file,unit)
    class(energy),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=2) :: smearing
    real(rp) :: degauss
    logical  :: fixed_fermi_level
    real(rp) :: en_f_ffl
    logical  :: fixed_spin_moment
    real(rp) :: m_fsm
    real(rp) :: en_min,en_max
    ! Namelist
    namelist /energy/ smearing, degauss, fixed_fermi_level, en_f_ffl, &
     fixed_spin_moment, m_fsm, en_min, en_max

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_energy.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    smearing = obj%smearing
    degauss = obj%degauss * obj%u%convert_energy('from','hau')
    fixed_fermi_level = obj%fixed_fermi_level
    en_f_ffl = obj%en_f_ffl * obj%u%convert_energy('from','hau')
    fixed_spin_moment = obj%fixed_spin_moment
    m_fsm = obj%m_fsm
    en_min = obj%en_min * obj%u%convert_energy('from','hau')
    en_max = obj%en_max * obj%u%convert_energy('from','hau')

    write(unit_rt,nml=energy)
    call TBKOSTER_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_energy.txt'), if tag (default: .true.) the namelist opening and
  !> closing tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(energy),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    real(rp),dimension(obj%h%nh,obj%k%nx,obj%a%nsl) :: en_k_2
    ! Local variables
    integer :: ip, isl, ik, ih,ia_band1,ia,io,ie
    integer,dimension(:),allocatable         :: iband2io

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_energy.txt'
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
      write(unit_rt,'(a)') '&energy'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('smearing')
        write(unit_rt,'(a)') ' smearing = ''' // trim(obj%smearing) // ''''
      case('degauss')
        write(unit_rt,'(a)') ' degauss = ' // real2str(obj%degauss &
         * obj%u%convert_energy('from','hau'))
      case('fixed_fermi_level')
        write(unit_rt,'(a)') ' fixed_fermi_level = ' &
         // log2str(obj%fixed_fermi_level)
      case('en_f_ffl')
        write(unit_rt,'(a)') ' en_f_ffl = ' // real2str(obj%en_f_ffl &
         * obj%u%convert_energy('from','hau'))
      case('fixed_spin_moment')
        write(unit_rt,'(a)') ' fixed_spin_moment = ' &
         // log2str(obj%fixed_spin_moment)
      case('m_fsm')
        write(unit_rt,'(a)') ' m_fsm = ' // real2str(obj%m_fsm)
      case('en_min')
        write(unit_rt,'(a)') ' en_min = ' // real2str(obj%en_min &
         * obj%u%convert_energy('from','hau'))
      case('en_max')
        write(unit_rt,'(a)') ' en_max = ' // real2str(obj%en_max &
         * obj%u%convert_energy('from','hau'))
      case('en_k_min')
        write(unit_rt,'(a)') ' en_k_min = ' // real2str(minval(obj%en_k) &
         * obj%u%convert_energy('from','hau'))
      case('en_k_max')
        write(unit_rt,'(a)') ' en_k_max = ' // real2str(maxval(obj%en_k) &
         * obj%u%convert_energy('from','hau'))
      case('en_f')
        write(unit_rt,'(a)') ' en_f = ' // real2str(obj%en_f &
         * obj%u%convert_energy('from','hau'))
      case('en_band')
        write(unit_rt,'(a)') ' en_band = ' // real2str(obj%en_band &
         * obj%u%convert_energy('from','hau'))
      case('en_dc_eei')
        write(unit_rt,'(a)') ' en_dc_eei = ' // real2str(obj%en_dc_eei &
         * obj%u%convert_energy('from','hau'))
      case('en_dc_lcn')
        write(unit_rt,'(a)') ' en_dc_lcn = ' // real2str(obj%en_dc_lcn &
         * obj%u%convert_energy('from','hau'))
      case('en_dc_pen')
        write(unit_rt,'(a)') ' en_dc_pen = ' // real2str(obj%en_dc_pen &
         * obj%u%convert_energy('from','hau'))
      case('en_in')
        write(unit_rt,'(a)') ' en = ' // real2str(obj%en_in &
         * obj%u%convert_energy('from','hau'))
      case('en_out')
        write(unit_rt,'(a)') ' en = ' // real2str(obj%en_out &
         * obj%u%convert_energy('from','hau'))
      case('delta_en')
        write(unit_rt,'(a)') ' delta_en = ' // real2str(obj%delta_en &
         * obj%u%convert_energy('from','hau'))
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
end module energy_mod
