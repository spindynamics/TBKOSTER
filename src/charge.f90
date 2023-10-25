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
!  charge.f90
!  TBKOSTER
module charge_mod
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use atom_mod
   use element_mod
   use math_mod, only: i_unit, rad2deg, sqrt_three, cart2sph, sph2cart, nm2rho, rho2nm
   use precision_mod, only: rp
   use string_mod, only: sl, cmplx2str, TBKOSTER_flush, int2str, log2str, lower, real2str
   implicit none
   private

   !> Derived type properties for i/o methods
   character(len=sl), dimension(*), parameter :: property_list = &
                                                 [character(len=sl) :: &
                                                  'q_mul_out', &
                                                  'rho_net_out' &
                                                  ]

   type, public :: charge
      ! Elements
      class(element), pointer :: e
      ! Atoms
      class(atom), pointer :: a

      ! Mulliken charge
      real(rp), dimension(:, :, :), allocatable :: q_mul_in, q_mul_out
      ! Orbital moment
      complex(rp), dimension(:, :, :), allocatable :: om
      ! Mulliken charge difference
      real(rp) :: delta_q_mul
      ! Net charge density
      complex(rp), dimension(:, :, :, :), allocatable :: rho_net_in, rho_net_out
      ! Net charge density difference
      real(rp) :: delta_rho_net
      ! Flag for an orbital-diagonal output net charge density
      logical :: rho_net_out_diagonal

   contains
      ! Destructor
      final :: destructor
      ! Procedures
      procedure :: add_charge_out_k
      procedure :: add_orbital_moment_k
      procedure :: calculate_delta_q
      procedure :: calculate_charge_in
      procedure :: calculate_mulliken_charge_in
      procedure :: calculate_mulliken_charge_out
      procedure :: calculate_charge_out
      procedure :: calculate_net_charge_in
      procedure :: calculate_net_charge_out
      procedure :: calculate_orbital_moment  ! calculate the orbital moment
      procedure :: initialize
      procedure :: is_converged
      procedure :: nullify_charge_out
      procedure :: nullify_orbital_moment
      procedure :: read_txt
      procedure :: set_q_in_to_q_out
      procedure :: write_mulliken_charge_analysis
      procedure :: write_orbital_moment_analysis
      procedure :: write_atom_mag_on_the_fly
      procedure :: write_txt
      procedure :: write_txt_formatted
   end type charge

   interface charge
      procedure :: constructor
   end interface charge

contains
   function constructor(a) result(obj)
      class(atom), target, intent(in) :: a
      type(charge) :: obj

      obj%e => a%e
      obj%a => a

      call obj%initialize()
   end function constructor

   subroutine destructor(obj)
      type(charge) :: obj

      if (allocated(obj%q_mul_in)) deallocate (obj%q_mul_in)
      if (allocated(obj%q_mul_out)) deallocate (obj%q_mul_out)
      if (allocated(obj%rho_net_in)) deallocate (obj%rho_net_in)
      if (allocated(obj%rho_net_out)) deallocate (obj%rho_net_out)
      if (allocated(obj%om)) deallocate (obj%om)
   end subroutine destructor

   subroutine add_charge_out_k(obj, charge_out)
      class(charge), intent(inout) :: obj
      class(charge), intent(in) :: charge_out

      obj%q_mul_out = obj%q_mul_out + charge_out%q_mul_out
      obj%rho_net_out = obj%rho_net_out + charge_out%rho_net_out
   end subroutine add_charge_out_k

   subroutine add_orbital_moment_k(obj, orbital_moment)
      class(charge), intent(inout) :: obj
      class(charge), intent(in) :: orbital_moment

      obj%om = obj%om + orbital_moment%om

   end subroutine add_orbital_moment_k

   subroutine calculate_delta_q(obj)
      class(charge), intent(inout) :: obj

      obj%delta_q_mul = maxval(abs(obj%q_mul_out - obj%q_mul_in))
      obj%delta_rho_net = maxval(abs(obj%rho_net_out - obj%rho_net_in))
   end subroutine calculate_delta_q

   subroutine calculate_charge_in(obj, ia)
      class(charge), intent(inout) :: obj
      integer, dimension(:), intent(in), optional :: ia
      integer, dimension(:), allocatable         :: ia_rt
      integer :: iia

      if (present(ia)) then
         ia_rt = ia
      else
         ia_rt = (/(iia, iia=1, obj%a%na)/)
      end if

      call obj%calculate_mulliken_charge_in(ia_rt)
      call obj%calculate_net_charge_in(ia_rt)

      deallocate (ia_rt)
   end subroutine calculate_charge_in

   subroutine calculate_mulliken_charge_in(obj, ia)
      class(charge), intent(inout) :: obj
      integer, dimension(:), intent(in) :: ia
      integer :: nia, iia, ie, ispin, sigma
      character(len=3) :: os
      real(rp), dimension(3) :: v_sph, v_cart

      nia = size(ia)

      select case (obj%a%ns)
      case (1, 4)
         do iia = 1, nia
            ie = obj%a%ia2ie(ia(iia))
            os = obj%e%os(ie)
            if (os(1:1) == 's') then
               obj%q_mul_in(ia(iia), 1, 0) = obj%e%q_s(ie)
            end if
            if (os(2:2) == 'p') then
               obj%q_mul_in(ia(iia), 2, 0) = obj%e%q_p(ie)
            end if
            if (os(3:3) == 'd') then
               obj%q_mul_in(ia(iia), 3, 0) = obj%e%q_d(ie)
            end if
         end do
      case (2)
         do iia = 1, nia
            ie = obj%a%ia2ie(ia(iia))
            os = obj%e%os(ie)
            do ispin = 1, 2
               sigma = 3 - 2*ispin ! sigma = +-1
               if (os(1:1) == 's') then
                  obj%q_mul_in(ia(iia), 1, ispin - 1) = obj%e%q_s(ie)/2 &
                                                        + sigma*obj%a%m(ia(iia), 1)/(2*10)
               end if
               if (os(2:2) == 'p') then
                  obj%q_mul_in(ia(iia), 2, ispin - 1) = obj%e%q_p(ie)/2 &
                                                        + sigma*obj%a%m(ia(iia), 1)/(2*10)
               end if
               if (os(3:3) == 'd') then
                  obj%q_mul_in(ia(iia), 3, ispin - 1) = obj%e%q_d(ie)/2 &
                                                        + sigma*obj%a%m(ia(iia), 1)/2
               end if
            end do
         end do
      end select

      if (obj%a%ns == 4) then
         do iia = 1, nia
            ie = obj%a%ia2ie(ia(iia))
            os = obj%e%os(ie)
            if (os(1:1) == 's') then
               ! this produce warning in temporary allocated arrays. To avoid, do a local copy
               v_sph(1) = obj%a%m(ia(iia), 1)/10
               v_sph(2) = obj%a%m(ia(iia), 2)
               v_sph(3) = obj%a%m(ia(iia), 3)
               ! obj%q_mul_in(ia(iia),1,1:3) = sph2cart((/obj%a%m(ia(iia),1)/10, &
               ! obj%a%m(ia(iia),2), obj%a%m(ia(iia),3)/))
               v_cart = sph2cart(v_sph)
               obj%q_mul_in(ia(iia), 1, 1:3) = v_cart(1:3)
            end if
            if (os(2:2) == 'p') then
               obj%q_mul_in(ia(iia), 2, 1:3) = sph2cart((/obj%a%m(ia(iia), 1)/10, &
                                                          obj%a%m(ia(iia), 2), obj%a%m(ia(iia), 3)/))
            end if
            if (os(3:3) == 'd') then
               ! this produce warning in temporary allocated arrays. To avoid, do a local copy
               v_sph(:) = obj%a%m(ia(iia), :)
               ! obj%q_mul_in(ia(iia),3,1:3) = sph2cart(obj%a%m(ia(iia),:))
               v_cart = sph2cart(v_sph)
               obj%q_mul_in(ia(iia), 3, 1:3) = v_cart(1:3)
            end if
         end do
      end if
   end subroutine calculate_mulliken_charge_in

   subroutine calculate_mulliken_charge_out(obj, ik, isl, nk, w, iaos2ih, f_k, v_k)
      ! INPUT
      class(charge), intent(inout) :: obj
      integer, intent(in) :: ik, isl
      integer, intent(in) :: nk
      real(rp), intent(in) :: w
      integer, dimension(:, :, :), intent(in) :: iaos2ih
      real(rp), dimension(:), intent(in) :: f_k
      complex(rp), dimension(:, :, :), intent(in) :: v_k
      ! LOCAL
      integer :: ia, ie, io, l, ispin, jspin, imat, jmat2, imat_ispin, imat_jspin, jmat, nn, nh
      complex(rp), dimension(3, 2, 2) :: rho
      complex(rp), dimension(2, 2) :: v_rho
      real(rp) :: n
      real(rp), dimension(3) :: m_cart

      nh = size(v_k, 2)

      obj%q_mul_out = 0.0_rp

      select case (obj%a%ns)
      case (1, 2)
         do ia = 1, obj%a%na
            ie = obj%a%ia2ie(ia)
            do io = 1, obj%e%no(ie)
               imat = iaos2ih(ia, io, 1)
               l = obj%e%o2l(obj%e%o(ie, io))

               do jmat = 1, nh
                  jmat2 = isl + (jmat - 1)*obj%a%ns
                  nn = ik + (jmat2 - 1)*nk
                  obj%q_mul_out(ia, l, isl - 1) = obj%q_mul_out(ia, l, isl - 1) &
                                                  + obj%a%g_s*w*f_k(nn) &
                                                  *real(v_k(1, imat, jmat)*conjg(v_k(2, imat, jmat)))
               end do
            end do
         end do
      case (4)
         do ia = 1, obj%a%na
            ie = obj%a%ia2ie(ia)
            rho = cmplx(0.0_rp, 0.0_rp, kind=rp)
            do ispin = 1, 2
               do jspin = 1, 2
                  do io = 1, obj%e%no(ie)
                     imat_ispin = iaos2ih(ia, io, ispin)
                     imat_jspin = iaos2ih(ia, io, jspin)
                     l = obj%e%o2l(obj%e%o(ie, io))

                     do jmat = 1, nh
                        nn = ik + (jmat - 1)*nk
                        rho(l, ispin, jspin) = rho(l, ispin, jspin) &
                                               + w*f_k(nn) &
                                               *(conjg(v_k(1, imat_jspin, jmat))*v_k(2, imat_ispin, jmat) &
                                                 + v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)))/2
                     end do
                  end do
               end do
            end do

            do l = 1, 3
               v_rho = rho(l, :, :)
               call rho2nm(v_rho, n, m_cart)
               obj%q_mul_out(ia, l, 0) = obj%q_mul_out(ia, l, 0) + n
               obj%q_mul_out(ia, l, 1:3) = obj%q_mul_out(ia, l, 1:3) + m_cart
            end do
         end do
      end select
   end subroutine calculate_mulliken_charge_out

   subroutine calculate_charge_out(obj, ik, isl, nk, w, nh, iaos2ih, f_k, v_k)
      ! INPUT
      class(charge), intent(inout) :: obj
      integer, intent(in) :: ik, isl
      integer, intent(in) :: nk
      real(rp), intent(in) :: w
      integer, intent(in) :: nh
      integer, dimension(obj%a%na, obj%e%no_max, obj%a%ns), intent(in) :: iaos2ih
      real(rp), dimension(:), intent(in) :: f_k
      complex(rp), dimension(:, :, :), intent(in) :: v_k
      call obj%calculate_mulliken_charge_out(ik, isl, nk, w, iaos2ih, f_k, v_k)
      call obj%calculate_net_charge_out(ik, isl, nk, w, iaos2ih, f_k, v_k(1, :, :))
   end subroutine calculate_charge_out

   subroutine calculate_net_charge_in(obj, ia)
      class(charge), intent(inout) :: obj
      integer, dimension(:), intent(in) :: ia
      integer :: nia, iia, ie, io, l, ispin, sigma, isp
      complex(rp), dimension(2, 2) :: rho
      real(rp) :: n
      real(rp), dimension(3) :: m_cart

      nia = size(ia)

      do iia = 1, nia
         ie = obj%a%ia2ie(ia(iia))
         do io = 1, obj%e%no(ie)
            l = obj%e%o2l(obj%e%o(ie, io))
            do isp = 1, obj%a%nsp
               if (l == 1) then
                  obj%rho_net_in(ia(iia), io, io, isp) = obj%e%q_s(ie)/(obj%a%nsp)
               elseif (l == 2) then
                  obj%rho_net_in(ia(iia), io, io, isp) = obj%e%q_p(ie)/(obj%a%nsp*3)
               elseif (l == 3) then
                  obj%rho_net_in(ia(iia), io, io, isp) = obj%e%q_d(ie)/(obj%a%nsp*5)
               end if
            end do
         end do
      end do

      select case (obj%a%ns)
      case (2)
         do iia = 1, nia
            ie = obj%a%ia2ie(ia(iia))
            do io = 1, obj%e%no(ie)
               l = obj%e%o2l(obj%e%o(ie, io))
               do ispin = 1, 2
                  sigma = 3 - 2*ispin ! sigma = +-1
                  if (l == 1) then
                     obj%rho_net_in(ia(iia), io, io, ispin) &
                        = obj%rho_net_in(ia(iia), io, io, ispin) &
                          + sigma*obj%a%m(ia(iia), 1)/(2*10)
                  elseif (l == 2) then
                     obj%rho_net_in(ia(iia), io, io, ispin) &
                        = obj%rho_net_in(ia(iia), io, io, ispin) &
                          + sigma*obj%a%m(ia(iia), 1)/(2*10*3)
                  elseif (l == 3) then
                     obj%rho_net_in(ia(iia), io, io, ispin) &
                        = obj%rho_net_in(ia(iia), io, io, ispin) &
                          + sigma*obj%a%m(ia(iia), 1)/(2*5)
                  end if
               end do
            end do
         end do
      case (4)
         do iia = 1, nia
            ie = obj%a%ia2ie(ia(iia))
            do io = 1, obj%e%no(ie)
               l = obj%e%o2l(obj%e%o(ie, io))
               rho(1, 1) = obj%rho_net_in(ia(iia), io, io, 1)
               rho(2, 2) = obj%rho_net_in(ia(iia), io, io, 2)
               rho(1, 2) = cmplx(0.0_rp, 0.0_rp, kind=rp)
               rho(2, 1) = cmplx(0.0_rp, 0.0_rp, kind=rp)
               call rho2nm(rho, n, m_cart)
               if (l == 1) then
                  m_cart = sph2cart((/obj%a%m(ia(iia), 1)/10, obj%a%m(ia(iia), 2), &
                                      obj%a%m(ia(iia), 3)/))
               elseif (l == 2) then
                  m_cart = sph2cart((/obj%a%m(ia(iia), 1)/(10*3), obj%a%m(ia(iia), 2), &
                                      obj%a%m(ia(iia), 3)/))
               elseif (l == 3) then
                  m_cart = sph2cart((/obj%a%m(ia(iia), 1)/5, obj%a%m(ia(iia), 2), &
                                      obj%a%m(ia(iia), 3)/))
               end if
               call nm2rho(n, m_cart, rho)
               obj%rho_net_in(ia(iia), io, io, 1) = rho(1, 1)
               obj%rho_net_in(ia(iia), io, io, 2) = rho(2, 2)
               obj%rho_net_in(ia(iia), io, io, 3) = rho(1, 2)
               obj%rho_net_in(ia(iia), io, io, 4) = rho(2, 1)
            end do
         end do
      end select
   end subroutine calculate_net_charge_in

   subroutine calculate_net_charge_out(obj, ik, isl, nk, w, iaos2ih, f_k, v_k)
      ! INPUT
      class(charge), intent(inout) :: obj
      integer, intent(in) :: ik, isl
      integer, intent(in) :: nk
      real(rp), intent(in) :: w
      integer, dimension(:, :, :), intent(in) :: iaos2ih
      real(rp), dimension(:), intent(in) :: f_k
      complex(rp), dimension(:, :), intent(in) :: v_k
      ! LOCAL
      integer :: ia, ie, io1, io2, ispin, jspin, imat, jmat, kmat, kmat2, &
                 imat_ispin, jmat_jspin, nn, nh

      nh = size(v_k, 1)

      obj%rho_net_out = cmplx(0.0_rp, 0.0_rp, kind=rp)

      select case (obj%a%ns)
      case (1, 2)
         if (obj%rho_net_out_diagonal) then
            do ia = 1, obj%a%na
               ie = obj%a%ia2ie(ia)
               do io1 = 1, obj%e%no(ie)
                  imat = iaos2ih(ia, io1, 1)
                  jmat = imat
                  do kmat = 1, nh
                     kmat2 = isl + (kmat - 1)*obj%a%ns
                     nn = ik + (kmat2 - 1)*nk
                     obj%rho_net_out(ia, io1, io1, isl) &
                        = obj%rho_net_out(ia, io1, io1, isl) &
                          + obj%a%g_s*w*f_k(nn)*v_k(imat, kmat)*conjg(v_k(jmat, kmat))
                  end do
               end do
            end do
         else !if(.not. obj%rho_net_out_diagonal)
            do ia = 1, obj%a%na
               ie = obj%a%ia2ie(ia)
               do io1 = 1, obj%e%no(ie)
                  imat = iaos2ih(ia, io1, 1)
                  do io2 = 1, obj%e%no(ie)
                     jmat = iaos2ih(ia, io2, 1)
                     do kmat = 1, nh
                        kmat2 = isl + (kmat - 1)*obj%a%ns
                        nn = ik + (kmat2 - 1)*nk
                        obj%rho_net_out(ia, io1, io2, isl) &
                           = obj%rho_net_out(ia, io1, io2, isl) &
                             + obj%a%g_s*w*f_k(nn)*v_k(imat, kmat)*conjg(v_k(jmat, kmat))
                     end do
                  end do
               end do
            end do
         end if
      case (4)
         if (obj%rho_net_out_diagonal) then
            do ia = 1, obj%a%na
               ie = obj%a%ia2ie(ia)
               do io1 = 1, obj%e%no(ie)
                  do ispin = 1, 2
                     imat_ispin = iaos2ih(ia, io1, ispin)
                     do jspin = 1, 2
                        jmat_jspin = iaos2ih(ia, io1, jspin)
                        do jmat = 1, nh
                           nn = ik + (jmat - 1)*nk
                           obj%rho_net_out(ia, io1, io1, obj%a%iss2is(ispin, jspin)) &
                              = obj%rho_net_out(ia, io1, io1, obj%a%iss2is(ispin, jspin)) &
                                + f_k(nn)*w*v_k(imat_ispin, jmat)*conjg(v_k(jmat_jspin, jmat))
                        end do
                     end do
                  end do
               end do
            end do
         else !if(.not. obj%rho_net_out_diagonal)
            do ia = 1, obj%a%na
               ie = obj%a%ia2ie(ia)
               do io1 = 1, obj%e%no(ie)
                  do io2 = 1, obj%e%no(ie)
                     do ispin = 1, 2
                        imat_ispin = iaos2ih(ia, io1, ispin)
                        do jspin = 1, 2
                           jmat_jspin = iaos2ih(ia, io2, jspin)
                           do jmat = 1, nh
                              nn = ik + (jmat - 1)*nk
                              obj%rho_net_out(ia, io1, io2, obj%a%iss2is(ispin, jspin)) &
                                 = obj%rho_net_out(ia, io1, io2, obj%a%iss2is(ispin, jspin)) &
                                   + f_k(nn)*w*v_k(imat_ispin, jmat)*conjg(v_k(jmat_jspin, jmat))
                           end do
                        end do
                     end do
                  end do
               end do
            end do
         end if
      end select
   end subroutine calculate_net_charge_out

   subroutine initialize(obj)
      class(charge), intent(inout) :: obj

      if (allocated(obj%q_mul_in)) deallocate (obj%q_mul_in)
      allocate (obj%q_mul_in(obj%a%na, 3, 0:obj%a%ns - 1))
      if (allocated(obj%q_mul_out)) deallocate (obj%q_mul_out)
      allocate (obj%q_mul_out(obj%a%na, 3, 0:obj%a%ns - 1))
      if (allocated(obj%rho_net_in)) deallocate (obj%rho_net_in)
      allocate (obj%rho_net_in(obj%a%na, obj%e%no_max, obj%e%no_max, obj%a%ns))
      if (allocated(obj%rho_net_out)) deallocate (obj%rho_net_out)
      allocate (obj%rho_net_out(obj%a%na, obj%e%no_max, obj%e%no_max, obj%a%ns))
      if (allocated(obj%om)) deallocate (obj%om)
      ! if(obj%a%ns==4) allocate(obj%om(obj%a%na,3,3))
      allocate (obj%om(obj%a%na, 3, 3))

      obj%q_mul_in = 0.0_rp
      obj%q_mul_out = 0.0_rp
      obj%delta_q_mul = 0.0_rp
      obj%rho_net_in = cmplx(0.0_rp, 0.0_rp, kind=rp)
      obj%rho_net_out = cmplx(0.0_rp, 0.0_rp, kind=rp)
      obj%delta_rho_net = 0.0_rp
      obj%rho_net_out_diagonal = .true.
      obj%om = cmplx(0.0_rp, 0.0_rp, kind=rp)
   end subroutine initialize

   !> Test if the Mulliken charge and the net charge density are converged
   function is_converged(obj, delta_q) result(l)
      class(charge), intent(in) :: obj
      real(rp), intent(in) :: delta_q
      logical :: l

      l = obj%delta_q_mul < delta_q .and. obj%delta_rho_net < delta_q
   end function is_converged

   subroutine nullify_charge_out(obj)
      class(charge), intent(inout) :: obj

      obj%q_mul_out = 0.0_rp
      obj%rho_net_out = cmplx(0.0_rp, 0.0_rp, kind=rp)
   end subroutine nullify_charge_out

   subroutine nullify_orbital_moment(obj)
      class(charge), intent(inout) :: obj

      obj%om = cmplx(0.0_rp, 0.0_rp, kind=rp)
   end subroutine nullify_orbital_moment
   
   !> Read object in text format from file (default: 'in_charge.txt')
   subroutine read_txt(obj, file)
      class(charge), intent(inout) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable :: file_rt
      integer :: iostatus
      logical :: isopen
      ! Namelist variables
      real(rp), dimension(obj%a%na, 3, 0:obj%a%ns - 1) :: q_mul
      complex(rp), dimension(obj%a%na, obj%e%no_max, obj%e%no_max, obj%a%ns) &
         :: rho_net
      ! Namelist
      namelist /charge/ q_mul, rho_net

      if (present(file)) then
         file_rt = trim(file)
      else
         file_rt = 'in_charge.txt'
      end if

      inquire (unit=10, opened=isopen)
      if (isopen) then
         write (error_unit, '(a)') 'charge%read_txt() : Unit 10 is already open'
         error stop
      else
         open (unit=10, file=file_rt, action='read', iostat=iostatus, status='old')
      end if
      if (iostatus /= 0) then
         write (error_unit, *) 'charge%read_txt(): file ', file_rt, ' not found'
         error stop
      end if

      read (10, nml=charge)

      obj%q_mul_in = q_mul
      obj%rho_net_in = rho_net

      close (unit=10)
      !deallocate(file_rt)
   end subroutine read_txt

   subroutine set_q_in_to_q_out(obj)
      class(charge), intent(inout) :: obj

      obj%q_mul_in = obj%q_mul_out
      obj%rho_net_in = obj%rho_net_out
   end subroutine set_q_in_to_q_out

   subroutine write_mulliken_charge_analysis(obj, unit, intent)
      class(charge), intent(in) :: obj
      integer, intent(in), optional :: unit
      integer                     :: unit_rt
      character(len=*), intent(in), optional :: intent
      character(len=:), allocatable         :: intent_rt
      character(len=:), allocatable :: line, header, fmt1, fmt2, fmt3
      integer :: ia
      real(rp) :: n_s, n_p, n_d, n
      real(rp) :: n_tot
      real(rp) :: m_s, m_p, m_d, m
      real(rp), dimension(3) :: m_cart, m_cart_s, m_cart_p, m_cart_d
      real(rp), dimension(3) :: m_sph, m_sph_s, m_sph_p, m_sph_d
      real(rp), dimension(3) :: m_tot, v_sph
      real(rp) :: m_r_tot

      !write(output_unit,*) "====> Entering write_mulliken_charge_analysis"

      if (present(unit)) then
         unit_rt = unit
      else
         unit_rt = output_unit
      end if
      if (present(intent)) then
         intent_rt = intent
      else
         intent_rt = 'out'
      end if

      select case (obj%a%ns)
      case (1, 2)
         line = '|'//repeat('-', 66)//'|'
         header = '|'//repeat(' ', 6)//'|' &
                  //repeat(' ', 7)//'s'//repeat(' ', 6)//'|' &
                  //repeat(' ', 7)//'p'//repeat(' ', 6)//'|' &
                  //repeat(' ', 7)//'d'//repeat(' ', 6)//'|' &
                  //repeat(' ', 5)//'total'//repeat(' ', 4)//'|'
      case (4)
         line = '|'//repeat('-', 72)//'|'
         header = '|'//repeat(' ', 12)//'|' &
                  //repeat(' ', 7)//'s'//repeat(' ', 6)//'|' &
                  //repeat(' ', 7)//'p'//repeat(' ', 6)//'|' &
                  //repeat(' ', 7)//'d'//repeat(' ', 6)//'|' &
                  //repeat(' ', 5)//'total'//repeat(' ', 4)//'|'
      end select
      fmt1 = '(a,f12.7,a,f12.7,a,f12.7,a,f12.7,a)'
      fmt2 = '(a,f12.7,a)'
      fmt3 = '(a,f12.7,a,f12.7,a,f12.7,a)'

      write (unit_rt, '(a)') line
      write (unit_rt, '(a)') header
      write (unit_rt, '(a)') line

      n_tot = 0.0_rp
      m_tot = 0.0_rp
      m_r_tot = 0.0_rp

      select case (obj%a%ns)
      case (1)
         select case (intent)
         case ('in')
            do ia = 1, obj%a%na
               n_s = obj%q_mul_in(ia, 1, 0)
               n_p = obj%q_mul_in(ia, 2, 0)
               n_d = obj%q_mul_in(ia, 3, 0)
               n = n_s + n_p + n_d
               n_tot = n_tot + n

               write (unit_rt, fmt1) '| n('//int2str(ia)//') | ', &
                  n_s, ' | ', n_p, ' | ', n_d, ' | ', n, ' |'
               write (unit_rt, '(a)') line
            end do
         case ('out')
            do ia = 1, obj%a%na
               n_s = obj%q_mul_out(ia, 1, 0)
               n_p = obj%q_mul_out(ia, 2, 0)
               n_d = obj%q_mul_out(ia, 3, 0)
               n = n_s + n_p + n_d
               n_tot = n_tot + n

               write (unit_rt, fmt1) '| n('//int2str(ia)//') | ', &
                  n_s, ' | ', n_p, ' | ', n_d, ' | ', n, ' |'
               write (unit_rt, '(a)') line
            end do
         end select
         write (unit_rt, fmt2) '| n_tot = ', n_tot, repeat(' ', 45)//'|'
         write (unit_rt, '(a)') line
      case (2)
         select case (intent)
         case ('in')
            do ia = 1, obj%a%na
               n_s = obj%q_mul_in(ia, 1, 0) + obj%q_mul_in(ia, 1, 1)
               n_p = obj%q_mul_in(ia, 2, 0) + obj%q_mul_in(ia, 2, 1)
               n_d = obj%q_mul_in(ia, 3, 0) + obj%q_mul_in(ia, 3, 1)
               n = n_s + n_p + n_d
               n_tot = n_tot + n

               m_s = obj%q_mul_in(ia, 1, 0) - obj%q_mul_in(ia, 1, 1)
               m_p = obj%q_mul_in(ia, 2, 0) - obj%q_mul_in(ia, 2, 1)
               m_d = obj%q_mul_in(ia, 3, 0) - obj%q_mul_in(ia, 3, 1)
               m = m_s + m_p + m_d
               m_tot(1) = m_tot(1) + m
               m_r_tot = m_r_tot + abs(m_s) + abs(m_p) + abs(m_d)

               write (unit_rt, fmt1) '| n('//int2str(ia)//') | ', &
                  n_s, ' | ', n_p, ' | ', n_d, ' | ', n, ' |'
               write (unit_rt, fmt1) '| m('//int2str(ia)//') | ', &
                  m_s, ' | ', m_p, ' | ', m_d, ' | ', m, ' |'
               write (unit_rt, '(a)') line
            end do
         case ('out')
            do ia = 1, obj%a%na
               n_s = obj%q_mul_out(ia, 1, 0) + obj%q_mul_out(ia, 1, 1)
               n_p = obj%q_mul_out(ia, 2, 0) + obj%q_mul_out(ia, 2, 1)
               n_d = obj%q_mul_out(ia, 3, 0) + obj%q_mul_out(ia, 3, 1)
               n = n_s + n_p + n_d
               n_tot = n_tot + n

               m_s = obj%q_mul_out(ia, 1, 0) - obj%q_mul_out(ia, 1, 1)
               m_p = obj%q_mul_out(ia, 2, 0) - obj%q_mul_out(ia, 2, 1)
               m_d = obj%q_mul_out(ia, 3, 0) - obj%q_mul_out(ia, 3, 1)
               m = m_s + m_p + m_d
               !          obj%a%m(ia,:)=m  ! update m
               m_tot(1) = m_tot(1) + m
               m_r_tot = m_r_tot + abs(m_s) + abs(m_p) + abs(m_d)

               write (unit_rt, fmt1) '| n('//int2str(ia)//') | ', &
                  n_s, ' | ', n_p, ' | ', n_d, ' | ', n, ' |'
               write (unit_rt, fmt1) '| m('//int2str(ia)//') | ', &
                  m_s, ' | ', m_p, ' | ', m_d, ' | ', m, ' |'
               write (unit_rt, '(a)') line
            end do
         end select
         write (unit_rt, fmt2) '| n_tot   = ', n_tot, repeat(' ', 43)//'|'
         write (unit_rt, fmt2) '| m_tot   = ', m_tot(1), repeat(' ', 43)//'|'
         write (unit_rt, fmt2) '| m_r_tot = ', m_r_tot, repeat(' ', 43)//'|'
         write (unit_rt, '(a)') line
      case (4)
         select case (intent)
         case ('in')
            do ia = 1, obj%a%na
               n_s = obj%q_mul_in(ia, 1, 0)
               n_p = obj%q_mul_in(ia, 2, 0)
               n_d = obj%q_mul_in(ia, 3, 0)
               n = n_s + n_p + n_d
               n_tot = n_tot + n

               !m_cart = sum(obj%q_mul_in(ia,:,1:3),1)
               m_cart_s = obj%q_mul_in(ia, 1, 1:3)
               m_cart_p = obj%q_mul_in(ia, 2, 1:3)
               m_cart_d = obj%q_mul_in(ia, 3, 1:3)
               m_cart = m_cart_s + m_cart_p + m_cart_d
               !write(output_unit,*) "===>"
               !write(output_unit,*) "m_cart_s =", m_cart_s
               !write(output_unit,*) "m_cart_p =", m_cart_p
               !write(output_unit,*) "m_cart_d =", m_cart_d
               !write(output_unit,*) "m_cart   =", m_cart
               m_sph = cart2sph(m_cart)
               !write(output_unit,*) "===>"
               m_sph_s = cart2sph(m_cart_s)
               m_sph_p = cart2sph(m_cart_p)
               m_sph_d = cart2sph(m_cart_d)
               m_tot = m_tot + m_cart
               m_r_tot = m_r_tot + norm2(m_cart)

               write (unit_rt, fmt1) '|    n('//int2str(ia)//')    | ', &
                  n_s, ' | ', n_p, ' | ', n_d, ' | ', n, ' |'
               write (unit_rt, fmt1) '|   m_r('//int2str(ia)//')   | ', &
                  m_sph_s(1), ' | ', m_sph_p(1), ' | ', m_sph_d(1), ' | ', &
                  m_sph(1), ' |'
               write (unit_rt, fmt1) '| m_theta('//int2str(ia)//') | ', &
                  m_sph_s(2), ' | ', m_sph_p(2), ' | ', m_sph_d(2), ' | ', &
                  m_sph(2), ' |'
               write (unit_rt, fmt1) '|  m_phi('//int2str(ia)//')  | ', &
                  m_sph_s(3), ' | ', m_sph_p(3), ' | ', m_sph_d(3), ' | ', &
                  m_sph(3), ' |'
               write (unit_rt, fmt1) '|   m_x('//int2str(ia)//')   | ', &
                  m_cart_s(1), ' | ', m_cart_p(1), ' | ', m_cart_d(1), ' | ', &
                  m_cart(1), ' |'
               write (unit_rt, fmt1) '|   m_y('//int2str(ia)//')   | ', &
                  m_cart_s(2), ' | ', m_cart_p(2), ' | ', m_cart_d(2), ' | ', &
                  m_cart(2), ' |'
               write (unit_rt, fmt1) '|   m_z('//int2str(ia)//')   | ', &
                  m_cart_s(3), ' | ', m_cart_p(3), ' | ', m_cart_d(3), ' | ', &
                  m_cart(3), ' |'
               write (unit_rt, '(a)') line
            end do
         case ('out')
            do ia = 1, obj%a%na
               n_s = obj%q_mul_out(ia, 1, 0)
               n_p = obj%q_mul_out(ia, 2, 0)
               n_d = obj%q_mul_out(ia, 3, 0)
               n = n_s + n_p + n_d
               n_tot = n_tot + n

               m_cart_s = obj%q_mul_out(ia, 1, 1:3)
               m_cart_p = obj%q_mul_out(ia, 2, 1:3)
               m_cart_d = obj%q_mul_out(ia, 3, 1:3)
               m_cart = m_cart_s + m_cart_p + m_cart_d
               m_sph = cart2sph(m_cart)
               !          obj%a%m(ia,:)=m_cart  ! update m
               m_sph_s = cart2sph(m_cart_s)
               m_sph_p = cart2sph(m_cart_p)
               m_sph_d = cart2sph(m_cart_d)
               m_tot = m_tot + m_cart
               m_r_tot = m_r_tot + norm2(m_cart)

               write (unit_rt, fmt1) '|    n('//int2str(ia)//')    | ', &
                  n_s, ' | ', n_p, ' | ', n_d, ' | ', n, ' |'
               write (unit_rt, fmt1) '|   m_r('//int2str(ia)//')   | ', &
                  m_sph_s(1), ' | ', m_sph_p(1), ' | ', m_sph_d(1), ' | ', &
                  m_sph(1), ' |'
               write (unit_rt, fmt1) '| m_theta('//int2str(ia)//') | ', &
                  m_sph_s(2)*rad2deg, ' | ', m_sph_p(2)*rad2deg, ' | ', &
                  m_sph_d(2)*rad2deg, ' | ', m_sph(2)*rad2deg, ' |'
               write (unit_rt, fmt1) '|  m_phi('//int2str(ia)//')  | ', &
                  m_sph_s(3)*rad2deg, ' | ', m_sph_p(3)*rad2deg, ' | ', &
                  m_sph_d(3)*rad2deg, ' | ', m_sph(3)*rad2deg, ' |'
               write (unit_rt, fmt1) '|   m_x('//int2str(ia)//')   | ', &
                  m_cart_s(1), ' | ', m_cart_p(1), ' | ', m_cart_d(1), ' | ', &
                  m_cart(1), ' |'
               write (unit_rt, fmt1) '|   m_y('//int2str(ia)//')   | ', &
                  m_cart_s(2), ' | ', m_cart_p(2), ' | ', m_cart_d(2), ' | ', &
                  m_cart(2), ' |'
               write (unit_rt, fmt1) '|   m_z('//int2str(ia)//')   | ', &
                  m_cart_s(3), ' | ', m_cart_p(3), ' | ', m_cart_d(3), ' | ', &
                  m_cart(3), ' |'
               write (unit_rt, '(a)') line
            end do
         end select
         write (unit_rt, fmt2) '| n_tot   = ', n_tot, repeat(' ', 49)//'|'
         write (unit_rt, fmt3) '| m_tot   = ', &
            m_tot(1), ', ', &
            m_tot(2), ', ', &
            m_tot(3), repeat(' ', 21)//'|'
         write (unit_rt, fmt2) '| m_r_tot = ', m_r_tot, repeat(' ', 49)//'|'
         write (unit_rt, '(a)') line
      end select
      call TBKOSTER_flush(unit_rt)
   end subroutine write_mulliken_charge_analysis

   subroutine write_orbital_moment_analysis(obj, unit)
      class(charge), intent(in) :: obj
      integer, intent(in), optional :: unit
      integer                     :: unit_rt
      character(len=:), allocatable :: line, header, fmt1, fmt2, fmt3
      integer :: ia
      real(rp), dimension(3) :: om_cart, om_cart_s, om_cart_p, om_cart_d
      real(rp), dimension(3) :: om_sph, om_sph_s, om_sph_p, om_sph_d
      real(rp), dimension(3) :: om_tot, v_sph
      real(rp) :: om_r_tot

      !write(output_unit,*) "====> Entering write_mulliken_charge_analysis"

      if (present(unit)) then
         unit_rt = unit
      else
         unit_rt = output_unit
      end if

      line = '|'//repeat('-', 72)//'|'
      header = '|'//repeat(' ', 12)//'|' &
               //repeat(' ', 7)//'s'//repeat(' ', 6)//'|' &
               //repeat(' ', 7)//'p'//repeat(' ', 6)//'|' &
               //repeat(' ', 7)//'d'//repeat(' ', 6)//'|' &
               //repeat(' ', 5)//'total'//repeat(' ', 4)//'|'

      fmt1 = '(a,f12.7,a,f12.7,a,f12.7,a,f12.7,a)'
      fmt2 = '(a,f12.7,a)'
      fmt3 = '(a,f12.7,a,f12.7,a,f12.7,a)'

      write (unit_rt, '(a)') line
      write (unit_rt, '(a)') '           Orbital moment analysis'
      write (unit_rt, '(a)') line
      write (unit_rt, '(a)') header
      write (unit_rt, '(a)') line

      om_tot = 0.0_rp
      om_r_tot = 0.0_rp

      do ia = 1, obj%a%na
         om_cart_s = obj%om(ia, 1:3, 1)
         om_cart_p = obj%om(ia, 1:3, 2)
         om_cart_d = obj%om(ia, 1:3, 3)
         om_cart = om_cart_s + om_cart_p + om_cart_d
         om_sph = cart2sph(om_cart)
         !          obj%a%m(ia,:)=m_cart  ! update m
         om_sph_s = cart2sph(om_cart_s)
         om_sph_p = cart2sph(om_cart_p)
         om_sph_d = cart2sph(om_cart_d)
         om_tot = om_tot + om_cart
         om_r_tot = om_r_tot + norm2(om_cart)

         write (unit_rt, fmt1) '|   om_r('//int2str(ia)//')   | ', &
            om_sph_s(1), ' | ', om_sph_p(1), ' | ', om_sph_d(1), ' | ', &
            om_sph(1), ' |'
         write (unit_rt, fmt1) '| om_theta('//int2str(ia)//') | ', &
            om_sph_s(2)*rad2deg, ' | ', om_sph_p(2)*rad2deg, ' | ', &
            om_sph_d(2)*rad2deg, ' | ', om_sph(2)*rad2deg, ' |'
         write (unit_rt, fmt1) '|  om_phi('//int2str(ia)//')  | ', &
            om_sph_s(3)*rad2deg, ' | ', om_sph_p(3)*rad2deg, ' | ', &
            om_sph_d(3)*rad2deg, ' | ', om_sph(3)*rad2deg, ' |'
         write (unit_rt, fmt1) '|   om_x('//int2str(ia)//')   | ', &
            om_cart_s(1), ' | ', om_cart_p(1), ' | ', om_cart_d(1), ' | ', &
            om_cart(1), ' |'
         write (unit_rt, fmt1) '|   om_y('//int2str(ia)//')   | ', &
            om_cart_s(2), ' | ', om_cart_p(2), ' | ', om_cart_d(2), ' | ', &
            om_cart(2), ' |'
         write (unit_rt, fmt1) '|   om_z('//int2str(ia)//')   | ', &
            om_cart_s(3), ' | ', om_cart_p(3), ' | ', om_cart_d(3), ' | ', &
            om_cart(3), ' |'
         write (unit_rt, '(a)') line
      end do

      write (unit_rt, fmt3) '| om_tot   = ', &
         om_tot(1), ', ', &
         om_tot(2), ', ', &
         om_tot(3), repeat(' ', 21)//'|'
      write (unit_rt, fmt2) '| om_r_tot = ', om_r_tot, repeat(' ', 49)//'|'
      write (unit_rt, '(a)') line

      call TBKOSTER_flush(unit_rt)
   end subroutine write_orbital_moment_analysis

   subroutine write_atom_mag_on_the_fly(obj, file, intent)
      !    class(charge),intent(in) :: obj
      class(charge) :: obj
      integer, parameter                   :: unit_rt = 10
      character(len=*), intent(in), optional :: intent
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable         :: file_rt
      real(rp) :: m_s, m_p, m_d
      real(rp), dimension(3) :: m_cart, m_sph
      real(rp), dimension(:, :), allocatable :: m
      ! Namelist variables
      real(rp), dimension(obj%a%na, 3) :: r
      real(rp), dimension(obj%a%na, obj%a%nn_max, 3) :: rn
      real(rp), dimension(obj%a%na) :: lambda_pen
      ! Local variables
      integer :: ia, in, ip, itag
      character(len=sl), dimension(*), parameter :: property_rt = &
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
      file_rt = file
      open (unit=unit_rt, file=file_rt, action='write')
      allocate (m(obj%a%na, 3))
      m_cart = 0_rp
      m_sph = 0_rp
      m = 0_rp
      select case (obj%a%ns)
      case (1)
         m = 0_rp
      case (2)
         do ia = 1, obj%a%na
            m_s = obj%q_mul_out(ia, 1, 0) - obj%q_mul_out(ia, 1, 1)
            m_p = obj%q_mul_out(ia, 2, 0) - obj%q_mul_out(ia, 2, 1)
            m_d = obj%q_mul_out(ia, 3, 0) - obj%q_mul_out(ia, 3, 1)
            m(ia, 1) = m_s + m_p + m_d
         end do
      case (4)
         do ia = 1, obj%a%na
            m_cart = sum(obj%q_mul_out(ia, :, 1:3), 1)
            m_sph = cart2sph(m_cart)
            m(ia, :) = m_sph
         end do

      end select
      if (intent == 'out') then
         obj%a%m = m
      end if
      write (unit_rt, '(a)') '&atom'

      do ip = 1, size(property_rt)
         select case (lower(trim(property_rt(ip))))
         case ('ns')
            write (unit_rt, '(a)') ' ns = '//int2str(obj%a%ns)
         case ('na')
            write (unit_rt, '(a)') ' na = '//int2str(obj%a%na)
         case ('ntag')
            write (unit_rt, '(a)') ' ntag = '//int2str(obj%a%ntag)
         case ('stag')
            do itag = 1, obj%a%ntag
               write (unit_rt, '(a)') ' stag('//int2str(itag)//') = ' &
                  //int2str(obj%a%stag(itag))
            end do
         case ('tag')
            do itag = 1, obj%a%ntag
               write (unit_rt, '(a)') ' tag('//int2str(itag)//') = ''' &
                  //trim(obj%a%tag(itag))//''''
            end do
         case ('ia2ie')
            do ia = 1, obj%a%na
               write (unit_rt, '(a)') ' ia2ie('//int2str(ia)//') = ' &
                  //int2str(obj%a%ia2ie(ia))
            end do
         case ('nel')
            write (unit_rt, '(a)') ' nel = '//real2str(obj%a%nel)
         case ('pbc')
            write (unit_rt, '(a)') ' pbc = ' &
               //int2str(obj%a%pbc(1))//', ' &
               //int2str(obj%a%pbc(2))//', ' &
               //int2str(obj%a%pbc(3))
         case ('k_spiral')
            write (unit_rt, '(a)') ' k_spiral = ' &
               //real2str(obj%a%k_spiral(1))//', ' &
               //real2str(obj%a%k_spiral(2))//', ' &
               //real2str(obj%a%k_spiral(3))
         case ('r_coord')
            write (unit_rt, '(a)') ' r_coord = '''//obj%a%r_coord//''''
         case ('m_coord')
            write (unit_rt, '(a)') ' m_coord = '''//obj%a%m_coord//''''
         case ('r')
            do ia = 1, obj%a%na
               if (obj%a%r_coord == 'cartesian') then
                  r = obj%a%r*obj%a%u%convert_length('from', 'hau')
               else
                  r = obj%a%r
               end if
               write (unit_rt, '(a)') ' r('//int2str(ia)//',:) = ' &
                  //real2str(r(ia, 1))//', ' &
                  //real2str(r(ia, 2))//', ' &
                  //real2str(r(ia, 3))
            end do
         case ('p')
            do ia = 1, obj%a%na
               write (unit_rt, '(a)') ' p('//int2str(ia)//',:) = ' &
                  //real2str(obj%a%p(ia, 1))//', ' &
                  //real2str(obj%a%p(ia, 2))//', ' &
                  //real2str(obj%a%p(ia, 3))
            end do
         case ('m')
            do ia = 1, obj%a%na
               write (unit_rt, '(a)') ' m('//int2str(ia)//',:) = ' &
                  //real2str(obj%a%m(ia, 1))//', ' &
                  //real2str(obj%a%m(ia, 2)*rad2deg)//', ' &
                  //real2str(obj%a%m(ia, 3)*rad2deg)
            end do
         case ('nn')
            do ia = 1, obj%a%na
               write (unit_rt, '(a)') ' nn('//int2str(ia)//') = ' &
                  //int2str(obj%a%nn(ia))
            end do
         case ('nn_max')
            write (unit_rt, '(a)') ' nn_max = '//int2str(obj%a%nn_max)
         case ('ian2ia')
            do ia = 1, obj%a%na
               do in = 1, obj%a%nn(ia)
                  write (unit_rt, '(a)') ' ian2ia('//int2str(ia)//',' &
                     //int2str(in)//') = '//int2str(obj%a%ian2ia(ia, in))
               end do
            end do
         case ('rn')
            rn = obj%a%rn*obj%a%u%convert_length('from', 'hau')
            do ia = 1, obj%a%na
               do in = 1, obj%a%nn(ia)
                  write (unit_rt, '(a)') ' rn('//int2str(ia)//','//int2str(in) &
                     //',:) = ' &
                     //real2str(rn(ia, in, 1))//', ' &
                     //real2str(rn(ia, in, 2))//', ' &
                     //real2str(rn(ia, in, 3))
               end do
            end do
         case ('lambda_pen')
            lambda_pen = obj%a%lambda_pen*obj%a%u%convert_energy('from', 'hau')
            do ia = 1, obj%a%na
               write (unit_rt, '(a)') ' lambda_pen('//int2str(ia)//') = ' &
                  //real2str(lambda_pen(ia))
            end do
         case ('b_pen')
            do ia = 1, obj%a%na
               write (unit_rt, '(a)') ' b_pen('//int2str(ia)//',:) = ' &
                  //real2str(obj%a%b_pen(ia, 1))//', ' &
                  //real2str(obj%a%b_pen(ia, 2))//', ' &
                  //real2str(obj%a%b_pen(ia, 3))
            end do
         end select
      end do

      write (unit_rt, '(a)') ' /'
      call TBKOSTER_flush(unit_rt)
      close (unit_rt)

   end subroutine write_atom_mag_on_the_fly

   !> Write object in text format to unit (default: 10), if it's a file
   !> its name is set to file (default: 'out_charge.txt')
   subroutine write_txt(obj, file, unit)
      class(charge), intent(in) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable         :: file_rt
      integer, intent(in), optional :: unit
      integer                     :: unit_rt
      ! Namelist variables
      real(rp), dimension(obj%a%na, 3, 0:obj%a%ns - 1) :: q_mul
      complex(rp), dimension(obj%a%na, obj%e%no_max, obj%e%no_max, obj%a%ns) &
         :: rho_net
      ! Namelist
      namelist /charge/ q_mul, rho_net

      if (present(file)) then
         file_rt = file
      else
         file_rt = 'out_charge.txt'
      end if
      if (present(unit)) then
         unit_rt = unit
      else
         unit_rt = 10
      end if

      if (.not. present(unit)) then
         open (unit=unit_rt, file=file_rt, action='write')
      end if

      q_mul = obj%q_mul_out
      rho_net = obj%rho_net_out

      write (unit_rt, nml=charge)
      call TBKOSTER_flush(unit_rt)

      if (.not. present(unit)) then
         close (unit_rt)
      end if
      !deallocate(file_rt)
   end subroutine write_txt

   !> Write property (default: property_list) in text format to unit
   !> (default: 10), if it's a file its name is set to file (default:
   !> 'out_charge.txt'), if tag (default: .true.) the namelist opening and
   !> closing tags are written
   subroutine write_txt_formatted(obj, file, property, tag, unit)
      class(charge), intent(in) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable         :: file_rt
      character(len=*), dimension(:), intent(in), optional :: property
      character(len=:), dimension(:), allocatable         :: property_rt
      logical, intent(in), optional :: tag
      logical                     :: tag_rt
      integer, intent(in), optional :: unit
      integer                     :: unit_rt
      ! Local variables
      integer :: ia, io1, io2, ip, is, l

      if (present(file)) then
         file_rt = file
      else
         file_rt = 'out_charge.txt'
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

      if (.not. present(unit)) then
         open (unit=unit_rt, file=file_rt, action='write')
      end if
      if (tag_rt) then
         write (unit_rt, '(a)') '&charge'
      end if

      do ip = 1, size(property_rt)
         select case (lower(trim(property_rt(ip))))
         case ('q_mul_in')
            do ia = 1, obj%a%na
               do l = 1, 3
                  do is = 0, obj%a%ns - 1
                     write (unit_rt, '(a)') ' q_mul('//int2str(ia)//',' &
                        //int2str(l)//','//int2str(is)//') = ' &
                        //real2str(obj%q_mul_in(ia, l, is))
                  end do
               end do
            end do
         case ('q_mul_out')
            do ia = 1, obj%a%na
               do l = 1, 3
                  do is = 0, obj%a%ns - 1
                     write (unit_rt, '(a)') ' q_mul('//int2str(ia)//',' &
                        //int2str(l)//','//int2str(is)//') = ' &
                        //real2str(obj%q_mul_out(ia, l, is))
                  end do
               end do
            end do
         case ('delta_q_mul')
            write (unit_rt, '(a)') ' delta_q_mul = '//real2str(obj%delta_q_mul)
         case ('rho_net_in')
            do ia = 1, obj%a%na
               do io1 = 1, obj%e%no_max
                  do io2 = 1, obj%e%no_max
                     do is = 1, obj%a%ns
                        write (unit_rt, '(a)') ' rho_net('//int2str(ia)//',' &
                           //int2str(io1)//','//int2str(io2)//','//int2str(is) &
                           //') = '//cmplx2str(obj%rho_net_in(ia, io1, io2, is))
                     end do
                  end do
               end do
            end do
         case ('rho_net_out')
            do ia = 1, obj%a%na
               do io1 = 1, obj%e%no_max
                  do io2 = 1, obj%e%no_max
                     do is = 1, obj%a%ns
                        write (unit_rt, '(a)') ' rho_net('//int2str(ia)//',' &
                           //int2str(io1)//','//int2str(io2)//','//int2str(is) &
                           //') = '//cmplx2str(obj%rho_net_out(ia, io1, io2, is))
                     end do
                  end do
               end do
            end do
         case ('delta_rho_net')
            write (unit_rt, '(a)') ' delta_rho_net = '//real2str(obj%delta_rho_net)
         case ('rho_net_out_diagonal')
            write (unit_rt, '(a)') ' rho_net_out_diagonal = ' &
               //log2str(obj%rho_net_out_diagonal)
         end select
      end do

      if (tag_rt) then
         write (unit_rt, '(a)') ' /'
      end if

      call TBKOSTER_flush(unit_rt)

      if (.not. present(unit)) then
         close (unit_rt)
      end if
      !deallocate(file_rt,property_rt)
   end subroutine write_txt_formatted

   ! subroutine allocate_orbital_moment(obj)
   !   class(charge),intent(in) :: obj
   !   allocate(LMAT(3,9,9))
   !   LMAT = cmplx(0.0_rp,0.0_rp,kind=rp)
   !   allocate (orbital_moment_mul(obj%a%na,3,3))
   !   orbital_moment_mul = cmplx(0.0_rp,0.0_rp,kind=rp)
   ! end subroutine allocate_orbital_moment

   ! subroutine allocate_orbital_moment(obj)
   !   class(charge),intent(in) :: obj
   !   allocate(LMAT(3,9,9))
   !   LMAT = cmplx(0.0_rp,0.0_rp,kind=rp)
   !   allocate (orbital_moment_mul(obj%a%na,3,3))
   !   orbital_moment_mul = cmplx(0.0_rp,0.0_rp,kind=rp)
   ! end subroutine allocate_orbital_moment

   ! subroutine build_lmat(LMAT)
   !   ! INPUT
   !   ! OUTPUT
   !   complex(rp), dimension(3,9,9), intent(out) :: LMAT
   !   ! LOCAL
   !   integer :: ir,ia,ja
   !
   !   ir=1  !  Lx
   !   do ia=1,9
   !     do ja=1,9
   !       if(ia==3.and.ja==4) then
   !         LMAT(ir,ia,ja)=-i_unit
   !       elseif(ia==4.and.ja==3) then
   !         LMAT(ir,ia,ja)=i_unit
   !       elseif(ia==5.and.ja==7) then
   !         LMAT(ir,ia,ja)=-i_unit
   !       elseif(ia==6.and.ja==8) then
   !         LMAT(ir,ia,ja)=-i_unit
   !       elseif(ia==6.and.ja==9) then
   !         LMAT(ir,ia,ja)=-i_unit*sqrt_three
   !       elseif(ia==7.and.ja==5) then
   !         LMAT(ir,ia,ja)=i_unit
   !       elseif(ia==8.and.ja==6) then
   !         LMAT(ir,ia,ja)=i_unit
   !       elseif(ia==9.and.ja==6) then
   !         LMAT(ir,ia,ja)=i_unit*sqrt_three
   !       end if
   !     end do
   !   end do
   !
   !   ir=2  ! Ly
   !   do ia=1,9
   !     do ja=1,9
   !       if(ia==2.and.ja==4) then
   !         LMAT(ir,ia,ja)=i_unit
   !       elseif(ia==4.and.ja==2) then
   !         LMAT(ir,ia,ja)=-i_unit
   !       elseif(ia==5.and.ja==6) then
   !         LMAT(ir,ia,ja)=i_unit
   !       elseif(ia==6.and.ja==5) then
   !         LMAT(ir,ia,ja)=-i_unit
   !       elseif(ia==7.and.ja==8) then
   !         LMAT(ir,ia,ja)=-i_unit
   !       elseif(ia==7.and.ja==9) then
   !         LMAT(ir,ia,ja)=i_unit*sqrt_three
   !       elseif(ia==8.and.ja==7) then
   !         LMAT(ir,ia,ja)=i_unit
   !       elseif(ia==9.and.ja==7) then
   !         LMAT(ir,ia,ja)=-i_unit*sqrt_three
   !       end if
   !     end do
   !   end do
   !
   !   ir=3  ! Lz
   !   do ia=1,9
   !     do ja=1,9
   !       if(ia==2.and.ja==3) then
   !         LMAT(ir,ia,ja)=-i_unit
   !       elseif(ia==3.and.ja==2) then
   !         LMAT(ir,ia,ja)=i_unit
   !       elseif(ia==5.and.ja==8) then
   !         LMAT(ir,ia,ja)=2*i_unit
   !       elseif(ia==6.and.ja==7) then
   !         LMAT(ir,ia,ja)=i_unit
   !       elseif(ia==7.and.ja==6) then
   !         LMAT(ir,ia,ja)=-i_unit
   !       elseif(ia==8.and.ja==5) then
   !         LMAT(ir,ia,ja)=-2*i_unit
   !       end if
   !     end do
   !   end do
   ! end subroutine build_lmat

   subroutine calculate_orbital_moment(obj, ik, nk, iaos2ih, w, f_k, v_k)
      use math_mod, only: L_x, L_y, L_z
      ! INPUT
      class(charge), intent(inout) :: obj
      integer, intent(in) :: ik, nk
      integer, dimension(:, :, :), intent(in) :: iaos2ih
      real(rp), intent(in) :: w
      real(rp), dimension(:), intent(in) :: f_k
      complex(rp), dimension(:, :, :), intent(in) :: v_k
      ! LOCAL
      integer :: ia, ie, ir, io1, io2, l, kmat, ispin, nn, imat, jmat, nh
      complex(rp), dimension(9, 9) :: LMAT

      nh = size(v_k, 2)

      obj%om = cmplx(0.0_rp, 0.0_rp, kind=rp)
      if (obj%a%ns == 4) then
         do ia = 1, obj%a%na
            ie = obj%a%ia2ie(ia)
            do ir = 1, 3  ! x: ir=1 , y: i=2, z: ir=3
                select case (ir)
                    case (1)
                        LMAT = L_x
                    case (2)
                        LMAT = L_y
                    case (3)
                        LMAT = L_z
                end select
               do io1 = 1, obj%e%no(ie)
                  l = obj%e%o2l(obj%e%o(ie, io1))
                  do kmat = 1, nh
                     nn = ik + (kmat - 1)*nk

                     do ispin = 1, 2
                        imat = iaos2ih(ia, io1, ispin)
                        do io2 = 1, obj%e%no(ie)
                           jmat = iaos2ih(ia, io2, ispin)
                           obj%om(ia, ir, l) = obj%om(ia, ir, l) &
                                               + w*f_k(nn)*obj%a%g_s*conjg(v_k(1, imat, kmat))&
                                               *LMAT(io1, io2)*v_k(2, jmat, kmat)
                        end do
                     end do

                  end do
               end do  !fin de la boucle sur ia
            end do   !fin de la boucle sur ir
         end do   !fin de la boucle sur i
      end if
   end subroutine calculate_orbital_moment

   ! subroutine read_input_q_mul
   !   real(rp)  :: m_sph(3),m_cart(3)
   !   real(rp)  :: n,m,theta,phi
   !   integer   :: ia,l,nAt_file,obj%a%ns_file
   !   character :: cdummy
   !
   !   call open_input_file(mulliken_charge_File_in)
   !
   !   read(mulliken_charge_File_in%Unit,*) nAt_file,obj%a%ns_file
   !   if(nAt_file/=obj%a%na) then
   !     write(*,*) 'Error of nAt in mulliken_charge'
   !     stop
   !   end if
   !   if(obj%a%ns_file/=obj%a%ns) then
   !     write(*,*) 'Error of obj%a%ns in mulliken_charge'
   !     stop
   !   end if
   !
   !   if(obj%a%ns==1) then
   !     do ia=1,obj%a%na
   !       read(mulliken_charge_File_in%Unit,*)
   !       do l=1,3
   !         read(mulliken_charge_File_in%Unit,"(A1,2(1X,D23.16))") cdummy,n,n
   !         q_mul_in(ia,l,0)=n*2
   !       end do
   !     end do
   !   elseif(obj%a%ns==2) then
   !     do ia=1,obj%a%na
   !       read(mulliken_charge_File_in%Unit,*)
   !       do l=1,3
   !         read(mulliken_charge_File_in%Unit,"(A1,2(1X,D23.16))") cdummy,n,m
   !         q_mul_in(ia,l,0)=n
   !         q_mul_in(ia,l,1)=m
   !       end do
   !     end do
   !   elseif(obj%a%ns==4) then
   !     do ia=1,obj%a%na
   !       read(mulliken_charge_File_in%Unit,*)
   !       do l=1,3
   !         read(mulliken_charge_File_in%Unit, *) cdummy,n,m,theta,phi
   !         m_sph(1)=m
   !         m_sph(2)=theta*deg2rad
   !         m_sph(3)=phi*deg2rad
   !         m_cart=sph2cart(m_sph)
   !         q_mul_in(ia,l,0)=n
   !         q_mul_in(ia,l,1:)=m_cart(:)
   !       end do
   !     end do
   !   end if !end if obj%a%ns
   ! end subroutine

   ! subroutine read_input_rho_net
   !   integer   :: ia,ie,io1,io2,ispin,nAt_file,obj%a%ns_file
   !   real(rp)  :: rho_real,rho_imag
   !   integer   :: idummy
   !   character :: cdummy
   !   character(*), parameter :: fmt1 = '(I3,1X,I3)'
   !   character(*), parameter :: fmt2 = '(3X,ES23.16,1X,ES23.16,1X,A)'
   !   !character(*), parameter :: fmt1 = '(I3,1X,I3)'
   !   !character(*), parameter :: fmt2 = '(3X,F19.16,1X,F19.16,1X,A)'
   !
   !   call open_input_file(net_charge_File_in)
   !
   !   read(net_charge_File_in%Unit,*) nAt_file,obj%a%ns_file
   !   if(nAt_file/=obj%a%na) then
   !     write(*,*) 'Error of nAt in rho_net'
   !     stop
   !   end if
   !   if(obj%a%ns_file/=obj%a%ns) then
   !     write(*,*) 'Error of obj%a%ns in rho_net'
   !     stop
   !   end if
   !
   !   do ia=1,obj%a%na
   !     read(net_charge_File_in%Unit,*)
   !     ie = obj%a%ia2ie(ia)
   !     do io1=1,obj%e%no(ie)
   !       do io2=1,obj%e%no(ie)
   !         read(net_charge_File_in%Unit,fmt1,advance='no') idummy, idummy
   !         do ispin = 1,obj%a%ns
   !           read(net_charge_File_in%Unit,fmt2,advance='no') rho_real,rho_imag,&
   !            cdummy
   !           rho_net_in(ia,io1,io2,ispin)=cmplx(rho_real,rho_imag,kind=rp)
   !         end do
   !         read(net_charge_File_in%Unit,*)
   !       end do
   !     end do
   !   end do
   ! end subroutine read_input_rho_net
end module charge_mod
