!
! Copyright (C) 2019
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Pascal Thibaudeau <mailto:pascal.thibaudeau@cea.fr>.
!
! This software is a computer program whose purpose is TBKOSTER.
!
! This software is governed by the CeciLL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/ or redistribute the software under the terms of the CeciLL
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
! knowledge of the CeciLL license and that you accept its terms.
!
!  col2nc.f90
!  TBKOSTER
PROGRAM col2nc
   use, intrinsic :: iso_fortran_env, only: output_unit
   use precision_mod
   use string_mod
   implicit none
   integer, parameter :: unit_charge_col = 10
   integer, parameter :: unit_mag = 11
   integer, parameter :: ns_col = 2, ns_nc = 4
   integer :: iostatus, na, no_max
   character(len=*), parameter :: file_charge_col = 'out_charge.txt'
   character(len=*), parameter :: file_mag = 'mag.txt'
   real(rp), dimension(:, :, :), allocatable :: q_mul, q_mul_nc
   complex(rp), dimension(:, :, :, :), allocatable :: rho_net, rho_net_nc
   real(rp), dimension(:, :), allocatable :: mconfig
   namelist /charge/ q_mul, rho_net
   namelist /mag/ mconfig

   write (*, *) 'enter na and no_max'
   read (*, *) na, no_max

   allocate (q_mul(na, 3, 0:ns_col - 1))
   allocate (q_mul_nc(na, 3, 0:ns_nc - 1))
   allocate (rho_net(na, no_max, no_max, ns_col))
   allocate (rho_net_nc(na, no_max, no_max, ns_nc))
   allocate (mconfig(na, 2))

   open (unit_charge_col, file=file_charge_col, action='read', iostat=iostatus, status='old')
   open (unit_mag, file=file_mag, action='read', iostat=iostatus, status='old')
   read (unit_charge_col, nml=charge)
   read (unit_mag, nml=mag)
   close (unit_charge_col)
   close (unit_mag)

   call charge_col_to_ncol(na, ns_col, ns_nc, no_max, q_mul, q_mul_nc, rho_net, rho_net_nc)
   call rotate_charge_nc(na, ns_nc, no_max, q_mul_nc, rho_net_nc, mconfig)
   call write_charge_nc(na, ns_nc, no_max, q_mul_nc, rho_net_nc)

contains

   !> Read collinear charge and transform into non-collinear in text format from file (default: 'in_charge.txt')
   subroutine charge_col_to_ncol(na, ns_col, ns_nc, no_max, q_mul, q_mul_nc, rho_net, rho_net_nc)
      use precision_mod
      character(len=:), allocatable :: file_rt
      integer :: iostatus
      logical :: isopen
      integer na, ns_col, ns_nc, no_max
      integer :: ia, io1, io2, is, l
      real(rp) :: n, mz
      ! Namelist variables
      real(rp), dimension(na, 3, 0:ns_col - 1) :: q_mul
      real(rp), dimension(na, 3, 0:ns_nc - 1) :: q_mul_nc
      complex(rp), dimension(na, no_max, no_max, ns_col):: rho_net
      complex(rp), dimension(na, no_max, no_max, ns_nc):: rho_net_nc

      q_mul_nc = 0.0_rp
      rho_net_nc = cmplx(0.0_rp, 0.0_rp, kind=rp)

      do ia = 1, na
         do l = 1, 3
            n = 0.0_rp
            mz = 0.0_rp
            do is = 0, 1
               n = n + q_mul(ia, l, is)
               mz = mz + (1 - 2*is)*q_mul(ia, l, is)
            end do
            q_mul_nc(ia, l, 0) = n
            q_mul_nc(ia, l, 3) = mz  ! mx and my are zero.
         end do
      end do

      do ia = 1, na
         do io1 = 1, no_max
            do io2 = 1, no_max
               do is = 1, 2 ! is=3 and 4 are zero for collinear spin.
                  rho_net_nc(ia, io1, io2, is) = rho_net(ia, io1, io2, is)
               end do
            end do
         end do
      end do

   end subroutine charge_col_to_ncol

   subroutine rotate_charge_nc(na, ns_nc, no_max, q_mul_nc, rho_net_nc, mconfig)
      use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
      use precision_mod
      use math_mod, only: sph2cart, rotate_rho, pi
      ! INPUT OUTPUT
      integer :: na, ns_nc, no_max
      real(rp), dimension(na, 3, 0:ns_nc - 1) :: q_mul_nc
      real(rp), dimension(na, 2) :: mconfig
      complex(rp), dimension(na, no_max, no_max, ns_nc) :: rho_net_nc
      ! Local variable
      integer :: ia, io1, io2, is, l
      real(rp) :: mz, n
      real(rp), dimension(3) :: m_sph, m_cart
      complex(rp), dimension(2, 2) ::rho

      do ia = 1, na
         do l = 1, 3
            n = 0.0_rp
            mz = 0.0_rp
            if (q_mul_nc(ia, l, 3) > 0) then
               mz = q_mul_nc(ia, l, 3)
               m_sph(1:3) = (/mz, mconfig(ia, 1), mconfig(ia, 2)/)
            else
               mz = -q_mul_nc(ia, l, 3)
               m_sph(1:3) = (/mz, mconfig(ia, 1) + pi, mconfig(ia, 2) + pi/)
            end if
            m_cart = sph2cart(m_sph)
            q_mul_nc(ia, l, 1:3) = m_cart
         end do
      end do

      do ia = 1, na
         do io1 = 1, no_max
            do io2 = 1, no_max
               rho(1, 1) = rho_net_nc(ia, io1, io2, 1)
               rho(2, 2) = rho_net_nc(ia, io1, io2, 2)
               rho(1, 2) = cmplx(0.0_rp, 0.0_rp, kind=rp)
               rho(2, 1) = cmplx(0.0_rp, 0.0_rp, kind=rp)
               call rotate_rho(rho, mconfig(ia, :))
               rho_net_nc(ia, io1, io2, 1) = rho(1, 1)
               rho_net_nc(ia, io1, io2, 2) = rho(2, 2)
               rho_net_nc(ia, io1, io2, 3) = rho(1, 2)
               rho_net_nc(ia, io1, io2, 4) = rho(2, 1)
            end do
         end do
      end do

   end subroutine rotate_charge_nc

   subroutine write_charge_nc(na, ns_nc, no_max, q_mul_nc, rho_net_nc)
      use precision_mod
      use string_mod
      ! INPUT
      integer :: na, ns_nc, no_max
      real(rp), dimension(na, 3, 0:ns_nc - 1) :: q_mul_nc
      real(rp), dimension(na, 2) :: mconfig
      complex(rp), dimension(na, no_max, no_max, ns_nc) :: rho_net_nc
      ! Local variables
      integer :: ia, io1, io2, is, l
      integer, parameter :: unit_charge_nc = 12
      character(len=*), parameter :: file_charge_nc = 'in_charge.txt'

      open (unit_charge_nc, file=file_charge_nc, action='write')
      write (unit_charge_nc, '(a)') '&charge'
      do ia = 1, na
         do l = 1, 3
            do is = 0, ns_nc - 1
               write (unit_charge_nc, '(a)') ' q_mul('//int2str(ia)//',' &
                  //int2str(l)//','//int2str(is)//') = ' &
                  //real2str(q_mul_nc(ia, l, is))
            end do
         end do
      end do

      do ia = 1, na
         do io1 = 1, no_max
            do io2 = 1, no_max
               do is = 1, ns_nc
                  write (unit_charge_nc, '(a)') ' rho_net('//int2str(ia)//',' &
                     //int2str(io1)//','//int2str(io2)//','//int2str(is) &
                     //') = '//cmplx2str(rho_net_nc(ia, io1, io2, is))
               end do
            end do
         end do
      end do
      write (unit_charge_nc, '(a)') ' /'
      close (unit_charge_nc)
   end subroutine write_charge_nc

end program col2nc
