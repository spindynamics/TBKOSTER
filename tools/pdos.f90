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
!  pdos.f90
!  TBKOSTER
program pdos
   use, intrinsic :: iso_fortran_env, only: output_unit
   use precision_mod
   use string_mod
   implicit none
   integer, parameter :: unit_dos = 10
   integer, parameter :: unit_hamiltonian_tb = 11
   integer, parameter :: unit_energy_scf = 12
   integer :: iostatus, nh, ns, nsl, nen, na_dos
   character(len=*), parameter :: dir = 'dos/'
   character(len=*), parameter :: file_dos = 'out_dos.txt'
   character(len=*), parameter :: file_hamiltonian_tb = 'out_hamiltonian_tb.txt'
   character(len=*), parameter :: file_energy_scf = 'out_energy.txt'
   integer, dimension(:), allocatable :: ia
   real(rp), dimension(:, :), allocatable :: dos_tot
   complex(rp), dimension(:, :, :), allocatable :: dos_s, dos_p, dos_px, dos_py, dos_pz
   complex(rp), dimension(:, :, :), allocatable :: dos_d, dos_dxy, dos_dyz, dos_dzx, dos_dx2y2, dos_dz2r2
   real(rp) :: en_min, en_max, en_f
   namelist /dos/ nen, en_min, en_max, ia, na_dos, dos_tot, dos_s, dos_p, dos_px, dos_py, dos_pz, &
      dos_d, dos_dxy, dos_dyz, dos_dzx, dos_dx2y2, dos_dz2r2
   namelist /hamiltonian_tb/ nh, ns
   namelist /energy/ en_f

   open (unit_dos, file=dir//file_dos, action='read', iostat=iostatus, status='old')
   open (unit_hamiltonian_tb, file=dir//file_hamiltonian_tb, action='read', iostat=iostatus, status='old')

   read (unit_hamiltonian_tb, nml=hamiltonian_tb, iostat=iostatus)
   close (unit_hamiltonian_tb)

   select case (ns)
   case (1)
      nsl = 1
   case (2)
      nsl = 2
   case (4)
      nsl = 1
   end select

   allocate (ia(0))
   allocate (dos_tot(0, 0), dos_s(0, 0, 0), dos_p(0, 0, 0), dos_px(0, 0, 0), dos_py(0, 0, 0), dos_pz(0, 0, 0))
   allocate (dos_d(0, 0, 0), dos_dxy(0, 0, 0), dos_dyz(0, 0, 0), dos_dzx(0, 0, 0), dos_dx2y2(0, 0, 0), dos_dz2r2(0, 0, 0))
   read (unit_dos, nml=dos, iostat=iostatus)
   deallocate (ia)
   deallocate (dos_tot, dos_s, dos_p, dos_px, dos_py, dos_pz)
   deallocate (dos_d, dos_dxy, dos_dyz, dos_dzx, dos_dx2y2, dos_dz2r2)

   allocate (ia(na_dos))
   allocate (dos_tot(nen, nsl))
   allocate (dos_s(na_dos, nen, ns), dos_p(na_dos, nen, ns))
   allocate (dos_px(na_dos, nen, ns), dos_py(na_dos, nen, ns), dos_pz(na_dos, nen, ns))
   allocate (dos_d(na_dos, nen, ns), dos_dxy(na_dos, nen, ns), dos_dyz(na_dos, nen, ns), dos_dzx(na_dos, nen, ns), &
             dos_dx2y2(na_dos, nen, ns), dos_dz2r2(na_dos, nen, ns))
   rewind (unit_dos)
   read (unit_dos, nml=dos, iostat=iostatus)
   open (unit_energy_scf, file=file_energy_scf, action='read', iostat=iostatus, status='old')
   read (unit_energy_scf, nml=energy, iostat=iostatus)
   close (unit_dos)
   close (unit_energy_scf)

   call build_pdos(en_f, nen, na_dos, ns, nsl, ia, en_min, en_max, dos_tot, dos_s, dos_p, dos_px, dos_py, dos_pz, &
                   dos_d, dos_dxy, dos_dyz, dos_dzx, dos_dx2y2, dos_dz2r2)
contains
   subroutine build_pdos(en_f, nen, na_dos, ns, nsl, ia, en_min, en_max, &
                         dos_tot, dos_s, dos_p, dos_px, dos_py, dos_pz, &
                         dos_d, dos_dxy, dos_dyz, dos_dzx, dos_dx2y2, dos_dz2r2)
      use precision_mod
      use string_mod
      use math_mod, only: i_unit
      implicit none
      integer :: nen, na_dos, ns, nsl
      integer :: ien, ia_dos, isl, ispin_up, ispin_dn, g_s
      integer, parameter :: unit_dos_tot = 10
      integer, parameter :: unit_pdos_s = 11
      integer, parameter :: unit_pdos_p = 12
      integer, parameter :: unit_pdos_d = 13
      integer, parameter :: unit_pdos_spd = 14
      integer :: ia(na_dos)
      real(rp) :: en, den
      real(rp) :: en_min, en_max, en_f
      real(rp) :: dos_tot(nen, nsl)
      complex(rp) :: dos_s(na_dos, nen, ns), dos_p(na_dos, nen, ns), &
                     dos_px(na_dos, nen, ns), dos_py(na_dos, nen, ns), dos_pz(na_dos, nen, ns)
      complex(rp) :: dos_d(na_dos, nen, ns), dos_dxy(na_dos, nen, ns), dos_dyz(na_dos, nen, ns), dos_dzx(na_dos, nen, ns), &
                     dos_dx2y2(na_dos, nen, ns), dos_dz2r2(na_dos, nen, ns)
      character(len=*), parameter :: dir = 'dos/'
      character(len=*), parameter :: file_dos_tot = 'dos-tot.dat'
      character(len=*), parameter :: file_pdos_s = 'pdos-s.dat'
      character(len=*), parameter :: file_pdos_p = 'pdos-p.dat'
      character(len=*), parameter :: file_pdos_d = 'pdos-d.dat'
      character(len=*), parameter :: file_pdos_spd = 'pdos-spd.dat'
      character(len=80) :: fmt
      real(rp) ::   dos_tot_x, dos_tot_y, dos_tot_z, dos_tot_n
      real(rp) ::   dos_s_x, dos_s_y, dos_s_z, dos_s_n
      real(rp) ::   dos_p_x, dos_p_y, dos_p_z, dos_p_n
      real(rp) ::   dos_px_x, dos_px_y, dos_px_z, dos_px_n
      real(rp) ::   dos_py_x, dos_py_y, dos_py_z, dos_py_n
      real(rp) ::   dos_pz_x, dos_pz_y, dos_pz_z, dos_pz_n
      real(rp) ::   dos_d_x, dos_d_y, dos_d_z, dos_d_n
      real(rp) ::   dos_dxy_x, dos_dxy_y, dos_dxy_z, dos_dxy_n
      real(rp) ::   dos_dyz_x, dos_dyz_y, dos_dyz_z, dos_dyz_n
      real(rp) ::   dos_dzx_x, dos_dzx_y, dos_dzx_z, dos_dzx_n
      real(rp) ::   dos_dx2y2_x, dos_dx2y2_y, dos_dx2y2_z, dos_dx2y2_n
      real(rp) ::   dos_dz2r2_x, dos_dz2r2_y, dos_dz2r2_z, dos_dz2r2_n
      real(rp) ::   dos_spd_x, dos_spd_y, dos_spd_z, dos_spd_n

      den = (en_max - en_min)/nen
      if (ns == 1 .or. ns == 2) then
         if (ns == 1) then
            ispin_up = 1
            ispin_dn = 1
            g_s = 2
         elseif (ns == 2) then
            ispin_up = 1
            ispin_dn = 2
            g_s = 1
         end if
      end if

      open (unit=unit_dos_tot, file=dir//file_dos_tot, action='write')
      open (unit=unit_pdos_s, file=dir//file_pdos_s, action='write')
      open (unit=unit_pdos_p, file=dir//file_pdos_p, action='write')
      open (unit=unit_pdos_d, file=dir//file_pdos_d, action='write')
      open (unit=unit_pdos_spd, file=dir//file_pdos_spd, action='write')

      if (ns == 1 .or. ns == 4) then
         write (unit_dos_tot, *) '@# e   total dos'
         do ien = 1, nen
            en = en_min + (ien - 1)*den
            write (unit_dos_tot, '(2(a,1x))') real2str(en), real2str(dos_tot(ien, 1))
         end do
      elseif (ns == 2) then
         write (unit_dos_tot, *) '@# e   total dos up and down'
         do ien = 1, nen
            en = en_min + (ien - 1)*den
            write (unit_dos_tot, '(3(a,1x))') real2str(en), real2str(dos_tot(ien, 1)), real2str(dos_tot(ien, 2))
         end do
      end if
      write (unit_dos_tot, *) '@# ef'
      write (unit_dos_tot, '(2f8.3)') en_f, 0.0_rp
      write (unit_dos_tot, '(2f8.3)') en_f, 10.0_rp

      if (na_dos > 0) then
         if (ns == 1 .or. ns == 2) then
            do ia_dos = 1, na_dos
               write (unit_pdos_s, *) '@#dos of site', ia(ia_dos)
               write (unit_pdos_s, "(a)", advance="no") "@# s(↑) s(↓)"
               write (unit_pdos_p, *) '@#dos of site', ia(ia_dos)
               write (unit_pdos_p, "(a)", advance="no") "@# p(↑) p(↓) px(↑) px(↓) py(↑) py(↓) pz(↑) pz(↓)"
               write (unit_pdos_d, *) '@#dos of site', ia(ia_dos)
               write (unit_pdos_d, "(a)", advance="no") "@# d(↑) d(↓) dxy(↑) dxy(↓) dyz(↑) dyz(↓) dxz(↑) "&
                  "dxz(↓) dy2y2(↑) dx2y2(↓) dz2(↑) dz2(↓)"
               write (unit_pdos_spd, *) '@#dos of site', ia(ia_dos)
               write (unit_pdos_spd, "(a)", advance="no") "@# s(↑) s(↓) p(↑) p(↓) d(↑) d(↓)"
               do ien = 1, nen
                  en = en_min + (ien - 1)*den
                  write (unit_pdos_s, '(3(a,1x))') real2str(en), real2str(real(dos_s(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_s(ia_dos, ien, ispin_dn)/g_s))
                  write (unit_pdos_p, '(9(a,1x))') real2str(en), real2str(real(dos_p(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_p(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_px(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_px(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_py(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_py(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_pz(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_pz(ia_dos, ien, ispin_dn)/g_s))
                  write (unit_pdos_d, '(13(a,1x))') real2str(en), real2str(real(dos_d(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_d(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_dxy(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_dxy(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_dyz(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_dyz(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_dzx(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_dzx(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_dx2y2(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_dx2y2(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_dz2r2(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_dz2r2(ia_dos, ien, ispin_dn)/g_s))

                  write (unit_pdos_spd, '(7(a,1x))') real2str(en), real2str(real(dos_s(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_s(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_p(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_p(ia_dos, ien, ispin_dn)/g_s)), &
                     real2str(real(dos_d(ia_dos, ien, ispin_up)/g_s)), &
                     real2str(real(dos_d(ia_dos, ien, ispin_dn)/g_s))
               end do
            end do
         elseif (ns == 4) then
            do ia_dos = 1, na_dos
               write (unit_pdos_s, *) '@#dos of site', ia(ia_dos)
               write (unit_pdos_s, "(a)", advance="no") "@# s(n) s(x) s(y) s(z)"
               write (unit_pdos_p, *) '@#dos of site', ia(ia_dos)
               write (unit_pdos_p, "(a)", advance="no") "@# p(n) p(x) p(y) p(z) px(n) px(x) px(y) px(z) "&
                  "py(n) py(x) py(y) py(z) pz(n) pz(x) pz(y) pz(z)"
               write (unit_pdos_d, *) '@#dos of site', ia(ia_dos)
               write (unit_pdos_d, "(a)", advance="no") "@# d(n) d(x) d(y) d(z) dxy(n) dxy(x) dxy(y) dxy(z) "&
                  "dyz(n) dyz(x) dyz(y) dyz(z) dxz(n) dxz(x) dxz(y) dxz(z) dx2y2(n) dx2y2(x) dx2y2(y) dx2y2(z) "&
                  "dz2(n) dz2(x) dz2(y) dz2(z)"
               write (unit_pdos_spd, *) '@#dos of site', ia(ia_dos)
               write (unit_pdos_spd, "(a)", advance="no") "@# s(n) s(x) s(y) s(z) p(n) p(x) p(y) p(z) "&
                  "d(n) d(x) d(y) d(z)"

               do ien = 1, nen
                  en = en_min + (ien - 1)*den

                  dos_s_n = dos_s(ia_dos, ien, 1) + dos_s(ia_dos, ien, 2)
                  dos_s_x = dos_s(ia_dos, ien, 3) + dos_s(ia_dos, ien, 4)
                  dos_s_y = i_unit*(dos_s(ia_dos, ien, 3) - dos_s(ia_dos, ien, 4))
                  dos_s_z = dos_s(ia_dos, ien, 1) - dos_s(ia_dos, ien, 2)

                  dos_p_n = dos_p(ia_dos, ien, 1) + dos_p(ia_dos, ien, 2)
                  dos_p_x = dos_p(ia_dos, ien, 3) + dos_p(ia_dos, ien, 4)
                  dos_p_y = i_unit*(dos_p(ia_dos, ien, 3) - dos_p(ia_dos, ien, 4))
                  dos_p_z = dos_p(ia_dos, ien, 1) - dos_p(ia_dos, ien, 2)

                  dos_px_n = dos_px(ia_dos, ien, 1) + dos_px(ia_dos, ien, 2)
                  dos_px_x = dos_px(ia_dos, ien, 3) + dos_px(ia_dos, ien, 4)
                  dos_px_y = i_unit*(dos_px(ia_dos, ien, 3) - dos_px(ia_dos, ien, 4))
                  dos_px_z = dos_px(ia_dos, ien, 1) - dos_px(ia_dos, ien, 2)

                  dos_py_n = dos_py(ia_dos, ien, 1) + dos_py(ia_dos, ien, 2)
                  dos_py_x = dos_py(ia_dos, ien, 3) + dos_py(ia_dos, ien, 4)
                  dos_py_y = i_unit*(dos_py(ia_dos, ien, 3) - dos_py(ia_dos, ien, 4))
                  dos_py_z = dos_py(ia_dos, ien, 1) - dos_py(ia_dos, ien, 2)

                  dos_pz_n = dos_pz(ia_dos, ien, 1) + dos_pz(ia_dos, ien, 2)
                  dos_pz_x = dos_pz(ia_dos, ien, 3) + dos_pz(ia_dos, ien, 4)
                  dos_pz_y = i_unit*(dos_pz(ia_dos, ien, 3) - dos_pz(ia_dos, ien, 4))
                  dos_pz_z = dos_pz(ia_dos, ien, 1) - dos_pz(ia_dos, ien, 2)

                  dos_d_n = dos_d(ia_dos, ien, 1) + dos_d(ia_dos, ien, 2)
                  dos_d_x = dos_d(ia_dos, ien, 3) + dos_d(ia_dos, ien, 4)
                  dos_d_y = i_unit*(dos_d(ia_dos, ien, 3) - dos_d(ia_dos, ien, 4))
                  dos_d_z = dos_d(ia_dos, ien, 1) - dos_d(ia_dos, ien, 2)

                  dos_dxy_n = dos_dxy(ia_dos, ien, 1) + dos_dxy(ia_dos, ien, 2)
                  dos_dxy_x = dos_dxy(ia_dos, ien, 3) + dos_dxy(ia_dos, ien, 4)
                  dos_dxy_y = i_unit*(dos_dxy(ia_dos, ien, 3) - dos_dxy(ia_dos, ien, 4))
                  dos_dxy_z = dos_dxy(ia_dos, ien, 1) - dos_dxy(ia_dos, ien, 2)

                  dos_dyz_n = dos_dyz(ia_dos, ien, 1) + dos_dyz(ia_dos, ien, 2)
                  dos_dyz_x = dos_dyz(ia_dos, ien, 3) + dos_dyz(ia_dos, ien, 4)
                  dos_dyz_y = i_unit*(dos_dyz(ia_dos, ien, 3) - dos_dyz(ia_dos, ien, 4))
                  dos_dyz_z = dos_dyz(ia_dos, ien, 1) - dos_dyz(ia_dos, ien, 2)

                  dos_dzx_n = dos_dzx(ia_dos, ien, 1) + dos_dzx(ia_dos, ien, 2)
                  dos_dzx_x = dos_dzx(ia_dos, ien, 3) + dos_dzx(ia_dos, ien, 4)
                  dos_dzx_y = i_unit*(dos_dzx(ia_dos, ien, 3) - dos_dzx(ia_dos, ien, 4))
                  dos_dzx_z = dos_dzx(ia_dos, ien, 1) - dos_dzx(ia_dos, ien, 2)

                  dos_dx2y2_n = dos_dx2y2(ia_dos, ien, 1) + dos_dx2y2(ia_dos, ien, 2)
                  dos_dx2y2_x = dos_dx2y2(ia_dos, ien, 3) + dos_dx2y2(ia_dos, ien, 4)
                  dos_dx2y2_y = i_unit*(dos_dx2y2(ia_dos, ien, 3) - dos_dx2y2(ia_dos, ien, 4))
                  dos_dx2y2_z = dos_dx2y2(ia_dos, ien, 1) - dos_dx2y2(ia_dos, ien, 2)

                  dos_dz2r2_n = dos_dz2r2(ia_dos, ien, 1) + dos_dz2r2(ia_dos, ien, 2)
                  dos_dz2r2_x = dos_dz2r2(ia_dos, ien, 3) + dos_dz2r2(ia_dos, ien, 4)
                  dos_dz2r2_y = i_unit*(dos_dz2r2(ia_dos, ien, 3) - dos_dz2r2(ia_dos, ien, 4))
                  dos_dz2r2_z = dos_dz2r2(ia_dos, ien, 1) - dos_dz2r2(ia_dos, ien, 2)

                  dos_spd_n = dos_s_n + dos_p_n + dos_d_n
                  dos_spd_x = dos_s_x + dos_p_x + dos_d_x
                  dos_spd_y = dos_s_y + dos_p_y + dos_d_y
                  dos_spd_z = dos_s_z + dos_p_z + dos_d_z

                  write (unit_pdos_s, '(5(a,1x))') real2str(en), real2str(dos_s_n), &
                     real2str(dos_s_x), real2str(dos_s_y), real2str(dos_s_z)
                  write (unit_pdos_p, '(17(a,1x))') real2str(en), real2str(dos_p_n), &
                     real2str(dos_p_x), real2str(dos_p_y), real2str(dos_p_z), &
                     real2str(dos_px_n), real2str(dos_px_x), &
                     real2str(dos_px_y), real2str(dos_px_z), &
                     real2str(dos_py_n), real2str(dos_py_x), &
                     real2str(dos_py_y), real2str(dos_py_z), &
                     real2str(dos_pz_n), real2str(dos_pz_x), &
                     real2str(dos_pz_y), real2str(dos_pz_z)
                  write (unit_pdos_d, '(25(a,1x))') real2str(en), real2str(dos_d_n), real2str(dos_d_x), &
                     real2str(dos_d_y), real2str(dos_d_z), &
                     real2str(dos_dxy_n), real2str(dos_dxy_x), &
                     real2str(dos_dxy_y), real2str(dos_dxy_z), &
                     real2str(dos_dyz_n), real2str(dos_dyz_x), &
                     real2str(dos_dyz_y), real2str(dos_dyz_z), &
                     real2str(dos_dzx_n), real2str(dos_dzx_x), &
                     real2str(dos_dzx_y), real2str(dos_dzx_z), &
                     real2str(dos_dx2y2_n), real2str(dos_dx2y2_x), &
                     real2str(dos_dx2y2_y), real2str(dos_dx2y2_z), &
                     real2str(dos_dz2r2_n), real2str(dos_dz2r2_x), &
                     real2str(dos_dz2r2_y), real2str(dos_dz2r2_z)
                  write (unit_pdos_spd, '(13(a,1x))') real2str(en), real2str(dos_s_n), real2str(dos_s_x), &
                     real2str(dos_s_y), real2str(dos_s_z), &
                     real2str(dos_p_n), real2str(dos_p_x), &
                     real2str(dos_p_y), real2str(dos_p_z), &
                     real2str(dos_d_n), real2str(dos_d_x), &
                     real2str(dos_d_y), real2str(dos_d_z)
               end do
            end do
         end if
         write (unit_pdos_s, *) '@# ef'
         write (unit_pdos_s, '(2f8.3)') en_f, 0.0_rp
         write (unit_pdos_s, '(2f8.3)') en_f, 10.0_rp
         write (unit_pdos_p, *) '@# ef'
         write (unit_pdos_p, '(2f8.3)') en_f, 0.0_rp
         write (unit_pdos_p, '(2f8.3)') en_f, 10.0_rp
         write (unit_pdos_d, *) '@# ef'
         write (unit_pdos_d, '(2f8.3)') en_f, 0.0_rp
         write (unit_pdos_d, '(2f8.3)') en_f, 10.0_rp
         write (unit_pdos_spd, *) '@# ef'
         write (unit_pdos_spd, '(2f8.3)') en_f, 0.0_rp
         write (unit_pdos_spd, '(2f8.3)') en_f, 10.0_rp

      end if
   end subroutine build_pdos

end program pdos
