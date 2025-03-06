!
! Copyright (C) 2019
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
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
!  bands.f90
!  TBKOSTER
program bands
   use, intrinsic :: iso_fortran_env, only: output_unit
   use precision_mod
   use string_mod
   implicit none
   integer, parameter :: unit_mesh = 11
   integer, parameter :: unit_mesh_out = 111
   integer, parameter :: unit_hamiltonian_tb = 12
   integer, parameter :: unit_band_in = 13
   integer, parameter :: unit_band_out = 14
   integer, parameter :: unit_log_out = 15
   integer :: iostatus, nx, nh, ns, nsl, no_max, na_band
   logical :: isopen, close_file
   character(len=*), parameter :: dir = 'band/'
   character(len=*), parameter :: file_band_in = 'in_band.txt'
   character(len=*), parameter :: file_band_out = 'out_band.txt'
   character(len=*), parameter :: file_mesh_in = 'in_mesh.txt'
   character(len=*), parameter :: file_mesh_out = 'out_mesh.txt'
   character(len=*), parameter :: file_hamiltonian_tb = 'out_hamiltonian_tb.txt'
   character(len=*), parameter :: file_log_out = 'out_log.txt'
   character(len=9) :: x_coord
   character(len=4) :: type, type_case
   character(len=80) :: line
   character(len=11) :: proj
   integer, dimension(:), allocatable :: ia_band, iband2io
   real(rp), dimension(:, :, :), allocatable :: en_k
   real(rp), dimension(:), allocatable :: w
   real(rp), dimension(:, :), allocatable :: x
   real(rp) :: v_factor
   real(rp), dimension(3, 3) :: v, vrec
   real(rp), dimension(:, :, :, :, :), allocatable :: w_band_site
   real(rp), dimension(:, :, :, :), allocatable :: w_band_spin
   real(rp), dimension(:, :, :, :), allocatable :: w_band_orb
   real(rp) :: en_min, en_max, en_f, Eref
   namelist /lattice/ v_factor, v, vrec
   namelist /mesh/ type
   namelist /mesh_out/ type, nx, x_coord, x, w
   namelist /hamiltonian_tb/ nh, ns
   namelist /band/ proj, na_band, ia_band
   namelist /band_out/ en_k, iband2io, w_band_site, w_band_spin, w_band_orb

   ! inquire(unit=unit_energy_in,opened=isopen)
   ! if (isopen) then
   !   write(*,'(a)') 'energy%read_txt() : unit 14 is already open'
   !   close(unit_energy_in)
   ! else
   !   open(unit=unit_energy_in,file=dir//file_energy_in,action='read',iostat=iostatus,status='old')
   ! end if
   ! if(iostatus /= 0) then
   !   write(*,*) 'energy%read_txt(): file ', file_energy_in, ' not found'
   !   error stop
   !  end if
   open (unit_band_in, file=dir//file_band_in, action='read', iostat=iostatus, status='old')
   open (unit_band_out, file=dir//file_band_out, action='read', iostat=iostatus, status='old')
   open (unit_mesh, file=dir//file_mesh_in, action='read', iostat=iostatus, status='old')
   open (unit_mesh_out, file=dir//file_mesh_out, action='read', iostat=iostatus, status='old')
   open (unit_hamiltonian_tb, file=dir//file_hamiltonian_tb, action='read', iostat=iostatus, status='old')
   open (unit_log_out, file=dir//file_log_out, action='read', iostat=iostatus, status='old')

   read (unit_log_out, nml=lattice, iostat=iostatus)
   close (unit_log_out)
   ! read k mesh in (test type of mesh)
   read (unit_mesh, nml=mesh, iostat=iostatus)
   type_case = lower(type)
   close (unit_mesh)
   ! read k mesh out
   allocate (x(0, 0), w(0))
   read (unit_mesh_out, nml=mesh_out, iostat=iostatus)
   deallocate (x, w)
   select case (type_case)
   case ('path')
      write (*, *) 'path is suited for band structure plot'
   case ('list')
      write (*, *) 'list is suited for 3D/2D maps+vector field'
   case ('mp')
      write (*, *) 'monkhorst pack is not suited for bands analysis'
      stop
   end select

   allocate (x(nx, 3))
   allocate (w(nx))
   x = 0.0_rp
   w = 1.0_rp/nx
   rewind (unit_mesh_out)
   read (unit_mesh_out, nml=mesh_out)
   close (unit_mesh_out)
   x_coord = lower(x_coord)
   if (x_coord == 'direct') then
      write (*, *) 'x_coord must be cartesian'
      stop
   end if

   close (unit_mesh_out)

   ! read information about the hamiltonian size
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
   na_band = 0
   allocate (ia_band(0))
   read (unit_band_in, nml=band, iostat=iostatus)
   deallocate (ia_band)
   allocate (ia_band(na_band))
   rewind (unit_band_in)
   read (unit_band_in, nml=band, iostat=iostatus)
   close (unit_band_in)

   ! read eigenvalues and different types of projections: sites or spin or orbital
   allocate (en_k(nh, nx, nsl))
   read (unit_band_out, nml=band_out, iostat=iostatus)

   rewind (unit_band_out)
   if (TRIM(proj) == 'site' .AND. na_band > 0) then
      allocate (iband2io(na_band))
      iband2io = 0
      allocate (w_band_site(0, 0, 0, 0, 0))
      read (unit_band_out, nml=band_out, iostat=iostatus)
      no_max = maxval(iband2io)
      deallocate (w_band_site)
      rewind (unit_band_out)
      allocate (w_band_site(nh, nx, na_band, no_max, nsl))
      read (unit_band_out, nml=band_out, iostat=iostatus)
   elseif (TRIM(proj) == 'spin' .AND. na_band > 0) then
      allocate (w_band_spin(0, 0, 0, 0))
      read (unit_band_out, nml=band_out, iostat=iostatus)
      deallocate (w_band_spin)
      rewind (unit_band_out)
      allocate (w_band_spin(nh, nx, 0:na_band, 4))
      read (unit_band_out, nml=band_out, iostat=iostatus)
   elseif (TRIM(proj) == 'orbit' .AND. na_band > 0) then
      allocate (w_band_orb(0, 0, 0, 0))
      read (unit_band_out, nml=band_out, iostat=iostatus)
      deallocate (w_band_orb)
      rewind (unit_band_out)
      allocate (w_band_orb(nh, nx, 0:na_band, 4))
      read (unit_band_out, nml=band_out, iostat=iostatus)
   elseif (TRIM(proj) == 'spin,orbit' .AND. na_band > 0) then
      allocate (w_band_spin(0, 0, 0, 0))
      allocate (w_band_orb(0, 0, 0, 0))
      read (unit_band_out, nml=band_out, iostat=iostatus)
      deallocate (w_band_spin)
      deallocate (w_band_orb)
      rewind (unit_band_out)
      allocate (w_band_spin(nh, nx, 0:na_band, 4))
      allocate (w_band_orb(nh, nx, 0:na_band, 4))
      read (unit_band_out, nml=band_out, iostat=iostatus)
   end if

   close (unit_band_out)
   ! set the Fermi energy to zero
   call get_fermi_scf(en_f)
   en_k = en_k - en_f

   select case (type_case)
   case ('path')
      close_file = .false.
      call build_band_path(x, en_k, nh, nx, nsl)
      if (TRIM(proj) == 'site' .AND. na_band > 0) then
         call build_band_path_site_orb(x, en_k, w_band_site, iband2io, na_band, no_max, nh, nx, nsl)
      elseif (TRIM(proj) == 'spin') then
         close_file = .true.
         call build_band_path_spinorb(close_file, x, en_k, w_band_spin, na_band, nh, nx)
      elseif (TRIM(proj) == 'orbit') then
         close_file = .true.
         call build_band_path_spinorb(close_file, x, en_k, w_band_orb, na_band, nh, nx)
      elseif (TRIM(proj) == 'spin,orbit') then
         call build_band_path_spinorb(close_file, x, en_k, w_band_spin, na_band, nh, nx)
         close_file = .true.
         call build_band_path_spinorb(close_file, x, en_k, w_band_orb, na_band, nh, nx)
      end if
   case ('list')
      write (*, *) 'list'
      close_file = .false.
      Eref = 0.0_rp
      if (TRIM(proj) == 'spin') then
         close_file = .true.
         call build_vector_field(close_file, x, en_k, Eref, w_band_spin, na_band, nh, nx)
      elseif (TRIM(proj) == 'orbit') then
         close_file = .true.
         call build_vector_field(close_file, x, en_k, Eref, w_band_orb, na_band, nh, nx)
      elseif (TRIM(proj) == 'spin,orbit') then
         call build_vector_field(close_file, x, en_k, Eref, w_band_spin, na_band, nh, nx)
         close_file = .true.
         call build_vector_field(close_file, x, en_k, Eref, w_band_orb, na_band, nh, nx)
      end if
   end select

contains

   subroutine build_band_path(x, en_k, nh, nx, nsl)
      use precision_mod
      use string_mod
      implicit none
      integer, intent(in) :: nh, nx, nsl
      real(rp), intent(in) :: x(nx, 3), en_k(nh, nx, nsl)

      integer :: ih, ix, isl
      integer, parameter :: unit_band = 10
      real(rp) :: sk
      character(len=*), parameter :: dir = 'band/'
      character(len=*), parameter :: file_band = 'band.dat'
      character(len=80) :: fmt

      open (unit=unit_band, file=dir//file_band, action='write')
      fmt = trim('('//int2str(nh + 1)//'f14.7'//')')

      do isl = 1, nsl
         write (unit_band, *) '@# k   band (ev)'
         sk = 0.0_rp
         write (unit_band, fmt) sk, (en_k(ih, 1, isl), ih=1, nh)
         do ix = 2, nx
            sk = sk + sqrt(sum((x(ix, :) - x(ix - 1, :))**2))
            write (unit_band, fmt) sk, (en_k(ih, ix, isl), ih=1, nh)
         end do
      end do
      write (unit_band, *) '@# k   ef=0'
      write (unit_band, '(2f8.3)') 0.0_rp, 0.0_rp
      write (unit_band, '(2f8.3)') sk, 0.0_rp
      close (unit_band)
   end subroutine build_band_path

   subroutine build_band_path_site_orb(x, en_k, w_band_site, iband2io, na_band, no_max, nh, nx, nsl)
      use precision_mod
      use string_mod
      implicit none
      integer, intent(in) :: na_band, no_max, nh, nx, nsl
      integer, intent(in) :: iband2io(na_band)
      real(rp), intent(in) :: x(nx, 3), en_k(nh, nx, nsl)
      real(rp), intent(in) :: w_band_site(nh, nx, na_band, no_max, nsl)
      integer :: ih, ix, isl, ia_band, io
      integer, parameter :: unit_band = 10
      real(rp) :: sk
      character(len=*), parameter :: dir = 'band/'
      character(len=*), parameter :: file_band_weight = 'band_weight_site_orb.dat'
      character(len=80) :: fmt

      open (unit=unit_band, file=dir//file_band_weight, action='write')

      do isl = 1, nsl
         write (unit_band, *) '@# k   band (ev)'
         if (isl > 1) then
            write (unit_band, *)
            write (unit_band, *)
         end if
         do ia_band = 1, na_band
            if (ia_band > 1) then
               write (unit_band, *)
               write (unit_band, *)
            end if
            write (unit_band, *) '@# atom no', ia_band, 'number of orbitals', iband2io(ia_band)
            sk = 0.0_rp
            fmt = trim('('//int2str(3 + iband2io(ia_band))//'f12.7'//')')
            do ih = 1, nh
               write (unit_band, fmt) sk, en_k(ih, 1, isl), sum(w_band_site(ih, 1, ia_band, :, isl)), &
                  (w_band_site(ih, 1, ia_band, io, isl), io=1, iband2io(ia_band))
            end do
            do ix = 2, nx
               sk = sk + sqrt(sum((x(ix, :) - x(ix - 1, :))**2))
               do ih = 1, nh
                  write (unit_band, fmt) sk, en_k(ih, ix, isl), sum(w_band_site(ih, ix, ia_band, :, isl)), &
                     (w_band_site(ih, ix, ia_band, io, isl), io=1, iband2io(ia_band))
               end do
            end do
         end do
      end do
      !    write (unit_band, *) '@# k   ef=0'
      !    write (unit_band, '(2f8.3)') 0.0_rp, 0.0_rp
      !    write (unit_band, '(2f8.3)') sk, 0.0_rp
      close (unit_band)
   end subroutine build_band_path_site_orb

   subroutine build_band_path_spinorb(close_file, x, en_k, w_band, na_band, nh, nx)
      use precision_mod
      use string_mod
      implicit none
      integer, intent(in) :: na_band, nh, nx
      real(rp), intent(in) :: x(nx, 3), en_k(nh, nx, 1)
      real(rp), intent(in) :: w_band(nh, nx, 0:na_band, 4)
      logical :: isopen
      integer :: ih, ix, ii, ia_band
      logical :: close_file
      integer, parameter :: unit_band = 10
      real(rp) :: sk
      character(len=*), parameter :: dir = 'band/'
      character(len=*), parameter :: file_band_weight = 'band_weight_spinorb.dat'
      character(len=80) :: fmt

      inquire (unit=unit_band, opened=isopen)
      if (.not. isopen) then
         open (unit=unit_band, file=dir//file_band_weight, action='write')
      end if
      write (unit_band, *) '@# k   band (ev)'
      do ia_band = 0, na_band
         if (ia_band > 0) then
            write (unit_band, *)
            write (unit_band, *)
         end if
         if (ia_band == 0) then
            write (unit_band, *) '@# total spin/orb', ia_band
         else
            write (unit_band, *) '@# atom no', ia_band
         end if
         sk = 0.0_rp
         fmt = trim('(6f12.7)')
         do ih = 1, nh
            write (unit_band, fmt) sk, en_k(ih, 1, 1), (w_band(ih, 1, ia_band, ii), ii=1, 4)
         end do
         do ix = 2, nx
            sk = sk + sqrt(sum((x(ix, :) - x(ix - 1, :))**2))
            do ih = 1, nh
               write (unit_band, fmt) sk, en_k(ih, ix, 1), (w_band(ih, ix, ia_band, ii), ii=1, 4)
            end do
         end do
      end do
      !    write (unit_band, *) '@# k   ef=0'
      !    write (unit_band, '(2f8.3)') 0.0_rp, 0.0_rp
      !    write (unit_band, '(2f8.3)') sk, 0.0_rp
      if (close_file) close (unit_band)
   end subroutine build_band_path_spinorb

   subroutine build_vector_field(close_file, x, en_k, Eref, w_band, na_band, nh, nx)
      use precision_mod
      use string_mod
      implicit none
      integer, intent(in) :: na_band, nh, nx
      real(rp), intent(in) :: x(nx, 3), en_k(nh, nx, 1)
      real(rp), intent(in) :: w_band(nh, nx, 0:na_band, 4)
      real(rp), intent(in) :: Eref
      logical :: isopen
      integer :: ih, ix, ii, ia_band
      logical :: close_file
      integer, parameter :: unit_vf = 10, unit_fermi = 11
      real(rp) :: average0, average(4)
      character(len=*), parameter :: dir = 'band/'
      character(len=*), parameter :: file_fermi = 'fermi.dat'
      character(len=*), parameter :: file_vector_field = 'vector_field.dat'
      character(len=80) :: fmt

      inquire (unit=unit_fermi, opened=isopen)
      if (.not. isopen) then
         open (unit=unit_fermi, file=dir//file_fermi, action='write')
      end if

      fmt = trim('(4f12.7)')
      do ix = 1, nx
         average0 = 0.0_rp
         do ih = 1, nh
            average0 = average0 + delta_function(en_k(ih, ix, 1) - Eref)
         end do
         write (unit_fermi, fmt) x(ix, 1), x(ix, 2), x(ix, 3), average0
      end do

      if (close_file) close (unit_fermi)

      inquire (unit=unit_vf, opened=isopen)
      if (.not. isopen) then
         open (unit=unit_vf, file=dir//file_vector_field, action='write')
      end if

      fmt = trim('(7f12.7)')
      do ix = 1, nx
         average(1:4) = 0.0_rp
         do ih = 1, nh
            average(2:4) = average(2:4) + w_band(ih, ix, 0, 2:4)*delta_function(en_k(ih, ix, 1) - Eref)
         end do
         average(1) = norm2(average(2:4))
         write (unit_vf, fmt) x(ix, 1), x(ix, 2), x(ix, 3), (average(ii), ii=1, 4)
      end do
      !    write (unit_band, *) '@# k   ef=0'
      !    write (unit_band, '(2f8.3)') 0.0_rp, 0.0_rp
      !    write (unit_band, '(2f8.3)') sk, 0.0_rp
      if (close_file) close (unit_vf)
   end subroutine build_vector_field

   subroutine get_fermi_scf(en_f)
      use precision_mod
      implicit none
      integer, parameter :: unit_energy_scf = 10
      integer :: iostatus
      character(len=*), parameter :: file_energy_scf = 'out_energy.txt'
      real(rp) :: en_f
      namelist /energy/ en_f

      open (unit_energy_scf, file=file_energy_scf, action='read', iostat=iostatus, status='old')
      read (unit_energy_scf, nml=energy, iostat=iostatus)

      close (unit_energy_scf)
   end subroutine get_fermi_scf

   real(rp) function delta_function(x)
      use precision_mod
    !
    ! --> 'fd': derivative of the Fermi-Dirac smearing. 0.5/(1.0+cosh(x))
    !
      real(rp), intent(in) :: x
    ! output: the value of the function
    ! input: the point where to compute the function
    ! local variable
      real(rp), parameter :: smearing = 0.02_rp  !smearing in eV
    ! Fermi-Dirac smearing
      if (abs(x) <= 36.0_rp) then
         delta_function = 1.0_rp/(smearing*(2.0_rp + exp(-x/smearing) + exp(+x/smearing)))
         ! in order to avoid problems for large values of x in the e
      else
         delta_function = 0.0_rp
      end if
   end function delta_function

end program bands
