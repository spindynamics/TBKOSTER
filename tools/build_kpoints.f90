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
!  build_kpoints.f90
!  TBKOSTER
program build_kpoints
   use precision_mod, only: rp
   implicit none
   integer, dimension(3) :: gx
   real(rp), dimension(:, :), allocatable :: x
   real(rp), dimension(:), allocatable :: w
   real(rp) :: w_cumul, v_factor, size
   integer :: ix, nx

   write (*, *) 'enter nkx, nky, nkz :'
   read (*, *) gx(1), gx(2), gx(3)

   write (*, *) 'enter v_factor (celldm) :'
   read (*, *) v_factor

   write (*, *) 'enter the size of the kmesh:'
   read (*, *) size

   nx = gx(1)*gx(2)*gx(3)
   write (*, *)
   write (*, *) 'number of kpoints : ', nx

   allocate (w(nx))
   w(:) = 1.0_rp/nx

   x = build_kmesh(gx)
   do ix = 1, nx
      x(ix, :) = size*x(ix, :)/v_factor
   end do

   call write_kmesh(x, w, nx)

contains

   function build_kmesh(gx) result(x)
      use precision_mod
      integer, dimension(3), intent(in) :: gx
      real(rp), dimension(:, :), allocatable :: x
      integer :: ix, igx1, igx2, igx3
      real(rp) :: x1, x2, x3

      ! if (size(x,1)/=product(gx)) then
      !    write(*,*) 'size of x incompatible'
      !    stop
      !endif
      allocate (x(product(gx), 3))

      ix = 1
      do igx1 = 1, gx(1)
         x1 = (2*igx1 - gx(1) - 1.0_rp)/(2.0_rp*gx(1))
         do igx2 = 1, gx(2)
            x2 = (2*igx2 - gx(2) - 1.0_rp)/(2.0_rp*gx(2))
            do igx3 = 1, gx(3)
               x3 = (2*igx3 - gx(3) - 1.0_rp)/(2.0_rp*gx(3))
               x(ix, :) = (/x1, x2, x3/)
               ix = ix + 1
            end do
         end do
      end do
   end function build_kmesh

   !> Write property (default: see source code) in text format to unit
   !> (default: 10), if it's a file its name is set to file (default:
   !> 'out_mesh.txt'), if tag (default: .true.) the namelist opening and closing
   !> tags are written
   subroutine write_kmesh(x, w, nx)
      use string_mod, only: int2str, real2str
      integer, intent(in):: nx
      real(rp), intent(in) :: x(nx, 3), w(nx)
      ! Local variables
      integer, parameter :: unit_mesh = 10
      integer :: ip, ix, ixs
      character(len=*), parameter :: file_kmesh = 'in_mesh.txt'
      character(len=9), parameter :: x_coord = 'cartesian'

      open (unit=unit_mesh, file=file_kmesh, action='write')
      write (unit_mesh, '(a)') '&mesh'
      write (unit_mesh, '(a)') ' type = ''list'''
      write (unit_mesh, '(a)') ' nx = '//int2str(nx)
      write (unit_mesh, '(a)') ' x_coord = '''//trim(x_coord)//''''
      do ix = 1, nx
         write (unit_mesh, '(a)') ' x('//int2str(ix)//',:) = ' &
            //real2str(x(ix, 1))//', ' &
            //real2str(x(ix, 2))//', ' &
            //real2str(x(ix, 3))
      end do
      do ix = 1, nx
         write (unit_mesh, '(a)') ' w('//int2str(ix)//') = '//real2str(w(ix))
      end do
      write (unit_mesh, '(a)') '/'
      close (unit_mesh)
   end subroutine write_kmesh

end program build_kpoints

