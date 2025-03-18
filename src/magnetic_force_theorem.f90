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
!  magnetic_force_theorem.f90
!  TBKOSTER
module magnetic_force_theorem_mod
   use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
   use atom_mod
   use element_mod
   use energy_mod
   use hamiltonian_tb_mod
   use math_mod, only: i_unit, delta_function, fermi_function, theta_function, deg2rad, rad2deg
   use mesh_mod
#if defined(OpenMP_Fortran_FOUND)
   use omp_lib
#endif
   use precision_mod, only: rp
   use string_mod, only: TBKOSTER_flush, int2str, lower, real2str, sl, cmplx2str
   use units_mod
   implicit none
   private

   !> Derived type properties for i/o methods
   character(len=sl), dimension(*), parameter :: property_list = &
                                                 [character(len=sl) :: &
                                                  'na_mft', &
                                                  'ia', &
                                                  'Eref', &
                                                  'mconfig' &
                                                  ]

   type, public :: magnetic_force_theorem
      ! Units
      class(units), pointer :: u
      ! Elements
      class(element), pointer :: e
      ! Atom
      class(atom), pointer :: a
      ! Reciprocal space mesh
      class(mesh), pointer :: k
      ! Hamiltonian
      class(hamiltonian_tb), pointer :: h
      ! Energy
      class(energy), pointer :: en
      ! number of local mft atomic site number
      integer :: na_mft
      ! Local mft atomic site index
      integer, dimension(:), allocatable :: ia
      ! number of spin angles
      integer :: nangle, nconfig
      ! mft angles in case of mae calculation
      real(rp), dimension(:, :), allocatable :: angle
      ! atomic magnetic configurations
      real(rp), dimension(:, :, :), allocatable :: mconfig
      ! reference energy (default value is Ef)
      real(rp) :: Eref
      ! Total mft
      real(rp) :: mft_tot
      ! Local mft
      complex(rp), dimension(:), allocatable :: mft_s
      complex(rp), dimension(:), allocatable :: mft_p, mft_px, mft_py, mft_pz
      complex(rp), dimension(:), allocatable :: mft_d, mft_dxy, mft_dyz, mft_dzx, &
                                                mft_dx2y2, mft_dz2r2

   contains
      ! Destructor
      final :: destructor
      ! Procedures
      procedure :: add_mft_k
      procedure :: add_mft_local_k
      procedure :: initialize
      procedure :: read_txt
      procedure :: write_txt
      procedure :: write_txt_formatted
   end type magnetic_force_theorem

   ! Constructor
   interface magnetic_force_theorem
      procedure :: constructor
   end interface magnetic_force_theorem

contains
   function constructor(en) result(obj)
      class(energy), target, intent(in) :: en
      type(magnetic_force_theorem) :: obj

      obj%u => en%u
      obj%e => en%e
      obj%a => en%a
      obj%k => en%k
      obj%h => en%h
      obj%en => en
   end function constructor

   subroutine destructor(obj)
      type(magnetic_force_theorem) :: obj

      if (allocated(obj%ia)) deallocate (obj%ia)
      if (allocated(obj%mft_s)) deallocate (obj%mft_s)
      if (allocated(obj%mft_p)) deallocate (obj%mft_p)
      if (allocated(obj%mft_px)) deallocate (obj%mft_px)
      if (allocated(obj%mft_py)) deallocate (obj%mft_py)
      if (allocated(obj%mft_pz)) deallocate (obj%mft_pz)
      if (allocated(obj%mft_d)) deallocate (obj%mft_d)
      if (allocated(obj%mft_dxy)) deallocate (obj%mft_dxy)
      if (allocated(obj%mft_dyz)) deallocate (obj%mft_dyz)
      if (allocated(obj%mft_dzx)) deallocate (obj%mft_dzx)
      if (allocated(obj%mft_dx2y2)) deallocate (obj%mft_dx2y2)
      if (allocated(obj%mft_dz2r2)) deallocate (obj%mft_dz2r2)
   end subroutine destructor

   subroutine add_mft_k(obj, ik, isl, Eref)
      use, intrinsic :: iso_fortran_env, only: output_unit
      implicit none
      ! INPUT
      class(magnetic_force_theorem), intent(inout) :: obj
      integer, intent(in) :: ik, isl
      real(rp), intent(in) :: Eref
      ! LOCAL
      real(rp) :: en, ff
      integer :: ien, jmat, jmat2, jj, nn

      select case (obj%a%ns)
      case (1, 2)
         write (*, *) 'for mft ns should be 4'
         stop
      case (4)
         do jmat = 1, obj%h%nh
            nn = ik + (jmat - 1)*obj%k%nx
            ff = theta_function(-(obj%en%en_k(nn) - Eref)/obj%en%degauss, obj%en%smearing) &
                 *(obj%en%en_k(nn) - Eref)
            obj%mft_tot = obj%mft_tot + ff*obj%k%w(ik)
         end do
      end select
   end subroutine add_mft_k

   subroutine add_mft_local_k(obj, ik, isl, v_k, Eref)
      ! INPUT
      class(magnetic_force_theorem), intent(inout) :: obj
      integer, intent(in) :: ik, isl
      real(rp), intent(in) :: Eref
      complex(rp), dimension(2, obj%h%nh, obj%h%nh), intent(in) :: v_k
      !        LOCAL
      real(rp) :: ff, en
      complex(rp) :: w_s, w_p, w_d
      complex(rp) :: w_px, w_py, w_pz
      complex(rp) :: w_dxy, w_dyz, w_dzx, w_dx2y2, w_dz2r2
      integer :: ia, ie, io, ia_mft, ispin, jspin, imat, jmat, jmat2, jj, nn
      integer :: ispin_up, ispin_dn, imat_ispin, imat_jspin

      if (obj%na_mft > 0) then
         select case (obj%a%ns)
         case (1, 2)
            write (*, *) 'for mft ns should be 4'
            stop
         case (4)

            ispin_up = 1
            ispin_dn = 2

            ! 1=s; 2=px; 3=py; 4=pz; 5=dxy; 6=dyz; 7=dzx; 8=dx^2-y^2; 9=d(3z^2-r^2)

            do ia_mft = 1, obj%na_mft
               ia = obj%ia(ia_mft)
               ie = obj%a%ia2ie(ia)

               do jmat = 1, obj%h%nh

                  nn = ik + (jmat - 1)*obj%k%nx

                  ff = theta_function(-(obj%en%en_k(nn) - Eref)/obj%en%degauss, obj%en%smearing) &
                       *(obj%en%en_k(nn) - Eref)

                  do io = 1, obj%e%no(ie)
                     do ispin = 1, 2

                        imat_ispin = obj%h%iaos2ih(ia, io, ispin)

                        !do jspin=1,2
                        jspin = ispin !diagonal elements
                        imat_jspin = obj%h%iaos2ih(ia, io, jspin)

                        if (obj%e%o(ie, io) == 1) then
                           w_s = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                  + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                           obj%mft_s(ia_mft) = obj%mft_s(ia_mft) + w_s*ff*obj%k%w(ik)

                        elseif (obj%e%o(ie, io) >= 2 .and. obj%e%o(ie, io) <= 4) then
                           w_p = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                  + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                           obj%mft_p(ia_mft) = obj%mft_p(ia_mft) + w_p*ff*obj%k%w(ik)

                           if (obj%e%o(ie, io) == 2) then
                              w_px = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                      + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                              obj%mft_px(ia_mft) = obj%mft_px(ia_mft) + w_px*ff*obj%k%w(ik)
                           elseif (obj%e%o(ie, io) == 3) then
                              w_py = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                      + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                              obj%mft_py(ia_mft) = obj%mft_py(ia_mft) + w_py*ff*obj%k%w(ik)
                           elseif (obj%e%o(ie, io) == 4) then
                              w_pz = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                      + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                              obj%mft_pz(ia_mft) = obj%mft_pz(ia_mft) + w_pz*ff*obj%k%w(ik)
                           end if

                        elseif (obj%e%o(ie, io) >= 5 .and. obj%e%o(ie, io) <= 9) then
                           w_d = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                  + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                           obj%mft_d(ia_mft) = obj%mft_d(ia_mft) + w_d*ff*obj%k%w(ik)

                           if (obj%e%o(ie, io) == 5) then
                              w_dxy = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                       + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                              obj%mft_dxy(ia_mft) = obj%mft_dxy(ia_mft) + w_dxy*ff*obj%k%w(ik)
                           elseif (obj%e%o(ie, io) == 6) then
                              w_dyz = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                       + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                              obj%mft_dyz(ia_mft) = obj%mft_dyz(ia_mft) + w_dyz*ff*obj%k%w(ik)
                           elseif (obj%e%o(ie, io) == 7) then
                              w_dzx = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                       + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                              obj%mft_dzx(ia_mft) = obj%mft_dzx(ia_mft) + w_dzx*ff*obj%k%w(ik)
                           elseif (obj%e%o(ie, io) == 8) then
                              w_dx2y2 = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                         + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                              obj%mft_dx2y2(ia_mft) = obj%mft_dx2y2(ia_mft) + w_dx2y2*ff*obj%k%w(ik)
                           elseif (obj%e%o(ie, io) == 9) then
                              w_dz2r2 = (v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)) &
                                         + v_k(2, imat_jspin, jmat)*conjg(v_k(1, imat_ispin, jmat)))/2
                              obj%mft_dz2r2(ia_mft) = obj%mft_dz2r2(ia_mft) + w_dz2r2*ff*obj%k%w(ik)
                           end if

                        end if

                        !end do ! fin de la boucle sur jspin
                     end do ! fin de la boucle sur ispin
                  end do ! fin de la boucle sur io
               end do ! fin de la boucle sur jmat
            end do  ! fin de la boucle sur ia_dos

         end select

      end if
   end subroutine add_mft_local_k

   subroutine initialize(obj)
      class(magnetic_force_theorem), intent(inout) :: obj

      obj%mft_tot = 0.0_rp
      ! Local mft
      if (obj%na_mft > 0) then
         if (allocated(obj%mft_s)) deallocate (obj%mft_s)
         allocate (obj%mft_s(obj%na_mft))
         if (allocated(obj%mft_p)) deallocate (obj%mft_p)
         allocate (obj%mft_p(obj%na_mft))
         if (allocated(obj%mft_px)) deallocate (obj%mft_px)
         allocate (obj%mft_px(obj%na_mft))
         if (allocated(obj%mft_py)) deallocate (obj%mft_py)
         allocate (obj%mft_py(obj%na_mft))
         if (allocated(obj%mft_pz)) deallocate (obj%mft_pz)
         allocate (obj%mft_pz(obj%na_mft))
         if (allocated(obj%mft_d)) deallocate (obj%mft_d)
         allocate (obj%mft_d(obj%na_mft))
         if (allocated(obj%mft_dxy)) deallocate (obj%mft_dxy)
         allocate (obj%mft_dxy(obj%na_mft))
         if (allocated(obj%mft_dyz)) deallocate (obj%mft_dyz)
         allocate (obj%mft_dyz(obj%na_mft))
         if (allocated(obj%mft_dzx)) deallocate (obj%mft_dzx)
         allocate (obj%mft_dzx(obj%na_mft))
         if (allocated(obj%mft_dx2y2)) deallocate (obj%mft_dx2y2)
         allocate (obj%mft_dx2y2(obj%na_mft))
         if (allocated(obj%mft_dz2r2)) deallocate (obj%mft_dz2r2)
         allocate (obj%mft_dz2r2(obj%na_mft))

         obj%mft_s = 0.0_rp
         obj%mft_p = 0.0_rp
         obj%mft_px = 0.0_rp
         obj%mft_py = 0.0_rp
         obj%mft_pz = 0.0_rp
         obj%mft_d = 0.0_rp
         obj%mft_dxy = 0.0_rp
         obj%mft_dyz = 0.0_rp
         obj%mft_dzx = 0.0_rp
         obj%mft_dx2y2 = 0.0_rp
         obj%mft_dz2r2 = 0.0_rp
      end if
   end subroutine initialize

   !> Check the validity of mesh type
   subroutine check_calc(type)
      character(len=*) :: type
      if (type /= 'mae' &
          .and. type /= 'mconfig') then
         write (error_unit, *) 'mft%check_calc(): mft%calc must be one of: &
   &       ''mae'', ''mconfig'', '
         error stop
      end if
   end subroutine check_calc

   !> Check the validity of mesh type
   subroutine check_type(type)
      character(len=*) :: type

      if (type /= 'list' &
          .and. type /= 'mesh' &
          .and. type /= 'path') then
         write (error_unit, *) 'mft%check_type(): mft%type must be one of: &
   &       ''list'', ''mesh'', ''path'''
         error stop
      end if
   end subroutine check_type

   function build_path(nxa, angle_xs) result(x)
      integer, intent(in) :: nxa
      real(rp), dimension(:, :), intent(in) :: angle_xs
      real(rp), dimension((size(angle_xs, 1) - 1)*nxa + 1, 2) :: x
      integer :: ix, ixa, nxs, ixs
      nxs = size(angle_xs, 1)
      ix = 1
      do ixs = 1, nxs - 1
         do ixa = 0, nxa - 1
            x(ix, :) = angle_xs(ixs, :) + (angle_xs(ixs + 1, :) - angle_xs(ixs, :))*ixa/nxa
            ix = ix + 1
         end do
      end do
      x(ix, :) = angle_xs(nxs, :)
   end function build_path

   function build_mesh(nxa) result(x)
      integer, intent(in) :: nxa
      real(rp), dimension(nxa**2, 2) :: x
      real(rp) :: theta, phi
      real(rp) :: theta_max, phi_max
      integer :: ix, ixa, jxa
      theta_max = 180_rp
      phi_max = 360_rp
      ix = 1
      do ixa = 0, nxa - 1
      do jxa = 0, nxa - 1
         theta = ixa*theta_max/(nxa - 1)
         phi = jxa*phi_max/(nxa - 1)
         x(ix, :) = (/theta, phi/)
         ix = ix + 1
      end do
      end do
   end function build_mesh

   !> Read object in text format from file (default: 'in_dos.txt')
   subroutine read_txt(obj, file)
      class(magnetic_force_theorem), intent(inout) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable :: file_rt
      integer :: iostatus
      logical :: isopen
      ! Namelist variables
      character(len=7) :: calc
      character(len=4) :: type = 'list'
      integer  :: na_mft, iangle, nangle, nxa, iconfig, nconfig
      integer  :: iia
      integer, dimension(:), allocatable :: ia
      real(rp) :: Eref
      real(rp), dimension(:, :), allocatable  :: angle_xs
      real(rp), dimension(:, :), allocatable  :: angle
      real(rp), dimension(:, :, :), allocatable  :: mconfig
      real(rp) :: en_f, en
      ! Namelist
      namelist /mft/ Eref, calc, type, na_mft, ia, nangle, nxa, angle_xs, mconfig
      namelist /energy/ en_f, en

      if (present(file)) then
         file_rt = trim(file)
      else
         file_rt = 'in_mft.txt'
      end if
      inquire (unit=10, opened=isopen)
      if (isopen) then
         write (error_unit, '(a)') 'mft%read_txt() : Unit 10 is already open'
         error stop
      else
         open (unit=10, file=file_rt, action='read', iostat=iostatus, status='old')
      end if
      if (iostatus /= 0) then
         write (error_unit, *) 'mft%read_txt(): file ', file_rt, ' not found'
         error stop
      end if

      na_mft = 0
      open (unit=11, file='out_energy.txt', action='read')
      read (unit=11, nml=energy)
      Eref = en_f
      close (unit=11)
      read (10, nml=mft, iostat=iostatus)
      calc = lower(calc)
      call check_calc(calc)
      type = lower(type)
      call check_type(type)
      rewind (10)

      allocate (ia(0))
      read (10, nml=mft, iostat=iostatus)
      deallocate (ia)
      allocate (ia(na_mft))
      if (na_mft > 0 .and. na_mft < obj%a%na) then
         rewind (10)
         read (10, nml=mft)
      elseif (na_mft == obj%a%na) then
         do iia = 1, obj%a%na
            ia(iia) = iia
         end do
      end if
      allocate (angle(0, 0))
      allocate (angle_xs(0, 0))
      select case (calc)
      case ('mae')
         select case (type)
         case ('list')
            read (10, nml=mft, iostat=iostatus)
            deallocate (angle)
            deallocate (angle_xs)
            rewind (10)
            read (10, nml=mft, iostat=iostatus)
            nconfig = nxa
            nangle = nxa
            allocate (angle_xs(nangle, 2))
            rewind (10)
            read (10, nml=mft, iostat=iostatus)
            allocate (mconfig(obj%a%na, nconfig, 2))
            do iia = 1, obj%a%na
               mconfig(iia, :, :) = angle_xs(:, :)
            end do
         case ('mesh')
            read (10, nml=mft, iostat=iostatus)
            deallocate (angle)
            rewind (10)
            read (10, nml=mft, iostat=iostatus)
            nconfig = nxa**2
            allocate (angle(nconfig, 2))
            angle = build_mesh(nxa)
            allocate (mconfig(obj%a%na, nconfig, 2))
            do iia = 1, obj%a%na
               mconfig(iia, :, :) = angle(:, :)
            end do
         case ('path')
            read (10, nml=mft, iostat=iostatus)
            deallocate (angle)
            deallocate (angle_xs)
            rewind (10)
            read (10, nml=mft, iostat=iostatus)
            nconfig = (nangle - 1)*nxa + 1
            allocate (angle(nconfig, 2))
            allocate (angle_xs(nangle, 2))
            rewind (10)
            read (10, nml=mft, iostat=iostatus)
            angle = build_path(nxa, angle_xs)
            allocate (mconfig(obj%a%na, nconfig, 2))
            do iia = 1, obj%a%na
               mconfig(iia, :, :) = angle(:, :)
            end do
         end select
      case ('mconfig')
         allocate (mconfig(0, 0, 0))
         read (10, nml=mft, iostat=iostatus)
         rewind (10)
         read (10, nml=mft, iostat=iostatus)
         nconfig = nxa
         deallocate (mconfig)
         rewind (10)
         allocate (mconfig(obj%a%na, nconfig, 2))
         read (10, nml=mft, iostat=iostatus)
      end select

      do iangle = 1, nangle
         angle(iangle, :) = angle(iangle, :)*deg2rad
      end do

      do iia = 1, obj%a%na
      do iconfig = 1, nconfig
         mconfig(iia, iconfig, :) = mconfig(iia, iconfig, :)*deg2rad
      end do
      end do

      obj%na_mft = na_mft
      obj%nangle = nangle
      obj%nconfig = nconfig
      !obj%nxa=nxa
      obj%angle = angle
      obj%mconfig = mconfig
      obj%Eref = Eref*obj%u%convert_energy('to', 'hau')
      call move_alloc(ia, obj%ia)
      call move_alloc(angle, obj%angle)
      call move_alloc(mconfig, obj%mconfig)
      call obj%initialize()
      close (unit=10)
      !deallocate(file_rt)
   end subroutine read_txt

   !> Write object in text format to unit (default: 10), if it's a file
   !> its name is set to file (default: 'out_dos.txt')
   subroutine write_txt(obj, file, unit)
      class(magnetic_force_theorem), intent(in) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable         :: file_rt
      integer, intent(in), optional :: unit
      integer                     :: unit_rt
      ! Namelist variables
      character(len=6) :: calc
      character(len=4) :: type
      integer  :: na_mft, nangle, nxa
      integer  :: iia, iangle
      integer, dimension(:), allocatable :: ia
      real(rp) :: Eref
      real(rp), dimension(:, :), allocatable  :: angle_xs
      real(rp), dimension(:, :, :), allocatable  :: mconfig
      ! Namelist
      namelist /mft/ Eref, calc, type, na_mft, ia, nangle, nxa, angle_xs, mconfig

      if (present(file)) then
         file_rt = file
      else
         file_rt = 'out_mft.txt'
      end if
      if (present(unit)) then
         unit_rt = unit
      else
         unit_rt = 10
      end if

      if (.not. present(unit)) then
         open (unit=unit_rt, file=file_rt, action='write')
      end if

      na_mft = obj%na_mft
      ia = obj%ia
      nangle = obj%nangle
      Eref = obj%Eref

      write (unit_rt, nml=mft)
      call TBKOSTER_flush(unit_rt)

      if (.not. present(unit)) then
         close (unit_rt)
      end if
      !deallocate(file_rt)
   end subroutine write_txt

   !> Write property (default: property_list) in text format to unit
   !> (default: 10), if it's a file its name is set to file (default:
   !> 'out_dos.txt'), if tag (default: .true.) the namelist opening and closing
   !> tags are written
   subroutine write_txt_formatted(obj, file, property, tag, config, unit)
      class(magnetic_force_theorem), intent(in) :: obj
      character(len=*), intent(in), optional :: file
      character(len=:), allocatable         :: file_rt
      character(len=*), dimension(:), intent(in), optional :: property
      character(len=:), dimension(:), allocatable         :: property_rt
      logical, intent(in), optional :: tag
      logical                     :: tag_rt
      integer, intent(in), optional :: unit, config
      integer                     :: unit_rt
      ! Local variables
      integer :: ia_mft, iangle, ip, ia, iconfig
      real(rp) :: mft_sum

      if (present(file)) then
         file_rt = file
      else
         file_rt = 'out_mft.txt'
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
      if (present(config)) then
         iconfig = config
      else
         iconfig = 0
      end if

      if (.not. present(unit)) then
         open (unit=unit_rt, file=file_rt, action='write')
      end if
      if (tag_rt) then
         write (unit_rt, '(a)') '&mft_out'
      end if

      do ip = 1, size(property_rt)
         select case (lower(trim(property_rt(ip))))
         case ('Eref')
            write (unit_rt, '(a)') ' Eref = '//real2str(obj%Eref)
         case ('mconfig')
            if (iconfig /= 0) then
               ! write(unit_rt,'(a)') ' nconfig = ' // int2str(obj%nconfig)
               !write(unit_rt,'(a)') ' iconfig = ' // int2str(iconfig)
               do ia = 1, obj%a%na
                  write (unit_rt, '(a)') ' mconfig('//int2str(ia)//',:)= ' &
                     //real2str(obj%mconfig(ia, iconfig, 1)*rad2deg)//' , ' &
                     //real2str(obj%mconfig(ia, iconfig, 2)*rad2deg)
               end do
            else
               !write(unit_rt,'(a)') ' nconfig = ' // int2str(obj%nconfig)
            end if

         case ('na_mft')
            write (unit_rt, '(a)') ' na_mft = '//int2str(obj%na_mft)
         case ('ia')
            do ia_mft = 1, obj%na_mft
               write (unit_rt, '(a)') ' ia('//int2str(ia_mft)//') = ' &
                  //int2str(obj%ia(ia_mft))
            end do
         case ('mft')
            mft_sum = 0.0_rp

            do ia_mft = 1, obj%na_mft

               write (unit_rt, '(a)') ' mft_s('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_s(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_p('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_p(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_px('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_px(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_py('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_py(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_pz('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_pz(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_d('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_d(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_dxy('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_dxy(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_dyz('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_dyz(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_dzx('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_dzx(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_dx2y2('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_dx2y2(ia_mft)*obj%u%convert_energy('from', 'hau'))

               write (unit_rt, '(a)') ' mft_dz2r2('//int2str(ia_mft)//') = ' &
                  //cmplx2str(obj%mft_dz2r2(ia_mft)*obj%u%convert_energy('from', 'hau'))

               mft_sum = mft_sum + obj%mft_s(ia_mft) + obj%mft_p(ia_mft) + obj%mft_d(ia_mft)
            end do
            write (unit_rt, '(a)') ' mft_sum= '//real2str(mft_sum*obj%u%convert_energy('from', 'hau'))
            write (unit_rt, '(a)') ' mft_tot= '//real2str(obj%mft_tot*obj%u%convert_energy('from', 'hau'))
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
end module magnetic_force_theorem_mod
