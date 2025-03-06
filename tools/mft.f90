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
! professionals having in-depth computer knowledge. users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The fact that you are presently reading this means that you have had
! knowledge of the CeciLL license and that you accept its terms.
!
!  mft.f90
!  TBKOSTER
program pmft
   use, intrinsic :: iso_fortran_env, ONLY: output_unit
   use precision_mod
   use string_mod
   implicit none
   integer, parameter :: unit_mft_in = 10
   integer, parameter :: unit_mft_out = 11
   integer, parameter :: unit_log = 12
   integer :: iostatus, iatom
   character(len=*), parameter :: dir = 'mft/'
   character(len=*), parameter :: file_mft_in = 'in_mft.txt'
   character(len=*), parameter :: file_mft_out = 'out_mft'
   character(len=*), parameter :: file_out_log = 'out_log.txt'
   integer :: ne, ns, na, na_mft, nangle, nxa, nconfig, iconfig
   integer, dimension(:), allocatable :: ia, ia2ie
   character(len=2), dimension(:), allocatable :: symbol
   character(len=7) :: calc
   character(len=4) :: type
   real(rp), dimension(3, 3) :: v
   real(rp), dimension(:, :), allocatable:: r
   real(rp) :: v_factor, mft_tot, mft_sum
   complex(rp), dimension(:), allocatable :: mft_s, mft_p, mft_px, mft_py, mft_pz
   complex(rp), dimension(:), allocatable :: mft_d, mft_dxy, mft_dyz, mft_dzx, mft_dx2y2, mft_dz2r2
   real(rp), dimension(:), allocatable :: mft_tot_i, mft_sum_i
   real(rp), dimension(:, :, :), allocatable::  mconfig_i
   real(rp), dimension(:, :), allocatable::  mconfig
   complex(rp), dimension(:, :), allocatable :: mft_s_i, mft_p_i, mft_px_i, mft_py_i, mft_pz_i
   complex(rp), dimension(:, :), allocatable :: mft_d_i, mft_dxy_i, mft_dyz_i, mft_dzx_i, &
                                                mft_dx2y2_i, mft_dz2r2_i
   namelist /element/ ne, symbol
   namelist /lattice/ v_factor, v
   namelist /atom/ ns, na, ia2ie, r
   namelist /mft/ calc, type, na_mft, ia, nangle, nxa
   namelist /mft_out/ mconfig, mft_s, mft_p, mft_px, mft_py, mft_pz, &
      mft_d, mft_dxy, mft_dyz, mft_dzx, mft_dx2y2, mft_dz2r2, mft_tot, mft_sum

   open (unit_log, file=dir//file_out_log, action='read', iostat=iostatus, status='old')
   allocate (symbol(0))
   read (unit_log, nml=element, iostat=iostatus)
   deallocate (symbol)
   allocate (symbol(ne))
   rewind (unit_log)
   read (unit_log, nml=element, iostat=iostatus)
   rewind (unit_log)
   read (unit_log, nml=lattice, iostat=iostatus)
   rewind (unit_log)
   allocate (r(0, 0), ia2ie(0))
   read (unit_log, nml=atom, iostat=iostatus)
   deallocate (ia2ie)
   allocate (ia2ie(na))
   rewind (unit_log)
   read (unit_log, nml=atom, iostat=iostatus)
   deallocate (r)
   allocate (r(na, 3))
   rewind (unit_log)
   read (unit_log, nml=atom, iostat=iostatus)
   close (unit_log)

   open (unit_mft_in, file=dir//file_mft_in, action='read', iostat=iostatus, status='old')
   allocate (ia(0))
   read (unit_mft_in, nml=mft, iostat=iostatus)
   deallocate (ia)
   allocate (ia(na_mft))
   if (na_mft == na) then
      do iatom = 1, na
         ia(iatom) = iatom
      end do
   else
      rewind (unit_mft_in)
      read (unit_mft_in, nml=mft, iostat=iostatus)
   end if

   rewind (unit_mft_in)
   read (unit_mft_in, nml=mft, iostat=iostatus)
   close (unit_mft_in)

   if (calc == 'mae') then
      if (type == 'list') then
         nconfig = nxa
      elseif (type == 'mesh') then
         nconfig = nxa**2
      elseif (type == 'path') then
         nconfig = (nangle - 1)*nxa + 1
      end if
   elseif (calc == 'mconfig') then
      nconfig = nxa
   end if
   allocate (mconfig(na, 2), mconfig_i(na, nconfig, 2))
   allocate (mft_tot_i(nconfig), mft_sum_i(nconfig), mft_s_i(nconfig, na_mft), mft_p_i(nconfig, na_mft), &
             mft_px_i(nconfig, na_mft), mft_py_i(nconfig, na_mft), mft_pz_i(nconfig, na_mft))
   allocate (mft_d_i(nconfig, na_mft), mft_dxy_i(nconfig, na_mft), mft_dyz_i(nconfig, na_mft), &
             mft_dzx_i(nconfig, na_mft), mft_dx2y2_i(nconfig, na_mft), mft_dz2r2_i(nconfig, na_mft))

   allocate (mft_s(na_mft), mft_p(na_mft), mft_px(na_mft), mft_py(na_mft), mft_pz(na_mft))
   allocate (mft_d(na_mft), mft_dxy(na_mft), mft_dyz(na_mft), mft_dzx(na_mft), mft_dx2y2(na_mft), mft_dz2r2(na_mft))
   do iconfig = 1, nconfig
      open (unit_mft_out, file=dir//file_mft_out//'_'//int2str(iconfig)//'.txt', action='read', iostat=iostatus, status='old')
      READ (unit_mft_out, nml=mft_out, iostat=iostatus)
      mconfig_i(:, iconfig, :) = mconfig(:, :)
      mft_tot_i(iconfig) = mft_tot
      mft_sum_i(iconfig) = mft_sum
      mft_s_i(iconfig, :) = mft_s(:)
      mft_p_i(iconfig, :) = mft_p(:)
      mft_d_i(iconfig, :) = mft_d(:)
      mft_px_i(iconfig, :) = mft_px(:)
      mft_py_i(iconfig, :) = mft_py(:)
      mft_pz_i(iconfig, :) = mft_pz(:)
      mft_dxy_i(iconfig, :) = mft_dxy(:)
      mft_dyz_i(iconfig, :) = mft_dyz(:)
      mft_dzx_i(iconfig, :) = mft_dzx(:)
      mft_dx2y2_i(iconfig, :) = mft_dx2y2(:)
      mft_dz2r2_i(iconfig, :) = mft_dz2r2(:)
      close (unit_mft_out)
   end do

   call build_mft(calc, type, na, ne, ia2ie, Symbol, v, r, na_mft, ia, nconfig, mconfig_i, &
                  mft_tot_i, mft_sum_i, mft_s_i, mft_p_i, mft_d_i, &
                  mft_px_i, mft_py_i, mft_pz_i, mft_dxy_i, mft_dyz_i, mft_dzx_i, mft_dx2y2_i, mft_dz2r2_i)

contains

  subroutine build_mft(calc, type, na, ne, ia2ie, Symbol, v, r, na_mft, ia, nconfig, mconfig_i, &
                      mft_tot_i, mft_sum_i, mft_s_i, mft_p_i, mft_d_i, &
                      mft_px_i, mft_py_i, mft_pz_i, mft_dxy_i, mft_dyz_i, mft_dzx_i, mft_dx2y2_i, mft_dz2r2_i)
    use precision_mod
    use string_mod
    use math_mod, ONLY: i_unit
    implicit none
    integer :: na, ne, nconfig, na_mft
    integer :: ia_mft, iconfig, iatom, imin
    integer, parameter :: unit_mae_angle = 10
    integer, parameter :: unit_mae_xyz = 11
    integer, parameter :: unit_mae_atom = 12
    integer, parameter :: unit_config = 13
    integer :: ia(na_mft), ia2ie(na)
    character(len=2), dimension(ne):: symbol
    character(len=7) :: calc
    character(len=4) :: type
    real(rp), dimension(na, nconfig, 2) :: mconfig_i
    real(rp), dimension(na, 3) :: r
    real(rp), dimension(3, 3) :: v
    real(rp) :: angle, mae, theta, phi, emin
    real(rp) :: mft_tot_i(nconfig), mft_sum_i(nconfig)
    complex(rp) :: mft_s_i(nconfig, na_mft), mft_p_i(nconfig, na_mft), &
                    mft_px_i(nconfig, na_mft), mft_py_i(nconfig, na_mft), mft_pz_i(nconfig, na_mft)
    complex(rp) :: mft_d_i(nconfig, na_mft), mft_dxy_i(nconfig, na_mft), mft_dyz_i(nconfig, na_mft), &
                    mft_dzx_i(nconfig, na_mft), mft_dx2y2_i(nconfig, na_mft), mft_dz2r2_i(nconfig, na_mft)
    character(len=*), parameter :: dir = 'mft/'
    character(len=*), parameter :: file_mae_theta = 'mae_angle.dat'
    character(len=*), parameter :: file_mae_xyz = 'mae.xyz'
    character(len=*), parameter :: file_mae_atom = 'mae_atom.dat'
    character(len=*), parameter :: file_config = 'mft_config.dat'
    character(len=80) :: fmt

    if (calc == 'mae') then

      if (type == 'list' .or. type == 'path') then
        open (unit=unit_mae_atom, file=dir//file_mae_atom, action='write')
        write (unit_mae_atom, *) '@# i mae(meV)'

        write (unit_mae_atom, *) '@# site mae s'
        do ia_mft = 1, na_mft
            write (unit_mae_atom, '(2(a,1X))') int2str(ia(ia_mft)), &
              real2str(1000*REAL(mft_s_i(nconfig, ia_mft) - mft_s_i(1, ia_mft)))
        end do

        write (unit_mae_atom, *) '@# site mae p'
        do ia_mft = 1, na_mft
            write (unit_mae_atom, '(5(a,1X))') int2str(ia(ia_mft)), &
              real2str(1000*REAL(mft_p_i(nconfig, ia_mft) - mft_p_i(1, ia_mft))), &
              real2str(1000*REAL(mft_px_i(nconfig, ia_mft) - mft_px_i(1, ia_mft))), &
              real2str(1000*REAL(mft_py_i(nconfig, ia_mft) - mft_py_i(1, ia_mft))), &
              real2str(1000*REAL(mft_pz_i(nconfig, ia_mft) - mft_pz_i(1, ia_mft)))
        end do
        write (unit_mae_atom, *) '@# site mae d'
        do ia_mft = 1, na_mft
            write (unit_mae_atom, '(7(a,1X))') int2str(ia(ia_mft)), &
              real2str(1000*REAL(mft_d_i(nconfig, ia_mft) - mft_d_i(1, ia_mft))), &
              real2str(1000*REAL(mft_dxy_i(nconfig, ia_mft) - mft_dxy_i(1, ia_mft))), &
              real2str(1000*REAL(mft_dyz_i(nconfig, ia_mft) - mft_dyz_i(1, ia_mft))), &
              real2str(1000*REAL(mft_dzx_i(nconfig, ia_mft) - mft_dzx_i(1, ia_mft))), &
              real2str(1000*REAL(mft_dx2y2_i(nconfig, ia_mft) - mft_dx2y2_i(1, ia_mft))), &
              real2str(1000*REAL(mft_dz2r2_i(nconfig, ia_mft) - mft_dz2r2_i(1, ia_mft)))
        end do

        close (unit_mae_atom)

        open (unit=unit_mae_angle, file=dir//file_mae_theta, action='write')
        Emin = minval(mft_tot_i, dim=1)
        imin = minloc(mft_tot_i, dim=1)
        angle = 0.0_rp
        mae = 0.0_rp
        !mae=1000*(mft_tot_i(1)-Emin)
        write (unit_mae_angle, *) '@# angle  MAE(meV)'
        write (unit_mae_angle, *) '@# theta= ', real2str(mconfig_i(1, imin, 1)), ' phi= ', real2str(mconfig_i(1, imin, 2))
        write (unit_mae_angle, *) '@# Emin= ', real2str(Emin)
        write (unit_mae_angle, '(3(a,1X))') real2str(angle), real2str(mae)
        do iconfig = 1, nconfig - 1
            angle = angle + &
                    sqrt((mconfig_i(1, iconfig + 1, 1) - mconfig_i(1, iconfig, 1))**2 &
                        + (mconfig_i(1, iconfig + 1, 2) - mconfig_i(1, iconfig, 2))**2)
            mae = 1000*(mft_tot_i(iconfig) - mft_tot_i(1))
            !mae=1000*(mft_tot_i(iconfig)-Emin)
            write (unit_mae_angle, '(3(a,1X))') real2str(angle), real2str(mae)
        end do

        if (na_mft == na) then
            open (unit=unit_mae_xyz, file=dir//file_mae_xyz, action='write')
            ! Number of atoms
            write (unit_mae_xyz, '(a)') int2str(na)
            ! Serie of key/value pairs
            write (unit_mae_xyz, '(a)') 'Lattice="' &
              //real2str(v(1, 1))//' ' &
              //real2str(v(1, 2))//' ' &
              //real2str(v(1, 3))//' ' &
              //real2str(v(2, 1))//' ' &
              //real2str(v(2, 2))//' ' &
              //real2str(v(2, 3))//' ' &
              //real2str(v(3, 1))//' ' &
              //real2str(v(3, 2))//' ' &
              //real2str(v(3, 3))//'" ' &
              //'Properties=species:S:1:pos:R:3:kinetic_energy:R:1'
            ! Atoms
            do iatom = 1, na
              write (unit_mae_xyz, '(a)') symbol(ia2ie(iatom))//' ' &
                  //real2str(r(iatom, 1))//' ' &
                  //real2str(r(iatom, 2))//' ' &
                  //real2str(r(iatom, 3))//' ' &
                  //real2str(1000*real(mft_s_i(nconfig, iatom) + mft_p_i(nconfig, iatom) + mft_d_i(nconfig, iatom) - &
                                      (mft_s_i(1, iatom) + mft_p_i(1, iatom) + mft_d_i(1, iatom))))//' '
            end do
            close (unit_mae_xyz)
        end if
      elseif (type == 'mesh') then
        open (unit=unit_mae_angle, file=dir//file_mae_theta, action='write')
        Emin = minval(mft_tot_i, dim=1)
        imin = minloc(mft_tot_i, dim=1)
        write (unit_mae_angle, *) '@# theta  phi  Energy(eV)'
        write (unit_mae_angle, *) '@# theta= ', real2str(mconfig_i(1, imin, 1)), ' phi= ', real2str(mconfig_i(1, imin, 2))
        write (unit_mae_angle, *) '@# Emin= ', real2str(Emin)
        do iconfig = 1, nconfig
            theta = mconfig_i(1, iconfig, 1)
            phi = mconfig_i(1, iconfig, 2)
            !mae=1000*(mft_tot_i(iconfig)-mft_tot_i(1))
            mae = 1000*(mft_tot_i(iconfig) - Emin)
            write (unit_mae_angle, '(3(a,1X))') real2str(theta), real2str(phi), real2str(mae)
        end do
      end if

    elseif (calc == 'mconfig') then
        open (unit=unit_config, file=dir//file_config, action='write')
        write (unit_config, *) '@# config  Energy(eV)'
        do iconfig = 1, nconfig
          write (unit_config, '(3(a,1X))') int2str(iconfig), real2str(mft_tot_i(iconfig))
        end do

    end if

  end subroutine build_mft
end program pmft

