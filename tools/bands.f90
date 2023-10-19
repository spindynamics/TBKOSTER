!
! Copyright (C) 2019
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
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
!  bands.f90
!  DyNaMol
PROGRAM bands
    USE, INTRINSIC :: iso_fortran_env, ONLY: output_unit
    USE precision_mod
    USE string_mod
    IMPLICIT NONE
    INTEGER, PARAMETER :: unit_mesh = 11
    INTEGER, PARAMETER :: unit_hamiltonian_tb = 12
    INTEGER, PARAMETER :: unit_band_in = 13 
    INTEGER, PARAMETER :: unit_band_out = 14
    INTEGER :: iostatus, nx, nh, ns, nsl, no_max, na_band
    logical :: isopen
    CHARACTER(len=*), PARAMETER :: dir = 'band/'
    CHARACTER(len=*), PARAMETER :: file_band_in = 'in_band.txt' 
    CHARACTER(len=*), PARAMETER :: file_band_out = 'out_band.txt'
    CHARACTER(len=*), PARAMETER :: file_mesh = 'out_mesh.txt'
    CHARACTER(len=*), PARAMETER :: file_hamiltonian_tb = 'out_hamiltonian_tb.txt'
    CHARACTER(len=9) :: x_coord
    CHARACTER(len=4) :: TYPE
    CHARACTER(len=80) :: line
    INTEGER, DIMENSION(:), ALLOCATABLE :: ia_band,iband2io
    REAL(rp), DIMENSION(:, :, :), ALLOCATABLE :: en_k
    REAL(rp), DIMENSION(:), ALLOCATABLE :: w
    REAL(rp), DIMENSION(:, :), ALLOCATABLE :: x
    REAL(rp), DIMENSION(:,:,:,:,:), allocatable :: w_en_band_local
    REAL(rp) :: en_min, en_max, en_f
    NAMELIST /mesh/ TYPE, nx, x_coord, x, w
    NAMELIST /hamiltonian_tb/ nh, ns
    NAMELIST /band/ na_band,ia_band 
    NAMELIST /band_out/en_k,iband2io,w_en_band_local

   ! inquire(unit=unit_energy_in,opened=isopen)
   ! if (isopen) then
   !   write(*,'(a)') 'energy%read_txt() : Unit 14 is already open'
   !   close(unit_energy_in)
   ! else
   !   open(unit=unit_energy_in,file=dir//file_energy_in,action='read',iostat=iostatus,status='old')
   ! end if
   ! if(iostatus /= 0) then
   !   write(*,*) 'energy%read_txt(): file ', file_energy_in, ' not found'
   !   error stop
   !  end if

    OPEN (unit_band_in, file=dir//file_band_in, action='read', iostat=iostatus, status='old')
    OPEN (unit_band_out, file=dir//file_band_out, action='read', iostat=iostatus, status='old')
    OPEN (unit_mesh, file=dir//file_mesh, action='read', iostat=iostatus, status='old')
    OPEN (unit_hamiltonian_tb, file=dir//file_hamiltonian_tb, action='read', iostat=iostatus, status='old')


    ALLOCATE (x(0, 0), w(0))
    READ (unit_mesh, nml=mesh, iostat=iostatus)
    DEALLOCATE (x, w)
    TYPE = lower(TYPE)
    SELECT CASE (TYPE)
    CASE ('list')
        ALLOCATE (x(nx, 3))
        ALLOCATE (w(nx))
        x = 0.0_rp
        w = 1.0_rp/nx
        REWIND (unit_mesh)
        READ (unit_mesh, nml=mesh)
        CLOSE (unit_mesh)
        x_coord = lower(x_coord)
    END SELECT

    READ (unit_hamiltonian_tb, nml=hamiltonian_tb, iostat=iostatus)
    CLOSE (unit_hamiltonian_tb)

    SELECT CASE (ns)
    CASE (1)
        nsl = 1
    CASE (2)
        nsl = 2
    CASE (4)
        nsl = 1
    END SELECT
    na_band=0
    ALLOCATE(ia_band(0))
    READ(unit_band_in, nml=band, iostat=iostatus)
    DEALLOCATE(ia_band)
    ALLOCATE(ia_band(na_band))
    REWIND(unit_band_in)
    READ(unit_band_in, nml=band, iostat=iostatus)
    CLOSE(unit_band_in) 

    ALLOCATE (en_k(nh, nx, nsl))
    READ (unit_band_out, nml=band_out, iostat=iostatus)

    REWIND(unit_band_out)
    if(na_band>0) then
       ALLOCATE(iband2io(na_band))   
       iband2io=0 
       ALLOCATE(w_en_band_local(0,0,0,0,0))
       READ (unit_band_out, nml=band_out, iostat=iostatus) 
       no_max=maxval(iband2io)
       DEALLOCATE(w_en_band_local)
       REWIND(unit_band_out)
       ALLOCATE(w_en_band_local(na_band,no_max,nx,nh,nsl))
       READ (unit_band_out, nml=band_out, iostat=iostatus) 
    end if
    CLOSE (unit_band_out)
    CALL get_Fermi_scf(en_f)
    en_k = en_k - en_f
    CALL build_band_path(x, en_k, nh, nx, nsl)
     if(na_band>0) then
     CALL build_band_path_weight(x, en_k,w_en_band_local,iband2io,na_band,no_max, nh, nx, nsl)
     endif
END PROGRAM bands

SUBROUTINE build_band_path(x, en_k, nh, nx, nsl)
    USE precision_mod
    USE string_mod
    IMPLICIT NONE
    INTEGER, intent(in) :: nh, nx, nsl
    REAL(rp), intent(in) :: x(nx, 3), en_k(nh, nx, nsl)

    INTEGER :: ih, ix, isl
    INTEGER, PARAMETER :: unit_band = 10
    REAL(rp) :: sk
    CHARACTER(len=*), PARAMETER :: dir = 'band/'
    CHARACTER(len=*), PARAMETER :: file_band = 'band.dat'
    CHARACTER(len=80) :: FMT

    OPEN (unit=unit_band, file=dir//file_band, action='write')
    FMT = TRIM('('//int2str(nh + 1)//'F14.7'//')')

    DO isl = 1, nsl
        WRITE (unit_band, *) '@# k   band (eV)'
        sk = 0.0_rp
        WRITE (unit_band, FMT) sk, (en_k(ih, 1, isl), ih=1, nh)
        DO ix = 2, nx
            sk = sk + SQRT(SUM((x(ix, :) - x(ix - 1, :))**2))
            WRITE (unit_band, FMT) sk, (en_k(ih, ix, isl), ih=1, nh)
        END DO
    END DO
    WRITE (unit_band, *) '@# k   EF=0'
    WRITE (unit_band, '(2F8.3)') 0.0_rp, 0.0_rp
    WRITE (unit_band, '(2F8.3)') sk, 0.0_rp
    CLOSE (unit_band)
END SUBROUTINE build_band_path


SUBROUTINE build_band_path_weight(x, en_k,w_en_band_local,iband2io,na_band,no_max, nh, nx, nsl)
    USE precision_mod
    USE string_mod
    IMPLICIT NONE
    INTEGER, intent(in) :: na_band,no_max,nh, nx, nsl
    INTEGER, intent(in) :: iband2io(na_band)
    REAL(rp), intent(in) :: x(nx, 3), en_k(nh, nx, nsl)
    REAL(rp), intent(in) :: w_en_band_local(na_band,no_max,nx,nh,nsl)   
    INTEGER :: ih, ix, isl,ia_band,io
    INTEGER, PARAMETER :: unit_band = 10
    REAL(rp) :: sk
    CHARACTER(len=*), PARAMETER :: dir = 'band/'
    CHARACTER(len=*), PARAMETER :: file_band_weight = 'band_weight.dat'
    CHARACTER(len=80) :: FMT

    OPEN (unit=unit_band, file=dir//file_band_weight, action='write')

    DO isl = 1, nsl
        WRITE (unit_band, *) '@# k   band (eV)'
        DO ia_band=1,na_band
        WRITE (unit_band, *) '@# atom no', ia_band,'number of orbitals', iband2io(ia_band)
        sk = 0.0_rp
           FMT = TRIM('('//int2str(3+ iband2io(ia_band))//'F12.7'//')')         	
              do ih=1,nh
               WRITE (unit_band, FMT) sk, en_k(ih, 1, isl),sum(w_en_band_local(ia_band,:,1,ih,isl)),&
                                     (w_en_band_local(ia_band,io,1,ih,isl),io=1,iband2io(ia_band))
              end do
              DO ix = 2, nx
                sk = sk + SQRT(SUM((x(ix, :) - x(ix - 1, :))**2))
                do ih=1,nh
                  WRITE (unit_band, FMT) sk, en_k(ih, ix, isl),sum(w_en_band_local(ia_band,:,ix,ih,isl)),&
                                         (w_en_band_local(ia_band,io,ix,ih,isl),io=1,iband2io(ia_band))
                end do
            END DO
        END DO
    END DO
!    WRITE (unit_band, *) '@# k   EF=0'
!    WRITE (unit_band, '(2F8.3)') 0.0_rp, 0.0_rp
!    WRITE (unit_band, '(2F8.3)') sk, 0.0_rp
    CLOSE (unit_band)
END SUBROUTINE build_band_path_weight


SUBROUTINE get_Fermi_scf(en_f)
    USE precision_mod
    implicit none
    INTEGER, PARAMETER :: unit_energy_scf = 10
    INTEGER :: iostatus
    CHARACTER(len=*), PARAMETER :: file_energy_scf = 'out_energy.txt'
    REAL(rp) :: en_f
    NAMELIST /energy/en_f

    OPEN (unit_energy_scf, file=file_energy_scf, action='read', iostat=iostatus, status='old')
    READ (unit_energy_scf, nml=energy, iostat=iostatus)
   
    CLOSE (unit_energy_scf)
END SUBROUTINE get_Fermi_scf
