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
    INTEGER, PARAMETER :: unit_energy = 12
    INTEGER, PARAMETER :: unit_hamiltonian_tb = 13
    INTEGER :: iostatus, nx, nh, ns, nsl
    CHARACTER(len=*), PARAMETER :: dir = 'band/'
    CHARACTER(len=*), PARAMETER :: file_energy = 'out_energy.txt'
    CHARACTER(len=*), PARAMETER :: file_mesh = 'out_mesh.txt'
    CHARACTER(len=*), PARAMETER :: file_hamiltonian_tb = 'out_hamiltonian_tb.txt'
    CHARACTER(len=9) :: x_coord
    CHARACTER(len=4) :: TYPE
    REAL(rp), DIMENSION(:, :, :), ALLOCATABLE :: en_k
    REAL(rp), DIMENSION(:), ALLOCATABLE :: w
    REAL(rp), DIMENSION(:, :), ALLOCATABLE :: x
    REAL(rp) :: en_min, en_max, en_f
    NAMELIST /energy/ en_min, en_max, en_f, en_k
    NAMELIST /mesh/ TYPE, nx, x_coord, x, w
    NAMELIST /hamiltonian_tb/ nh, ns

    OPEN (unit_energy, file=dir//file_energy, action='read', iostat=iostatus, status='old')
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


    ALLOCATE (en_k(nh, nx, nsl))
    READ (unit_energy, nml=energy, iostat=iostatus)
    CALL get_Fermi_scf(en_f)
    en_k = en_k - en_f
    CLOSE (unit_energy)
    CALL build_band_path(x, en_k, nh, nx, nsl)
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
    !write(number,'(I3)') nh+1
    FMT = TRIM('('//int2str(nh + 1)//'F23.16'//')')

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

SUBROUTINE get_Fermi_scf(en_f)
    USE precision_mod
    implicit none
    INTEGER, PARAMETER :: unit_energy_scf = 10
    INTEGER :: iostatus
    CHARACTER(len=*), PARAMETER :: file_energy_scf = 'out_energy.txt'
    REAL(rp) :: en_f
    NAMELIST /energy/ en_f

    OPEN (unit_energy_scf, file=file_energy_scf, action='read', iostat=iostatus, status='old')
    READ (unit_energy_scf, nml=energy, iostat=iostatus)
   !write(*,*) 'EF=',en_f
    CLOSE (unit_energy_scf)
END SUBROUTINE get_Fermi_scf
