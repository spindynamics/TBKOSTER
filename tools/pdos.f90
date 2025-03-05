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
PROGRAM pdos
  USE, INTRINSIC :: iso_fortran_env, ONLY: output_unit
  USE precision_mod
  USE string_mod
  IMPLICIT NONE
  INTEGER, PARAMETER :: unit_dos=10
  INTEGER, PARAMETER :: unit_hamiltonian_tb=11
  INTEGER, PARAMETER :: unit_energy_scf=12
  INTEGER :: iostatus,nh,ns,nsl,nen,na_dos
  CHARACTER(len=*),PARAMETER :: dir='dos/'
  CHARACTER(len=*),PARAMETER :: file_dos = 'out_dos.txt'
  CHARACTER(len=*),PARAMETER :: file_hamiltonian_tb = 'out_hamiltonian_tb.txt'
  CHARACTER(len=*),PARAMETER :: file_energy_scf = 'out_energy.txt'
  INTEGER, DIMENSION(:), ALLOCATABLE :: ia
  REAL(rp),DIMENSION(:,:),ALLOCATABLE :: dos_tot
  COMPLEX(rp),DIMENSION(:,:,:),ALLOCATABLE :: dos_s,dos_p,dos_px,dos_py,dos_pz
  COMPLEX(rp),DIMENSION(:,:,:),ALLOCATABLE :: dos_d,dos_dxy,dos_dyz,dos_dzx,dos_dx2y2,dos_dz2r2
  REAL(rp) :: en_min,en_max,en_f
  NAMELIST /dos/nen,en_min,en_max,ia,na_dos,dos_tot,dos_s,dos_p,dos_px,dos_py,dos_pz, &
       dos_d,dos_dxy,dos_dyz,dos_dzx,dos_dx2y2,dos_dz2r2
  NAMELIST /hamiltonian_tb/nh,ns
  NAMELIST /energy/en_f

  OPEN(unit_dos,file=dir // file_dos,action='read',iostat=iostatus,status='old')
  OPEN(unit_hamiltonian_tb,file=dir // file_hamiltonian_tb,action='read',iostat=iostatus,status='old')

  READ(unit_hamiltonian_tb,nml=hamiltonian_tb,iostat=iostatus)
  CLOSE(unit_hamiltonian_tb)

  SELECT CASE(ns)
  CASE(1)
     nsl=1
  CASE(2)
     nsl=2
  CASE(4)
     nsl=1
  END SELECT

  ALLOCATE(ia(0))
  ALLOCATE(dos_tot(0,0),dos_s(0,0,0),dos_p(0,0,0),dos_px(0,0,0),dos_py(0,0,0),dos_pz(0,0,0))
  ALLOCATE(dos_d(0,0,0),dos_dxy(0,0,0),dos_dyz(0,0,0),dos_dzx(0,0,0),dos_dx2y2(0,0,0),dos_dz2r2(0,0,0))
  READ(unit_dos,nml=dos,iostat=iostatus)
  DEALLOCATE(ia)
  DEALLOCATE(dos_tot,dos_s,dos_p,dos_px,dos_py,dos_pz)
  DEALLOCATE(dos_d,dos_dxy,dos_dyz,dos_dzx,dos_dx2y2,dos_dz2r2)

  ALLOCATE(ia(na_dos))
  ALLOCATE(dos_tot(nen,nsl))
  ALLOCATE(dos_s(na_dos,nen,ns),dos_p(na_dos,nen,ns),dos_px(na_dos,nen,ns),dos_py(na_dos,nen,ns),dos_pz(na_dos,nen,ns))
  ALLOCATE(dos_d(na_dos,nen,ns),dos_dxy(na_dos,nen,ns),dos_dyz(na_dos,nen,ns),dos_dzx(na_dos,nen,ns),&
       dos_dx2y2(na_dos,nen,ns),dos_dz2r2(na_dos,nen,ns))
  REWIND(unit_dos)
  READ(unit_dos,nml=dos,iostat=iostatus)
  OPEN(unit_energy_scf,file=file_energy_scf,action='read',iostat=iostatus,status='old')
  READ(unit_energy_scf,nml=energy,iostat=iostatus)
  CLOSE(unit_dos)
  CLOSE(unit_energy_scf)

  CALL build_pdos(en_f,nen,na_dos,ns,nsl,ia,en_min,en_max,dos_tot,dos_s,dos_p,dos_px,dos_py,dos_pz,&
       dos_d,dos_dxy,dos_dyz,dos_dzx,dos_dx2y2,dos_dz2r2)
END PROGRAM pdos

SUBROUTINE build_pdos(en_f,nen,na_dos,ns,nsl,ia,en_min,en_max,dos_tot,dos_s,dos_p,dos_px,dos_py,dos_pz,&
     dos_d,dos_dxy,dos_dyz,dos_dzx,dos_dx2y2,dos_dz2r2)
  USE precision_mod
  USE string_mod
  USE math_mod, ONLY: i_unit
  IMPLICIT NONE
  INTEGER :: nen,na_dos,ns,nsl
  INTEGER :: ien,ia_dos,isl,ispin_up,ispin_dn,g_s
  INTEGER, PARAMETER :: unit_dos_tot=10
  INTEGER, PARAMETER :: unit_pdos_s=11
  INTEGER, PARAMETER :: unit_pdos_p=12
  INTEGER, PARAMETER :: unit_pdos_d=13
  INTEGER, PARAMETER :: unit_pdos_spd=14
  INTEGER :: ia(na_dos)
  REAL(rp) :: en,den
  REAL(rp) :: en_min,en_max,en_f
  REAL(rp) :: dos_tot(nen,nsl)
  COMPLEX(rp) :: dos_s(na_dos,nen,ns),dos_p(na_dos,nen,ns),dos_px(na_dos,nen,ns),dos_py(na_dos,nen,ns),dos_pz(na_dos,nen,ns)
  COMPLEX(rp) :: dos_d(na_dos,nen,ns),dos_dxy(na_dos,nen,ns),dos_dyz(na_dos,nen,ns),dos_dzx(na_dos,nen,ns),&
       dos_dx2y2(na_dos,nen,ns),dos_dz2r2(na_dos,nen,ns)
  CHARACTER(len=*),PARAMETER :: dir='dos/'
  CHARACTER(len=*),PARAMETER :: file_dos_tot = 'dos-tot.dat'
  CHARACTER(len=*),PARAMETER :: file_pdos_s = 'pdos-s.dat'
  CHARACTER(len=*),PARAMETER :: file_pdos_p = 'pdos-p.dat'
  CHARACTER(len=*),PARAMETER :: file_pdos_d = 'pdos-d.dat'
  CHARACTER(len=*),PARAMETER :: file_pdos_spd = 'pdos-spd.dat'
  CHARACTER(len=80) :: FMT
  REAL(rp) ::   dos_tot_x,dos_tot_y,dos_tot_z,dos_tot_n
  REAL(rp) ::   dos_s_x,dos_s_y,dos_s_z,dos_s_n
  REAL(rp) ::   dos_p_x,dos_p_y,dos_p_z,dos_p_n
  REAL(rp) ::   dos_px_x,dos_px_y,dos_px_z,dos_px_n
  REAL(rp) ::   dos_py_x,dos_py_y,dos_py_z,dos_py_n
  REAL(rp) ::   dos_pz_x,dos_pz_y,dos_pz_z,dos_pz_n
  REAL(rp) ::   dos_d_x,dos_d_y,dos_d_z,dos_d_n
  REAL(rp) ::   dos_dxy_x,dos_dxy_y,dos_dxy_z,dos_dxy_n
  REAL(rp) ::   dos_dyz_x,dos_dyz_y,dos_dyz_z,dos_dyz_n
  REAL(rp) ::   dos_dzx_x,dos_dzx_y,dos_dzx_z,dos_dzx_n
  REAL(rp) ::   dos_dx2y2_x,dos_dx2y2_y,dos_dx2y2_z,dos_dx2y2_n
  REAL(rp) ::   dos_dz2r2_x,dos_dz2r2_y,dos_dz2r2_z,dos_dz2r2_n
  REAL(rp) ::   dos_spd_x,dos_spd_y,dos_spd_z,dos_spd_n

  den=(en_max-en_min)/nen
  IF(ns==1.OR.ns==2) THEN
     IF(ns==1) THEN
        ispin_up=1
        ispin_dn=1
        g_s=2
     ELSEIF(ns==2) THEN
        ispin_up=1
        ispin_dn=2
        g_s=1
     end if
  end if

  OPEN(unit=unit_dos_tot,file=dir //file_dos_tot,action='write')
  OPEN(unit=unit_pdos_s,file=dir //file_pdos_s,action='write')
  OPEN(unit=unit_pdos_p,file=dir //file_pdos_p,action='write')
  OPEN(unit=unit_pdos_d,file=dir //file_pdos_d,action='write')
  OPEN(unit=unit_pdos_spd,file=dir //file_pdos_spd,action='write')

  IF(ns==1.OR.ns==4) THEN
     WRITE(unit_dos_tot,*) '@# E   TOTAL DOS'
     DO ien=1,nen
        en=en_min+(ien-1)*den
        WRITE(unit_dos_tot,'(2(a,1X))') real2str(en),real2str(dos_tot(ien,1))
     END DO
  ELSEIF(ns==2) THEN
     WRITE(unit_dos_tot,*) '@# E   TOTAL DOS up and down'
     DO ien=1,nen
        en=en_min+(ien-1)*den
        WRITE(unit_dos_tot,'(3(a,1X))') real2str(en),real2str(dos_tot(ien,1)),real2str(dos_tot(ien,2))
     END DO
  end if
  WRITE(unit_dos_tot,*) '@# EF'
  WRITE(unit_dos_tot,'(2F8.3)') en_f,0.0_rp
  WRITE(unit_dos_tot,'(2F8.3)') en_f,10.0_rp

  IF(na_dos>0) THEN
     IF(ns==1.OR.ns==2) THEN
        DO ia_dos=1,na_dos
           WRITE(unit_pdos_s,*) '@#dos of site',ia(ia_dos),'s(up) s(dn)'
           WRITE(unit_pdos_p,*) '@#dos of site',ia(ia_dos),'p(up) p(dn) px(up) px(dn) py(up) py(dn) pz(up) pz(dn)'
           WRITE(unit_pdos_d,*) '@#dos of site',ia(ia_dos),'d(up) d(dn) dxy(up) dxy(dn) dyz(up) dyz(dn) dxz(up) &
                dxz(dn) dy2y2(up) dx2y2(dn) dz2(up) dz2(dn)  '
           WRITE(unit_pdos_spd,*) '@#dos of site',ia(ia_dos),'s(up) s(dn) p(up) p(dn) d(up) d(dn)'

           DO ien=1,nen
              en=en_min+(ien-1)*den
              WRITE(unit_pdos_s,'(3(a,1X))') real2str(en),real2str(REAL(dos_s(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_s(ia_dos,ien,ispin_dn)/g_s))
              WRITE(unit_pdos_p,'(9(a,1X))') real2str(en),real2str(REAL(dos_p(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_p(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_px(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_px(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_py(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_py(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_pz(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_pz(ia_dos,ien,ispin_dn)/g_s))
              WRITE(unit_pdos_d,'(13(a,1X))') real2str(en),real2str(REAL(dos_d(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_d(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_dxy(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_dxy(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_dyz(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_dyz(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_dzx(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_dzx(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_dx2y2(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_dx2y2(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_dz2r2(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_dz2r2(ia_dos,ien,ispin_dn)/g_s))

              WRITE(unit_pdos_spd,'(7(a,1X))') real2str(en),real2str(REAL(dos_s(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_s(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_p(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_p(ia_dos,ien,ispin_dn)/g_s)),&
                   real2str(REAL(dos_d(ia_dos,ien,ispin_up)/g_s)),&
                   real2str(REAL(dos_d(ia_dos,ien,ispin_dn)/g_s))
           END DO
        END DO
     ELSEIF(ns==4) THEN
        DO ia_dos=1,na_dos
           WRITE(unit_pdos_s,*) '@#dos of site',ia(ia_dos),'s(n) s(x) s(y) s(z)'
           WRITE(unit_pdos_p,*) '@#dos of site',ia(ia_dos),'p(n) p(x) p(y) p(z) px(n) &
                px(x) px(y) px(z) py(n) py(x) py(y) py(z) pz(n) pz(x) pz(y) pz(z)'
           WRITE(unit_pdos_d,*) '@#dos of site',ia(ia_dos),'d(n) d(x) d(y) d(z) dxy(n) dxy(x) &
                dxy(y) dxy(z) dyz(n) dyz(x) dyz(y) dyz(z) dxz(n) dxz(x) dxz(y) &
                dxz(z)  dx2y2(n) dx2y2(x) dx2y2(y) dx2y2(z) dz2(n) dz2(x) dz2(y) dz2(z)'
           WRITE(unit_pdos_spd,*) '@#dos of site',ia(ia_dos),'s(n) s(x) s(y) s(z) p(n) p(x) p(y) p(z) d(n) d(x) d(y) d(z)'

           DO ien=1,nen
              en=en_min+(ien-1)*den

              dos_s_n=dos_s(ia_dos,ien,1)+dos_s(ia_dos,ien,2)
              dos_s_x=dos_s(ia_dos,ien,3)+dos_s(ia_dos,ien,4)
              dos_s_y=i_unit*(dos_s(ia_dos,ien,3)-dos_s(ia_dos,ien,4))
              dos_s_z=dos_s(ia_dos,ien,1)-dos_s(ia_dos,ien,2)

              dos_p_n=dos_p(ia_dos,ien,1)+dos_p(ia_dos,ien,2)
              dos_p_x=dos_p(ia_dos,ien,3)+dos_p(ia_dos,ien,4)
              dos_p_y=i_unit*(dos_p(ia_dos,ien,3)-dos_p(ia_dos,ien,4))
              dos_p_z=dos_p(ia_dos,ien,1)-dos_p(ia_dos,ien,2)

              dos_px_n=dos_px(ia_dos,ien,1)+dos_px(ia_dos,ien,2)
              dos_px_x=dos_px(ia_dos,ien,3)+dos_px(ia_dos,ien,4)
              dos_px_y=i_unit*(dos_px(ia_dos,ien,3)-dos_px(ia_dos,ien,4))
              dos_px_z=dos_px(ia_dos,ien,1)-dos_px(ia_dos,ien,2)

              dos_py_n=dos_py(ia_dos,ien,1)+dos_py(ia_dos,ien,2)
              dos_py_x=dos_py(ia_dos,ien,3)+dos_py(ia_dos,ien,4)
              dos_py_y=i_unit*(dos_py(ia_dos,ien,3)-dos_py(ia_dos,ien,4))
              dos_py_z=dos_py(ia_dos,ien,1)-dos_py(ia_dos,ien,2)

              dos_pz_n=dos_pz(ia_dos,ien,1)+dos_pz(ia_dos,ien,2)
              dos_pz_x=dos_pz(ia_dos,ien,3)+dos_pz(ia_dos,ien,4)
              dos_pz_y=i_unit*(dos_pz(ia_dos,ien,3)-dos_pz(ia_dos,ien,4))
              dos_pz_z=dos_pz(ia_dos,ien,1)-dos_pz(ia_dos,ien,2)

              dos_d_n=dos_d(ia_dos,ien,1)+dos_d(ia_dos,ien,2)
              dos_d_x=dos_d(ia_dos,ien,3)+dos_d(ia_dos,ien,4)
              dos_d_y=i_unit*(dos_d(ia_dos,ien,3)-dos_d(ia_dos,ien,4))
              dos_d_z=dos_d(ia_dos,ien,1)-dos_d(ia_dos,ien,2)


              dos_dxy_n=dos_dxy(ia_dos,ien,1)+dos_dxy(ia_dos,ien,2)
              dos_dxy_x=dos_dxy(ia_dos,ien,3)+dos_dxy(ia_dos,ien,4)
              dos_dxy_y=i_unit*(dos_dxy(ia_dos,ien,3)-dos_dxy(ia_dos,ien,4))
              dos_dxy_z=dos_dxy(ia_dos,ien,1)-dos_dxy(ia_dos,ien,2)

              dos_dyz_n=dos_dyz(ia_dos,ien,1)+dos_dyz(ia_dos,ien,2)
              dos_dyz_x=dos_dyz(ia_dos,ien,3)+dos_dyz(ia_dos,ien,4)
              dos_dyz_y=i_unit*(dos_dyz(ia_dos,ien,3)-dos_dyz(ia_dos,ien,4))
              dos_dyz_z=dos_dyz(ia_dos,ien,1)-dos_dyz(ia_dos,ien,2)

              dos_dzx_n=dos_dzx(ia_dos,ien,1)+dos_dzx(ia_dos,ien,2)
              dos_dzx_x=dos_dzx(ia_dos,ien,3)+dos_dzx(ia_dos,ien,4)
              dos_dzx_y=i_unit*(dos_dzx(ia_dos,ien,3)-dos_dzx(ia_dos,ien,4))
              dos_dzx_z=dos_dzx(ia_dos,ien,1)-dos_dzx(ia_dos,ien,2)

              dos_dx2y2_n=dos_dx2y2(ia_dos,ien,1)+dos_dx2y2(ia_dos,ien,2)
              dos_dx2y2_x=dos_dx2y2(ia_dos,ien,3)+dos_dx2y2(ia_dos,ien,4)
              dos_dx2y2_y=i_unit*(dos_dx2y2(ia_dos,ien,3)-dos_dx2y2(ia_dos,ien,4))
              dos_dx2y2_z=dos_dx2y2(ia_dos,ien,1)-dos_dx2y2(ia_dos,ien,2)

              dos_dz2r2_n=dos_dz2r2(ia_dos,ien,1)+dos_dz2r2(ia_dos,ien,2)
              dos_dz2r2_x=dos_dz2r2(ia_dos,ien,3)+dos_dz2r2(ia_dos,ien,4)
              dos_dz2r2_y=i_unit*(dos_dz2r2(ia_dos,ien,3)-dos_dz2r2(ia_dos,ien,4))
              dos_dz2r2_z=dos_dz2r2(ia_dos,ien,1)-dos_dz2r2(ia_dos,ien,2)

              dos_spd_n=dos_s_n+dos_p_n+dos_d_n
              dos_spd_x=dos_s_x+dos_p_x+dos_d_x
              dos_spd_y=dos_s_y+dos_p_y+dos_d_y
              dos_spd_z=dos_s_z+dos_p_z+dos_d_z

              WRITE(unit_pdos_s,'(5(a,1X))') real2str(en),real2str(dos_s_n),&
                   real2str(dos_s_x),real2str(dos_s_y),real2str(dos_s_z)
              WRITE(unit_pdos_p,'(17(a,1X))') real2str(en),real2str(dos_p_n), &
                   real2str(dos_p_x),real2str(dos_p_y),real2str(dos_p_z),&
                   real2str(dos_px_n),real2str(dos_px_x),&
                   real2str(dos_px_y),real2str(dos_px_z), &
                   real2str(dos_py_n),real2str(dos_py_x), &
                   real2str(dos_py_y),real2str(dos_py_z), &
                   real2str(dos_pz_n),real2str(dos_pz_x), &
                   real2str(dos_pz_y),real2str(dos_pz_z)
              WRITE(unit_pdos_d,'(25(a,1x))') real2str(en),real2str(dos_d_n),real2str(dos_d_x), &
                   real2str(dos_d_y),real2str(dos_d_z),&
                   real2str(dos_dxy_n),real2str(dos_dxy_x), &
                   real2str(dos_dxy_y),real2str(dos_dxy_z), &
                   real2str(dos_dyz_n),real2str(dos_dyz_x), &
                   real2str(dos_dyz_y),real2str(dos_dyz_z), &
                   real2str(dos_dzx_n),real2str(dos_dzx_x), &
                   real2str(dos_dzx_y),real2str(dos_dzx_z), &
                   real2str(dos_dx2y2_n),real2str(dos_dx2y2_x), &
                   real2str(dos_dx2y2_y),real2str(dos_dx2y2_z), &
                   real2str(dos_dz2r2_n),real2str(dos_dz2r2_x), &
                   real2str(dos_dz2r2_y),real2str(dos_dz2r2_z)
              WRITE(unit_pdos_spd,'(13(a,1x))') real2str(en),real2str(dos_s_n),real2str(dos_s_x), &
                   real2str(dos_s_y),real2str(dos_s_z), &
                   real2str(dos_p_n),real2str(dos_p_x), &
                   real2str(dos_p_y),real2str(dos_p_z),&
                   real2str(dos_d_n),real2str(dos_d_x),  &
                   real2str(dos_d_y), real2str(dos_d_z)
           END DO
        END DO
     end if
     WRITE(unit_pdos_s,*) '@# EF'
     WRITE(unit_pdos_s,'(2F8.3)') en_f,0.0_rp
     WRITE(unit_pdos_s,'(2F8.3)') en_f,10.0_rp
     WRITE(unit_pdos_p,*) '@# EF'
     WRITE(unit_pdos_p,'(2F8.3)') en_f,0.0_rp
     WRITE(unit_pdos_p,'(2F8.3)') en_f,10.0_rp
     WRITE(unit_pdos_d,*) '@# EF'
     WRITE(unit_pdos_d,'(2F8.3)') en_f,0.0_rp
     WRITE(unit_pdos_d,'(2F8.3)') en_f,10.0_rp
     WRITE(unit_pdos_spd,*) '@# EF'
     WRITE(unit_pdos_spd,'(2F8.3)') en_f,0.0_rp
     WRITE(unit_pdos_spd,'(2F8.3)') en_f,10.0_rp


  end if
END SUBROUTINE build_pdos
