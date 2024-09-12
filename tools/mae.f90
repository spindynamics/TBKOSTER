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
!  bands.f90
!  TBKOSTER
PROGRAM pmae
  USE, INTRINSIC :: iso_fortran_env, ONLY: output_unit
  USE precision_mod
  USE string_mod
  IMPLICIT NONE
  INTEGER, PARAMETER :: unit_mae_in=10
  INTEGER, PARAMETER :: unit_mae_out=11
  INTEGER, PARAMETER :: unit_log=12
  INTEGER :: iostatus,iatom
  CHARACTER(len=*),PARAMETER :: dir='mae/'
  CHARACTER(len=*),PARAMETER :: file_mae_in = 'in_mae.txt'
  CHARACTER(len=*),PARAMETER :: file_mae_out = 'out_mae'
  CHARACTER(len=*),PARAMETER :: file_out_log = 'out_log.txt'
  INTEGER :: ne,ns,na,na_mae,nangle,iangle
  INTEGER, DIMENSION(:), ALLOCATABLE :: ia,ia2ie
  CHARACTER(len=2),dimension(:),allocatable :: symbol
  REAL(rp),DIMENSION(2):: angle1,angle2
  REAL(rp),DIMENSION(3,3):: v
  REAL(rp),DIMENSION(:,:), ALLOCATABLE:: r
  REAL(rp) :: v_factor,mae_tot,mae_sum
  COMPLEX(rp),DIMENSION(:),ALLOCATABLE :: mae_s,mae_p,mae_px,mae_py,mae_pz
  COMPLEX(rp),DIMENSION(:),ALLOCATABLE :: mae_d,mae_dxy,mae_dyz,mae_dzx,mae_dx2y2,mae_dz2r2
  REAL(rp),DIMENSION(:),ALLOCATABLE :: mae_tot_i,mae_sum_i
  COMPLEX(rp),DIMENSION(:,:),ALLOCATABLE :: mae_s_i,mae_p_i,mae_px_i,mae_py_i,mae_pz_i
  COMPLEX(rp),DIMENSION(:,:),ALLOCATABLE :: mae_d_i,mae_dxy_i,mae_dyz_i,mae_dzx_i,&
                                            mae_dx2y2_i,mae_dz2r2_i
   NAMELIST /element/ne,Symbol
   NAMELIST /lattice/v_factor,v
   NAMELIST /atom/ns,na,ia2ie,r
   NAMELIST /mae/nangle,angle1,angle2,na_mae,ia
   NAMELIST /mae_out/mae_tot,mae_sum,mae_s,mae_p,mae_px,mae_py,mae_pz, &
       mae_d,mae_dxy,mae_dyz,mae_dzx,mae_dx2y2,mae_dz2r2

  OPEN(unit_log,file=dir // file_out_log,action='read',iostat=iostatus,status='old')
  allocate(Symbol(0))
  READ(unit_log,nml=element,iostat=iostatus)
  deallocate(Symbol)
  allocate(Symbol(ne))
  rewind(unit_log)
  READ(unit_log,nml=element,iostat=iostatus)
  rewind(unit_log)
  READ(unit_log,nml=lattice,iostat=iostatus)
  rewind(unit_log)
  allocate(r(0,0),ia2ie(0))
  READ(unit_log,nml=atom,iostat=iostatus)
  deallocate(ia2ie)
  allocate(ia2ie(na))
  rewind(unit_log)
  READ(unit_log,nml=atom,iostat=iostatus)
  deallocate(r)
  allocate(r(na,3))
   rewind(unit_log)
  READ(unit_log,nml=atom,iostat=iostatus)
  close(unit_log)

 OPEN(unit_mae_in,file=dir // file_mae_in,action='read',iostat=iostatus,status='old')
  allocate(ia(0))
  READ(unit_mae_in,nml=mae,iostat=iostatus)
  deallocate(ia)
  allocate(ia(na_mae))
  if(na_mae==na) then 
     do iatom=1,na
       ia(iatom)=iatom
     end do
   else
     rewind(unit_mae_in)
     READ(unit_mae_in,nml=mae,iostat=iostatus)
   endif
   close(unit_mae_in)

  ALLOCATE(mae_tot_i(nangle),mae_sum_i(nangle),mae_s_i(nangle,na_mae),mae_p_i(nangle,na_mae),&
          mae_px_i(nangle,na_mae),mae_py_i(nangle,na_mae),mae_pz_i(nangle,na_mae))         
  ALLOCATE(mae_d_i(nangle,na_mae),mae_dxy_i(nangle,na_mae),mae_dyz_i(nangle,na_mae),&
          mae_dzx_i(nangle,na_mae),mae_dx2y2_i(nangle,na_mae),mae_dz2r2_i(nangle,na_mae))

  ALLOCATE(mae_s(na_mae),mae_p(na_mae),mae_px(na_mae),mae_py(na_mae),mae_pz(na_mae))
  ALLOCATE(mae_d(na_mae),mae_dxy(na_mae),mae_dyz(na_mae),mae_dzx(na_mae),mae_dx2y2(na_mae),mae_dz2r2(na_mae))

  do iangle=1,nangle
  OPEN(unit_mae_out,file=dir // file_mae_out//'_'//int2str(iangle)//'.txt',action='read',iostat=iostatus,status='old')
  READ(unit_mae_out,nml=mae_out,iostat=iostatus)
  mae_tot_i(iangle )=mae_tot
  mae_sum_i(iangle)=mae_sum
  mae_s_i(iangle,:)=mae_s(:)
  mae_p_i(iangle,:)=mae_p(:)
  mae_d_i(iangle,:)=mae_d(:)
  mae_px_i(iangle,:)=mae_px(:)
  mae_py_i(iangle,:)=mae_py(:)
  mae_pz_i(iangle,:)=mae_pz(:)
  mae_dxy_i(iangle,:)=mae_dxy(:)
  mae_dyz_i(iangle,:)=mae_dyz(:)
  mae_dzx_i(iangle,:)=mae_dzx(:)
  mae_dx2y2_i(iangle,:)=mae_dx2y2(:)
  mae_dz2r2_i(iangle,:)=mae_dz2r2(:)
  close(unit_mae_out)
  end do

  CALL build_mae(na,ne,ia2ie,Symbol,v,r,nangle,angle1,angle2,na_mae,ia,&
             mae_tot_i,mae_sum_i,mae_s_i,mae_p_i,mae_d_i,&
             mae_px_i,mae_py_i,mae_pz_i,mae_dxy_i,mae_dyz_i,mae_dzx_i,mae_dx2y2_i,mae_dz2r2_i)
END PROGRAM pmae

SUBROUTINE build_mae(na,ne,ia2ie,Symbol,v,r,nangle,angle1,angle2,na_mae,ia,&
             mae_tot_i,mae_sum_i,mae_s_i,mae_p_i,mae_d_i,&
             mae_px_i,mae_py_i,mae_pz_i,mae_dxy_i,mae_dyz_i,mae_dzx_i,mae_dx2y2_i,mae_dz2r2_i)
  USE precision_mod
  USE string_mod
  USE math_mod, ONLY: i_unit
  IMPLICIT NONE
  INTEGER :: na,ne,nangle,na_mae
  INTEGER :: ia_mae,iangle,iatom
  INTEGER, PARAMETER :: unit_mae_theta=10
  INTEGER, PARAMETER :: unit_mae_xyz=11
  INTEGER, PARAMETER :: unit_mae_atom=12
  INTEGER :: ia(na_mae),ia2ie(na)
  CHARACTER(len=2),dimension(ne):: Symbol
  REAL(rp), dimension(2):: angle1,angle2
  REAL(rp), dimension(na,3):: r 
  REAL(rp), dimension(3,3):: v
  REAL(rp) ::   angle
  REAL(rp) ::    mae_tot_i(nangle),mae_sum_i(nangle)
  COMPLEX(rp) :: mae_s_i(nangle,na_mae),mae_p_i(nangle,na_mae),&
                 mae_px_i(nangle,na_mae),mae_py_i(nangle,na_mae),mae_pz_i(nangle,na_mae)
  COMPLEX(rp) :: mae_d_i(nangle,na_mae),mae_dxy_i(nangle,na_mae),mae_dyz_i(nangle,na_mae),&
          mae_dzx_i(nangle,na_mae),mae_dx2y2_i(nangle,na_mae),mae_dz2r2_i(nangle,na_mae)
  CHARACTER(len=*),PARAMETER :: dir='mae/'
  CHARACTER(len=*),PARAMETER :: file_mae_theta = 'mae_angle.dat'
  CHARACTER(len=*),PARAMETER :: file_mae_xyz = 'mae.xyz'
  CHARACTER(len=*),PARAMETER :: file_mae_atom= 'mae_atom.dat'
  CHARACTER(len=80) :: FMT

  OPEN(unit=unit_mae_atom,file=dir //file_mae_atom,action='write')
   WRITE(unit_mae_atom,*) '@# i MAE(meV)'

         Write(unit_mae_atom,*) '@# site mae s'
         DO ia_mae=1,na_mae
         WRITE(unit_mae_atom,'(2(a,1X))') int2str(ia(ia_mae)), &
          real2str(1000*REAL(mae_s_i(nangle,ia_mae)-mae_s_i(1,ia_mae)))
         END DO

         Write(unit_mae_atom,*) '@# site mae p'
         DO ia_mae=1,na_mae
         WRITE(unit_mae_atom,'(5(a,1X))') int2str(ia(ia_mae)), &
          real2str(1000*REAL(mae_p_i(nangle,ia_mae)-mae_p_i(1,ia_mae))),&
          real2str(1000*REAL(mae_px_i(nangle,ia_mae)-mae_px_i(1,ia_mae))),&
          real2str(1000*REAL(mae_py_i(nangle,ia_mae)-mae_py_i(1,ia_mae))),&
          real2str(1000*REAL(mae_pz_i(nangle,ia_mae)-mae_pz_i(1,ia_mae)))
         END DO
         Write(unit_mae_atom,*) '@# site mae d'
         DO ia_mae=1,na_mae
         WRITE(unit_mae_atom,'(7(a,1X))') int2str(ia(ia_mae)), &
          real2str(1000*REAL(mae_d_i(nangle,ia_mae)-mae_d_i(1,ia_mae))),&
          real2str(1000*REAL(mae_dxy_i(nangle,ia_mae)-mae_dxy_i(1,ia_mae))),&
          real2str(1000*REAL(mae_dyz_i(nangle,ia_mae)-mae_dyz_i(1,ia_mae))),&
          real2str(1000*REAL(mae_dzx_i(nangle,ia_mae)-mae_dzx_i(1,ia_mae))),&
          real2str(1000*REAL(mae_dx2y2_i(nangle,ia_mae)-mae_dx2y2_i(1,ia_mae))),&
          real2str(1000*REAL(mae_dz2r2_i(nangle,ia_mae)-mae_dz2r2_i(1,ia_mae)))
         END DO

   close(unit_mae_atom)

 OPEN(unit=unit_mae_theta,file=dir //file_mae_theta,action='write')
 WRITE(unit_mae_theta,*) '@# angle  Energy(eV)'
       
       DO iangle=1,nangle
       angle=(iangle-1)*sqrt((angle2(1)-angle1(1))**2+(angle2(2)-angle1(2))**2)/(nangle-1)
        WRITE(unit_mae_theta,'(3(a,1X))') real2str(angle),real2str(1000*(mae_tot_i(iangle)-mae_tot_i(1)))
       END DO

  if (na_mae==na) then
      OPEN(unit=unit_mae_xyz,file=dir //file_mae_xyz,action='write')
       ! Number of atoms
    write(unit_mae_xyz,'(a)') int2str(na)
    ! Serie of key/value pairs
    write(unit_mae_xyz,'(a)') 'Lattice="' &
     // real2str(v(1,1)) // ' ' &
     // real2str(v(1,2)) // ' ' &
     // real2str(v(1,3)) // ' ' &
     // real2str(v(2,1)) // ' ' &
     // real2str(v(2,2)) // ' ' &
     // real2str(v(2,3)) // ' ' &
     // real2str(v(3,1)) // ' ' &
     // real2str(v(3,2)) // ' ' &
     // real2str(v(3,3)) // '" ' &
     // 'Properties=species:S:1:pos:R:3:kinetic_energy:R:1'
    ! Atoms
    do iatom=1,na
      write(unit_mae_xyz,'(a)') symbol(ia2ie(iatom)) // ' ' &
       // real2str(r(iatom,1)) // ' ' &
       // real2str(r(iatom,2)) // ' ' &
       // real2str(r(iatom,3)) // ' ' &
       // real2str(1000*real(mae_s_i(nangle,iatom)+mae_p_i(nangle,iatom)+mae_d_i(nangle,iatom)- &
                             (mae_s_i(1,iatom)+mae_p_i(1,iatom)+mae_d_i(1,iatom)))) // ' ' 
    enddo
    close(unit_mae_xyz)
  endif


END SUBROUTINE build_mae
