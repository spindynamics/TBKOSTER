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
!  dat2gnu.f90
!  TBKOSTER
    PROGRAM DAT2GNU
    implicit none
    real(8),dimension(3) :: k
    real(8),dimension(0:3) :: w
    integer :: ix,iy,iz,nx,ny,nz    
    
    write(*,*) 'enter nx ny and nz'
    read(*,*) nx,ny,nz
    
    open(unit=10,file='fermi.dat',action='read')
    open(unit=11,file='gnu-fermi.dat',action='write')
 !   read(10,*) 
    do ix=1, nx
    do iy=1, ny  
    do iz=1, nz   
    read(10,*) k(1),k(2),k(3),w(0)
    write(11,'(4F23.17)') k(1),k(2),k(3),w(0)
    end do
    end do   
    write(11,*)
    end do
    
    
    open(unit=12,file='vector_field.dat',action='read')
    open(unit=13,file='gnu-vector-field.dat',action='write')
 !   read(10,*) 
    do ix=1, nx
    do iy=1, ny  
    do iz=1, nz
    read(12,*) k(1),k(2),k(3),w(0), w(1),w(2),w(3)
    write(13,'(4F23.17)') k(1),k(2),w(1),w(2)
    end do
    end do
    write(13,*)
    end do
        
    
    
    END
    
    
    
