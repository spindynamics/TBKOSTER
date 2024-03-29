 !  TBKOSTER
    PROGRAM BUILD_KPOINTS
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
    
    
    
