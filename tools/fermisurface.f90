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
!  fermisurface.f90
!  TBKOSTER
program fermisurface
    use, intrinsic :: iso_fortran_env, only: output_unit
    use precision_mod
    use string_mod
    implicit none
    integer, parameter :: unit_log_out = 11
    integer, parameter :: unit_hamiltonian_tb = 12
    integer, parameter :: unit_band_in = 13 
    integer, parameter :: unit_band_out = 14
    integer :: iostatus, nx, nh, ns, isl,nsl,na_band,no_max
    logical :: isopen
    character(len=*), parameter :: dir = 'band/'
    character(len=*), parameter :: file_log_out = 'out_log.txt'
    character(len=*), parameter :: file_band_in = 'in_band.txt'
    character(len=*), parameter :: file_band_out = 'out_band.txt'
    character(len=*), parameter :: file_hamiltonian_tb = 'out_hamiltonian_tb.txt'
    character(len=80) :: file_fs
    character(len=9) :: x_coord
    character(len=4) :: type
    character(len=80) :: line
    character(len=11) :: proj
    integer,dimension(3) :: gx,dx
    integer :: first_band,last_band
    integer, dimension(:), allocatable :: ia_band,iband2io
    real(rp), dimension(:, :, :), allocatable :: en_k
    real(rp), dimension(:), allocatable :: w
    real(rp), dimension(:, :), allocatable :: x,x_shifted
    real(rp) :: v_factor
    real(rp), dimension(3, 3) :: v,vrec
    real(rp), dimension(:,:,:,:,:), allocatable :: w_band_site
    real(rp), dimension(:,:,:,:), allocatable :: w_band_spin
    real(rp), dimension(:,:,:,:), allocatable :: w_band_orb
    real(rp) :: en_min, en_max, en_f
    namelist /lattice/v_factor,v,vrec
    namelist /mesh_out/type,gx,dx,nx
    namelist /hamiltonian_tb/nh, ns
    namelist /band/proj,na_band,ia_band 
    namelist /band_out/en_k,iband2io,w_band_site,w_band_spin,w_band_orb

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
    open (unit_log_out, file=dir//file_log_out, action='read', iostat=iostatus, status='old')
    open (unit_hamiltonian_tb, file=dir//file_hamiltonian_tb, action='read', iostat=iostatus, status='old')

   
    read (unit_log_out, nml=lattice, iostat=iostatus)
    read (unit_log_out, nml=mesh_out, iostat=iostatus)
    close (unit_log_out)
    type = lower(type)
    select case (type)
    case ('mp')
        if(dx(1)==0.and.dx(2)==0.and.dx(3)==0) then
            x=build_monkhorst_pack(gx)
            x_shifted=build_shifted_monkhorst_pack(gx)
        else
         write(*,*) 'dx must be set to zero'
         stop
        endif
        
    case('list')
        write(*,*) 'fermi surface only works with monkhorst pack'
        stop
    case('path')
        write(*,*) 'fermi surface only works with monkhorst pack'
        stop
    end select

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
    na_band=0
    allocate(ia_band(0))
    read(unit_band_in, nml=band, iostat=iostatus)
    deallocate(ia_band)
    allocate(ia_band(na_band))
    rewind(unit_band_in)
    read(unit_band_in, nml=band, iostat=iostatus)
    close(unit_band_in) 

    allocate (en_k(nh, nx, nsl))
    read (unit_band_out, nml=band_out, iostat=iostatus)

    rewind(unit_band_out)
    if(TRIM(proj)=='site'.AND.na_band>0) then
       allocate(iband2io(na_band))   
       iband2io=0 
       allocate(w_band_site(0,0,0,0,0))
       read (unit_band_out, nml=band_out, iostat=iostatus) 
       no_max=maxval(iband2io)
       deallocate(w_band_site)
       rewind(unit_band_out)
       allocate(w_band_site(nh,nx,na_band,no_max,nsl))
       read (unit_band_out, nml=band_out, iostat=iostatus) 
    elseif(TRIM(proj)=='spin'.AND.na_band>0) then
        allocate(w_band_spin(0,0,0,0))
        read (unit_band_out, nml=band_out, iostat=iostatus) 
        deallocate(w_band_spin)
        rewind(unit_band_out)
        allocate(w_band_spin(nh,nx,na_band,4))
        read (unit_band_out, nml=band_out, iostat=iostatus)
    elseif(TRIM(proj)=='orbit'.AND.na_band>0) then
        allocate(w_band_orb(0,0,0,0))
        read (unit_band_out, nml=band_out, iostat=iostatus) 
        deallocate(w_band_orb)
        rewind(unit_band_out)
        allocate(w_band_orb(nh,nx,na_band,4))
        read (unit_band_out, nml=band_out, iostat=iostatus)
    elseif(TRIM(proj)=='spin,orbit'.AND.na_band>0) then
        allocate(w_band_spin(0,0,0,0))
        allocate(w_band_orb(0,0,0,0))
        read (unit_band_out, nml=band_out, iostat=iostatus) 
        deallocate(w_band_spin)
        deallocate(w_band_orb)
        rewind(unit_band_out)
        allocate(w_band_spin(nh,nx,na_band,4))
        allocate(w_band_orb(nh,nx,na_band,4))
        read (unit_band_out, nml=band_out, iostat=iostatus)
    end if

    close (unit_band_out)
    call get_fermi_scf(en_f)
    en_k = en_k - en_f
    
    call cross_fermi(nh, nsl,nx, 0.0_rp, en_k,first_band, last_band )
    write(*,*) first_band,last_band
    do isl=1,nsl
        if(ns==1) then
           file_fs=dir//'fermi.bxsf'
        elseif(ns==2.and.isl==1)then
            file_fs=dir//'fermi.up.bxsf'
        elseif(ns==2.and.isl==2)then
            file_fs=dir//'fermi.dn.bxsf'
        elseif(ns==4)then
            file_fs=dir//'fermi.updn.bxsf'      
        endif     
       call xsf_fs (file_fs,0.0_rp, nh, first_band, last_band, isl, nx, &
                   gx(1),gx(2),gx(3),vrec(1,:), vrec(2,:), vrec(3,:), en_k )
    end do
   contains

   
subroutine get_fermi_scf(en_f)
    use precision_mod
    implicit none
    integer, parameter :: unit_energy_scf = 10
    integer :: iostatus
    character(len=*), parameter :: file_energy_scf = 'out_energy.txt'
    real(rp) :: en_f
    namelist /energy/en_f

    open (unit_energy_scf, file=file_energy_scf, action='read', iostat=iostatus, status='old')
    read (unit_energy_scf, nml=energy, iostat=iostatus)
   
    close (unit_energy_scf)
end subroutine get_fermi_scf

subroutine cross_fermi(nh, nsl,nkfs, efermi, e_fs,first_band, last_band )
    use precision_mod
 !----------------------------------------------------------------------- 
 !     Write  header file here
 !
 IMPLICIT NONE
 INTEGER, INTENT(IN) :: nh, nsl, nkfs
 INTEGER, INTENT(OUT) :: first_band, last_band
 REAL(rp), INTENT(IN) :: efermi
 REAL(rp), INTENT(IN) ::  e_fs(:,:,:)
 ! local variables
 INTEGER :: ih, ik
 REAL(rp), parameter ::  deltaE=1.0_rp
 REAL(rp) ::  emin, emax

    ! find bands that cross the Fermi surface at +- deltaE
    first_band = 0
    last_band  = 0
    DO ih=1,nh
       emin = MINVAL(e_fs(ih,1:nkfs,1:nsl))
       emax = MAXVAL(e_fs(ih,1:nkfs,1:nsl))
       IF ( emin-deltaE < efermi .AND. emax+deltaE > efermi ) THEN
          ! include all bands crossing the Fermi energy 
          ! within a tolerance +/- deltaE
          IF ( first_band == 0 ) first_band = ih
          last_band = ih
       END IF
    END DO
end subroutine cross_fermi

function build_monkhorst_pack(gx) result(x)
    use precision_mod
    integer,dimension(3),intent(in) :: gx
    real(rp),dimension(:,:), allocatable :: x
    integer :: ix,igx1,igx2,igx3
    real(rp) :: x1,x2,x3

   ! if (size(x,1)/=product(gx)) then
    !    write(*,*) 'size of x incompatible'
    !    stop
    !endif
    allocate(x(product(gx),3))

    ix = 1
    do igx1=1,gx(1)
      x1 = (igx1-1)/(gx(1)-1.0_rp)+dx(1)/(2*(gx(1)-1.0_rp))
      do igx2=1,gx(2)
        x2 =(igx2-1)/(gx(2)-1.0_rp)+dx(2)/(2*(gx(2)-1.0_rp))
        do igx3=1,gx(3)
          x3 = (igx3-1)/(gx(3)-1.0_rp)+dx(3)/(2*(gx(3)-1.0_rp))
          x(ix,:) = (/x1,x2,x3/)
          ix = ix+1
        end do
      end do
    end do
  end function build_monkhorst_pack

  function build_shifted_monkhorst_pack(gx) result(x)
    use precision_mod
    integer,dimension(3),intent(in) :: gx
    real(rp),dimension(:,:), allocatable :: x
    integer :: ix,igx1,igx2,igx3
    real(rp) :: x1,x2,x3

   ! if (size(x,1)/=product(gx)) then
    !    write(*,*) 'size of x incompatible'
    !    stop
    !endif
    allocate(x(product(gx),3))

    ix = 1
    do igx1=1,gx(1)
       x1 = (igx1-1)/(gx(1)-1.0_rp)+dx(1)/(2*(gx(1)-1.0_rp))
       if(x1>0.5_rp) x1=x1-1.0_rp
      do igx2=1,gx(2)
        x2 =(igx2-1)/(gx(2)-1.0_rp)+dx(2)/(2*(gx(2)-1.0_rp))
        if(x2>0.5_rp) x2=x2-1.0_rp
        do igx3=1,gx(3)
          x3 = (igx3-1)/(gx(3)-1.0_rp)+dx(3)/(2*(gx(3)-1.0_rp))
          if(x3>0.5_rp) x3=x3-1.0_rp
          x(ix,:) = (/x1,x2,x3/)
          ix = ix+1
        end do
      end do
    end do
  end function build_shifted_monkhorst_pack


  SUBROUTINE xsf_fs (filename,Ef, nh, nh_first,nh_last,isl, nkfs, na,nb,nc, ga, gb, gc, e_fs )
    use precision_mod
 !----------------------------------------------------------------------- 
 !     Write  header file here
 !
 IMPLICIT NONE
 CHARACTER(LEN=*), INTENT(IN) :: filename
 INTEGER, INTENT(IN) :: nh,nh_first,nh_last, isl, nkfs, na, nb, nc
 REAL(rp), INTENT(IN) :: Ef, ga(:), gb(:), gc(:), e_fs(:,:,:)
 REAL(rp) :: x0=0.0_rp, y0=0.0_rp, z0=0.0_rp
 INTEGER :: ih, ik, fsunit = 4
 !
 OPEN (unit = fsunit, file = filename, status='unknown', form='formatted')
 WRITE(fsunit, '(" BEGIN_INFO")')
 WRITE(fsunit, '("   #")') 
 WRITE(fsunit, '("   # this is a Band-XCRYSDEN-Structure-File")')
 WRITE(fsunit, '("   # aimed at Visualization of Fermi Surface")')
 WRITE(fsunit, '("   #")') 
 WRITE(fsunit, '("   # Case:   ",A)')     TRIM(filename)
 WRITE(fsunit, '("   #")') 
 WRITE(fsunit, '(" Fermi Energy:    ", f12.4)') Ef
 WRITE(fsunit, '(" END_INFO")') 
  
 WRITE(fsunit, '(" BEGIN_BLOCK_BANDGRID_3D")')
 WRITE(fsunit, '(" band_energies")')
 WRITE(fsunit, '(" BANDGRID_3D_BANDS")')
 WRITE(fsunit, '(I5)')  nh_last-nh_first+1
 WRITE(fsunit, '(3I5)') na, nb, nc
 WRITE(fsunit, '(3f10.6)') x0,  y0,  z0
 WRITE(fsunit, '(3f10.6)') ga(1), ga(2), ga(3)
 WRITE(fsunit, '(3f10.6)') gb(1), gb(2), gb(3)
 WRITE(fsunit, '(3f10.6)') gc(1), gc(2), gc(3)
 DO ih=nh_first,nh_last
    WRITE(fsunit, '("BAND:", i4)') ih
    WRITE(fsunit, '(6f10.4)') (e_fs(ih,ik,isl),ik=1,nkfs)
 END DO
 WRITE(fsunit, '(" END_BANDGRID_3D")')
 WRITE(fsunit, '(" END_BLOCK_BANDGRID_3D")')
 !
 CLOSE  (unit = fsunit, status='keep' )
 !
END SUBROUTINE xsf_fs


end program fermisurface