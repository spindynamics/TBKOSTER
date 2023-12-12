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
program bands
    use, intrinsic :: iso_fortran_env, only: output_unit
    use precision_mod
    use string_mod
    implicit none
    integer, parameter :: unit_mesh = 11
    integer, parameter :: unit_hamiltonian_tb = 12
    integer, parameter :: unit_band_in = 13 
    integer, parameter :: unit_band_out = 14
    integer :: iostatus, nx, nh, ns, nsl, no_max, na_band
    logical :: isopen
    character(len=*), parameter :: dir = 'band/'
    character(len=*), parameter :: file_band_in = 'in_band.txt' 
    character(len=*), parameter :: file_band_out = 'out_band.txt'
    character(len=*), parameter :: file_mesh = 'out_mesh.txt'
    character(len=*), parameter :: file_hamiltonian_tb = 'out_hamiltonian_tb.txt'
    character(len=9) :: x_coord
    character(len=4) :: type
    character(len=80) :: line
    integer, dimension(:), allocatable :: ia_band,iband2io
    real(rp), dimension(:, :, :), allocatable :: en_k
    real(rp), dimension(:), allocatable :: w
    real(rp), dimension(:, :), allocatable :: x
    real(rp), dimension(:,:,:,:,:), allocatable :: w_en_band_local
    real(rp) :: en_min, en_max, en_f
    namelist /mesh/ type, nx, x_coord, x, w
    namelist /hamiltonian_tb/ nh, ns
    namelist /band/ na_band,ia_band 
    namelist /band_out/en_k,iband2io,w_en_band_local

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
    open (unit_mesh, file=dir//file_mesh, action='read', iostat=iostatus, status='old')
    open (unit_hamiltonian_tb, file=dir//file_hamiltonian_tb, action='read', iostat=iostatus, status='old')


    allocate (x(0, 0), w(0))
    read (unit_mesh, nml=mesh, iostat=iostatus)
    deallocate (x, w)
    type = lower(type)
    select case (type)
    case ('list')
        allocate (x(nx, 3))
        allocate (w(nx))
        x = 0.0_rp
        w = 1.0_rp/nx
        rewind (unit_mesh)
        read (unit_mesh, nml=mesh)
        close (unit_mesh)
        x_coord = lower(x_coord)
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
    if(na_band>0) then
       allocate(iband2io(na_band))   
       iband2io=0 
       allocate(w_en_band_local(0,0,0,0,0))
       read (unit_band_out, nml=band_out, iostat=iostatus) 
       no_max=maxval(iband2io)
       deallocate(w_en_band_local)
       rewind(unit_band_out)
       allocate(w_en_band_local(na_band,no_max,nx,nh,nsl))
       read (unit_band_out, nml=band_out, iostat=iostatus) 
    end if
    close (unit_band_out)
    call get_fermi_scf(en_f)
    en_k = en_k - en_f
    call build_band_path(x, en_k, nh, nx, nsl)
     if(na_band>0) then
     call build_band_path_weight(x, en_k,w_en_band_local,iband2io,na_band,no_max, nh, nx, nsl)
     endif
end program bands

subroutine build_band_path(x, en_k, nh, nx, nsl)
    use precision_mod
    use string_mod
    implicit none
    integer, intent(in) :: nh, nx, nsl
    real(rp), intent(in) :: x(nx, 3), en_k(nh, nx, nsl)

    integer :: ih, ix, isl
    integer, parameter :: unit_band = 10
    real(rp) :: sk
    character(len=*), parameter :: dir = 'band/'
    character(len=*), parameter :: file_band = 'band.dat'
    character(len=80) :: fmt

    open (unit=unit_band, file=dir//file_band, action='write')
    fmt = trim('('//int2str(nh + 1)//'f14.7'//')')

    do isl = 1, nsl
        write (unit_band, *) '@# k   band (ev)'
        sk = 0.0_rp
        write (unit_band, fmt) sk, (en_k(ih, 1, isl), ih=1, nh)
        do ix = 2, nx
            sk = sk + sqrt(sum((x(ix, :) - x(ix - 1, :))**2))
            write (unit_band, fmt) sk, (en_k(ih, ix, isl), ih=1, nh)
        end do
    end do
    write (unit_band, *) '@# k   ef=0'
    write (unit_band, '(2f8.3)') 0.0_rp, 0.0_rp
    write (unit_band, '(2f8.3)') sk, 0.0_rp
    close (unit_band)
end subroutine build_band_path


subroutine build_band_path_weight(x, en_k,w_en_band_local,iband2io,na_band,no_max, nh, nx, nsl)
    use precision_mod
    use string_mod
    implicit none
    integer, intent(in) :: na_band,no_max,nh, nx, nsl
    integer, intent(in) :: iband2io(na_band)
    real(rp), intent(in) :: x(nx, 3), en_k(nh, nx, nsl)
    real(rp), intent(in) :: w_en_band_local(na_band,no_max,nx,nh,nsl)   
    integer :: ih, ix, isl,ia_band,io
    integer, parameter :: unit_band = 10
    real(rp) :: sk
    character(len=*), parameter :: dir = 'band/'
    character(len=*), parameter :: file_band_weight = 'band_weight.dat'
    character(len=80) :: fmt
    
    open (unit=unit_band, file=dir//file_band_weight, action='write')
    
    do isl = 1, nsl
        write (unit_band, *) '@# k   band (ev)'
        if(isl>1) then
            write (unit_band, *) 
            write (unit_band, *) 
        endif   
        do ia_band=1,na_band
            if(ia_band>1) then
                write (unit_band, *) 
                write (unit_band, *) 
            endif
            write (unit_band, *) '@# atom no', ia_band,'number of orbitals', iband2io(ia_band)
            sk = 0.0_rp
            fmt = trim('('//int2str(3+ iband2io(ia_band))//'f12.7'//')')         	
            do ih=1,nh
                write (unit_band, fmt) sk, en_k(ih, 1, isl),sum(w_en_band_local(ia_band,:,1,ih,isl)),&
                (w_en_band_local(ia_band,io,1,ih,isl),io=1,iband2io(ia_band))
            end do
            do ix = 2, nx
                sk = sk + sqrt(sum((x(ix, :) - x(ix - 1, :))**2))
                do ih=1,nh
                    write (unit_band, fmt) sk, en_k(ih, ix, isl),sum(w_en_band_local(ia_band,:,ix,ih,isl)),&
                    (w_en_band_local(ia_band,io,ix,ih,isl),io=1,iband2io(ia_band))
                end do
            end do
        end do
    end do
    !    write (unit_band, *) '@# k   ef=0'
    !    write (unit_band, '(2f8.3)') 0.0_rp, 0.0_rp
    !    write (unit_band, '(2f8.3)') sk, 0.0_rp
    close (unit_band)
end subroutine build_band_path_weight


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
