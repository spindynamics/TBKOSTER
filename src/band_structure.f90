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
!  band_structure.f90
!  TBKOSTER
module band_structure_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_mod
  use element_mod
  use energy_mod
  use hamiltonian_tb_mod
  use math_mod, only: i_unit, delta_function, rho2nm
  use mesh_mod
#if defined(OpenMP_Fortran_FOUND)
  use omp_lib
#endif
  use precision_mod, only: rp,ip
  use string_mod, only: TBKOSTER_flush, int2str, lower, real2str, sl, cmplx2str
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'nen', &
   'en_min', &
   'en_max', &
   'na_dos', &
   'ia' &
   ]

  type,public :: band_structure
    ! Units
    class(units),pointer :: u
    ! Elements
    class(element),pointer :: e
    ! Atom
    class(atom),pointer :: a
    ! Reciprocal space mesh
    class(mesh),pointer :: k
    ! Hamiltonian
    class(hamiltonian_tb),pointer :: h
    ! Energy
    class(energy),pointer :: en

 !> @defgroup Local_band_energy_weight Local band energy weight-related
    !> variables
    !> @{
!>   Projection type ; options: 'none' (default), 'site', 'spin',
    !> 'orbit', 'spin-orbite'
    character(len=11) :: proj
    ! Local band energy atomic site number
    integer :: na_band,i_min,i_max
    ! Local band energy atomic site index
    integer(ip),dimension(:),allocatable ::  i_band,ia_band
    ! Lowest and highest band index
    ! integer :: iband_min, iband_max
    ! Site band weight
    real(rp), dimension(:,:,:,:,:), allocatable :: w_band_site
     ! Spin band weight
    real(rp), dimension(:,:,:,:), allocatable :: w_band_spin
     ! orbital band weight
    real(rp), dimension(:,:,:,:), allocatable :: w_band_orb
    !> @}


  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: save_proj_band_site
    procedure :: save_proj_band_spin
    procedure :: save_proj_band_orbit
    procedure :: initialize
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type band_structure

  ! Constructor
  interface band_structure
    procedure :: constructor
  end interface band_structure

contains
  function constructor(en) result(obj)
    class(energy),target,intent(in) :: en
    type(band_structure) :: obj

    obj%u => en%u
    obj%e => en%e
    obj%a => en%a
    obj%k => en%k
    obj%h => en%h
    obj%en => en
  end function constructor

  subroutine destructor(obj)
    type(band_structure) :: obj

    if(allocated(obj%ia_band))        deallocate(obj%ia_band)
    if(allocated(obj%w_band_site))    deallocate(obj%w_band_site)
    if(allocated(obj%w_band_spin))    deallocate(obj%w_band_spin)
    if(allocated(obj%w_band_orb))     deallocate(obj%w_band_orb)
  end subroutine destructor

  subroutine save_proj_band_site(obj,ik,isl,v_k)
    ! INPUT
    class(band_structure),intent(inout) :: obj
    integer :: ik,isl
    complex(rp),intent(in),dimension(2,obj%h%nh,obj%h%nh) :: v_k
    ! LOCAL
    integer :: ia_band1,ia,ie,io,ispin,jspin,imat,jmat,jmat2,imat_ispin,imat_jspin,jj,nn
   
  if(obj%na_band>0) then
    select case(obj%a%ns)
    case(1,2)
     do jmat=1,obj%h%nh
        jmat2=isl+(jmat-1)*obj%a%ns
        jj=ik+(jmat-1)*obj%k%nx
        nn=ik+(jmat2-1)*obj%k%nx
          do ia_band1=1,obj%na_band
            ia = obj%ia_band(ia_band1)
            ie = obj%a%ia2ie(ia)
            do io=1,obj%e%no(ie)
              imat=obj%h%iaos2ih(ia,io,1)
              obj%w_band_site(jmat,ik,ia_band1,io,isl) &
               = real(v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat)))
            end do
          end do
      end do
   
    case(4)
      do jmat=1,obj%h%nh
        jj=ik+(jmat-1)*obj%k%nx
          do ia_band1=1,obj%na_band
            ia = obj%ia_band(ia_band1)
            ie = obj%a%ia2ie(ia)
            do io=1,obj%e%no(ie)
              do ispin=1,2
                imat_ispin=obj%h%iaos2ih(ia,io,ispin)
                do jspin=1,2
                  imat_jspin=obj%h%iaos2ih(ia,io,jspin)
                  obj%w_band_site(jmat,ik,ia_band1,io,isl) &
                 = real(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)))
                end do
              end do
            end do
          end do
      end do
    end select

  endif

  end subroutine save_proj_band_site

  subroutine save_proj_band_spin(obj,ik,isl,v_k)
    ! INPUT
    class(band_structure),intent(inout) :: obj
    integer :: ik,isl
    complex(rp),intent(in),dimension(2,obj%h%nh,obj%h%nh) :: v_k
    ! LOCAL
    integer :: ih,ia_band1,ia,ie,io,ispin,jspin,imat,jmat,imat_ispin,imat_jspin
    complex(rp), dimension(2, 2) :: rho
    real(rp) :: n
    real(rp), dimension(3) :: m_cart


    do ih=1,obj%i_max-obj%i_min+1
      jmat=obj%i_band(ih)
      do ia= 1, obj%a%na
        ie = obj%a%ia2ie(ia)
        rho = cmplx(0.0_rp, 0.0_rp, kind=rp)
        do ispin = 1, 2
        do jspin = 1, 2
           do io = 1, obj%e%no(ie)
              imat_ispin = obj%h%iaos2ih(ia, io, ispin)
              imat_jspin = obj%h%iaos2ih(ia, io, jspin)
                 rho(ispin, jspin) = rho(ispin, jspin)+ &
                    obj%a%g_s*(conjg(v_k(1, imat_jspin, jmat))*v_k(2, imat_ispin, jmat)     &
                              + v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)))/2
           end do
        end do
        end do
        call rho2nm(rho, n, m_cart)
        obj%w_band_spin(ih,ik,0,2:4) = obj%w_band_spin(ih,ik,0,2:4)+m_cart
     end do
       m_cart=obj%w_band_spin(jmat,ik,0,2:4)
       obj%w_band_spin(ih,ik,0,1) = norm2(m_cart)
   end do

    do ih=1,obj%i_max-obj%i_min+1
      jmat=obj%i_band(ih)
     do ia_band1 = 1, obj%na_band
       ia = obj%ia_band(ia_band1)
       ie = obj%a%ia2ie(ia)
       rho = cmplx(0.0_rp, 0.0_rp, kind=rp)
       do ispin = 1, 2
       do jspin = 1, 2
          do io = 1, obj%e%no(ie)
             imat_ispin = obj%h%iaos2ih(ia, io, ispin)
             imat_jspin = obj%h%iaos2ih(ia, io, jspin)
                rho(ispin, jspin) = rho(ispin, jspin)+ &
                   obj%a%g_s*(conjg(v_k(1, imat_jspin, jmat))*v_k(2, imat_ispin, jmat)     &
                             + v_k(1, imat_ispin, jmat)*conjg(v_k(2, imat_jspin, jmat)))/2
          end do
       end do
       end do
       call rho2nm(rho, n, m_cart)
       obj%w_band_spin(ih,ik,ia_band1,1) = norm2(m_cart)
       obj%w_band_spin(ih,ik,ia_band1,2:4) = m_cart
     end do
   end do

  end subroutine save_proj_band_spin


  subroutine save_proj_band_orbit(obj,ik,isl,v_k)
    use math_mod, only: L_x, L_y, L_z
    ! INPUT
    class(band_structure),intent(inout) :: obj
    integer :: ik,isl
    complex(rp),intent(in),dimension(2,obj%h%nh,obj%h%nh) :: v_k
    ! LOCAL
    integer :: ih,ia_band1,ia,ie,io1,io2,l,ispin,jspin,imat,jmat,kmat,imat_ispin,imat_jspin,ir
    complex(rp), dimension(3, 2, 2) :: rho
    real(rp) :: n
    real(rp), dimension(3) :: m_cart
    complex(rp), dimension(3,9, 9) :: LMAT

    do ir = 1, 3  ! x: ir=1 , y: i=2, z: ir=3
      select case (ir)
           case (1)
            LMAT(ir,:,:) = L_x(:,:)
          case (2)
            LMAT(ir,:,:) = L_y(:,:)
          case (3)
            LMAT(ir,:,:) = L_z(:,:)
      end select  
    end do

    do ih=1,obj%i_max-obj%i_min+1
      kmat=obj%i_band(ih)
      obj%w_band_orb(ih,ik,0,:)=0.0D0
      do ia = 1, obj%a%na
          ie = obj%a%ia2ie(ia)
        do ir = 1, 3  ! x: ir=1 , y: i=2, z: ir=3
          do io1 = 1, obj%e%no(ie)
             l = obj%e%o2l(obj%e%o(ie, io1))
                do ispin = 1, 2
                   imat = obj%h%iaos2ih(ia, io1, ispin)
                   do io2 = 1, obj%e%no(ie)
                      jmat = obj%h%iaos2ih(ia, io2, ispin)
                      obj%w_band_orb(ih,ik,0,ir+1) = obj%w_band_orb(ih,ik,0,ir+1) &
                           +obj%a%g_s*conjg(v_k(1, imat, kmat))*LMAT(ir,io1, io2)*v_k(2, jmat, kmat)
                   end do
                end do  
          end do  
        end do   
      end do   !fin de la boucle sur ia
       m_cart=obj%w_band_orb(ih,ik,0,2:4)
       obj%w_band_orb(ih,ik,0,1) = norm2(m_cart)
    end do     ! fin de la boucle sur kmat

    do ih=1,obj%i_max-obj%i_min+1
      kmat=obj%i_band(ih)
      do ia_band1 = 1, obj%na_band
        obj%w_band_orb(ih,ik,ia_band1,:)=0.0D0
          ia = obj%ia_band(ia_band1)
          ie = obj%a%ia2ie(ia)
        do ir = 1, 3  ! x: ir=1 , y: i=2, z: ir=3
          do io1 = 1, obj%e%no(ie)
             l = obj%e%o2l(obj%e%o(ie, io1))
                do ispin = 1, 2
                   imat = obj%h%iaos2ih(ia, io1, ispin)
                   do io2 = 1, obj%e%no(ie)
                      jmat = obj%h%iaos2ih(ia, io2, ispin)
                      obj%w_band_orb(ih,ik,ia_band1,ir+1) = obj%w_band_orb(ih,ik,ia_band1,ir+1) &
                           +obj%a%g_s*conjg(v_k(1, imat, kmat))*LMAT(ir,io1, io2)*v_k(2, jmat, kmat)
                   end do
                end do  
          end do  
        end do   
        m_cart=obj%w_band_orb(ih,ik,ia_band1,2:4)
       obj%w_band_orb(ih,ik,ia_band1,1) = norm2(m_cart)
      end do   !fin de la boucle sur ia_band
    end do     ! fin de la boucle sur kmat

  end subroutine save_proj_band_orbit

  subroutine initialize(obj)
    class(band_structure),intent(inout) :: obj
       ! Local band weight
    if(TRIM(obj%proj)=='site') then
      if(obj%na_band>0) then
         allocate(obj%w_band_site(obj%h%nh,obj%k%nx,obj%na_band,obj%e%no_max,obj%a%nsl))
      end if
    elseif(TRIM(obj%proj)=='spin') then
        allocate(obj%w_band_spin(obj%i_max-obj%i_min+1,obj%k%nx,0:obj%na_band,4))
    elseif(TRIM(obj%proj)=='orbit') then
        allocate(obj%w_band_orb(obj%i_max-obj%i_min+1,obj%k%nx,0:obj%na_band,4))
    elseif(TRIM(obj%proj)=='spin,orbit') then
        allocate(obj%w_band_spin(obj%i_max-obj%i_min+1,obj%k%nx,0:obj%na_band,4))
        allocate(obj%w_band_orb(obj%i_max-obj%i_min+1,obj%k%nx,0:obj%na_band,4))
    endif
  end subroutine initialize


  !> Read object in text format from file (default: 'in_dos.txt')
  subroutine read_txt(obj,file)
    class(band_structure),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    character(len=11) :: proj
    integer(ip) :: i_min,i_max,na_band,ih
    integer(ip),dimension(:),allocatable :: i_band,ia_band
    ! Namelist
    namelist /band/proj,i_min,i_max,na_band,ia_band
    
    call initialize_proj(proj)
    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_band.txt'
    end if
    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'band%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'band%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    na_band=0
    i_min=1
    i_max=obj%h%nh
    allocate(ia_band(0))
    read(10,nml=band,iostat=iostatus)
    obj%proj=proj
    call check_proj(proj,obj%a%ns)
    allocate(i_band(i_max-i_min+1))
    do ih=1,i_max-i_min+1
     i_band(ih)=ih+i_min-1
    end do
    deallocate(ia_band)
    allocate(ia_band(na_band))
    if(na_band > 0) then
      rewind(10)
      read(10,nml=band,iostat=iostatus)
    end if
    obj%i_min=i_min
    obj%i_max=i_max
    obj%na_band=na_band
    call move_alloc(ia_band,obj%ia_band)
    call move_alloc(i_band,obj%i_band)
    call obj%initialize()
    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

   !> Check the validity of the magnetic penalization
  subroutine check_proj(proj,ns)
    ! INPUT
    character(len=*) :: proj
    integer :: ns
    proj=TRIM(proj)
    if(proj=='spin'.AND.ns/=4) then
      write(error_unit,*) 'if band%proj is ""spin"" ns must be 4'
      error stop
    endif
    if(proj=='orbit'.AND.ns/=4) then
      write(error_unit,*) 'if band%proj is ""orbit"" ns must be 4'
      error stop
    endif
    if(proj=='spin,orbit'.AND.ns/=4) then
      write(error_unit,*) 'if band%proj is ""spin,orbit"" ns must be 4'
      error stop
    endif
    if(proj /= 'none' &
     .and. proj /= 'site' &
     .and. proj /= 'spin' &
     .and. proj /= 'orbit' &
     .and. proj /= 'spin,orbit') then
      write(error_unit,*) 'band%check_proj(): band%proj must be &
       &one of: ''none'', ''site'', ''spin'', ''orbit'', ''spin,orbit'''
      error stop
    end if
  end subroutine check_proj

   !> Initialize projection type
  subroutine initialize_proj(proj)
    character(len=11),intent(out) :: proj

    proj = 'none'
  end subroutine initialize_proj

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_dos.txt')
  subroutine write_txt(obj,file,unit)
    class(band_structure),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=11) :: proj
    integer(ip) :: i_min,i_max,na_band
    integer(ip),dimension(:),allocatable :: ia_band
    ! Namelist
    namelist /band/proj,i_min,i_max,na_band,ia_band

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_band.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if
    i_min=obj%i_min
    i_max=obj%i_max
    na_band = obj%na_band
    ia_band = obj%ia_band

    write(unit_rt,nml=band)
    call TBKOSTER_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_dos.txt'), if tag (default: .true.) the namelist opening and closing
  !> tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(band_structure),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
   ! Namelist variables
    real(rp),dimension(obj%h%nh,obj%k%nx,obj%a%nsl) :: en_k_2
    ! Local variables
    integer :: ip, isl, ik, ih,ia_band1,ia,io,ie, ir
    integer,dimension(:),allocatable         :: iband2io

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_band.txt'
    end if
    if(present(property)) then
      property_rt = property
    else
      property_rt = property_list
    end if
    if(present(tag)) then
      tag_rt = tag
    else
      tag_rt = .true.
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if
    if(tag_rt) then
      write(unit_rt,'(a)') '&band_out'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('en_k')
        en_k_2 = obj%en%en_k_2 * obj%u%convert_energy('from','hau')
        do isl=1,obj%a%nsl
          do ik=1,obj%k%nx
            do ih=1,obj%h%nh
              write(unit_rt,'(a)') ' en_k(' // int2str(ih) // ',' &
               // int2str(ik) // ',' // int2str(isl) // ') = ' &
               // real2str(en_k_2(ih,ik,isl))
            end do
          end do
        end do
      
       if(TRIM(obj%proj)=='site'.AND.obj%na_band>0) then
        allocate(iband2io(obj%na_band))
        do ia_band1=1,obj%na_band
          ia = obj%ia_band(ia_band1)
          ie = obj%a%ia2ie(ia)
          iband2io(ia_band1)=obj%e%no(ie) 
          write(unit_rt,'(a)') ' iband2io(' // int2str(ia_band1) //') = '//int2str(iband2io(ia_band1))
         end do
        
        do isl=1,obj%a%nsl
         do ih=1,obj%h%nh
           do ik=1,obj%k%nx
              do ia_band1=1,obj%na_band
              ia = obj%ia_band(ia_band1)
              ie = obj%a%ia2ie(ia)
              do io=1,obj%e%no(ie)
                   write(unit_rt,'(a)') ' w_band_site(' // int2str(ih) // ','// int2str(ik) // ','  &
                   // int2str(ia_band1) // ',' // int2str(io) // ',' // int2str(isl) // ') = ' &
                  // real2str(obj%w_band_site(ih,ik,ia_band1,io,isl))
                end do
             end do
           end do
         end do
        end do
        elseif(TRIM(obj%proj)=='spin') then    
           do ih=1,obj%i_max-obj%i_min+1
             do ik=1,obj%k%nx
              do ia_band1=0,obj%na_band
               do ir=1,4
                     write(unit_rt,'(a)') ' w_band_spin(' // int2str(ih) // ','// int2str(ik) //','   &
                     // int2str(ia_band1) //','// int2str(ir) // ') = ' // real2str(obj%w_band_spin(ih,ik,ia_band1,ir))
               end do
              end do
             end do
           end do
           elseif(TRIM(obj%proj)=='orbit') then    
            do ih=1,obj%i_max-obj%i_min+1
              do ik=1,obj%k%nx
                do ia_band1=0,obj%na_band
                 do ir=1,4
                  write(unit_rt,'(a)') ' w_band_orb(' // int2str(ih) // ','// int2str(ik) //','   &
                     // int2str(ia_band1) //','// int2str(ir) // ') = ' // real2str(obj%w_band_orb(ih,ik,ia_band1,ir))
                 end do
                end do
              end do
            end do
           elseif(TRIM(obj%proj)=='spin,orbit') then    
              do ih=1,obj%i_max-obj%i_min+1
                do ik=1,obj%k%nx
                 do ia_band1=0,obj%na_band
                   do ir=1,4
                        write(unit_rt,'(a)') ' w_band_spin(' // int2str(ih) // ','// int2str(ik) //','   &
                        // int2str(ia_band1) //','// int2str(ir) // ') = ' // real2str(obj%w_band_spin(ih,ik,ia_band1,ir))
                    end do
                  end do
                end do
               end do
               do ih=1,obj%i_max-obj%i_min+1
                do ik=1,obj%k%nx
                 do ia_band1=0,obj%na_band
                   do ir=1,4
                        write(unit_rt,'(a)') ' w_band_orb(' // int2str(ih) // ','// int2str(ik) //','   &
                        // int2str(ia_band1) //','// int2str(ir) // ') = ' // real2str(obj%w_band_orb(ih,ik,ia_band1,ir))
                    end do
                  end do
                end do
               end do
           end if
      end select
    end do

    if(tag_rt) then
      write(unit_rt,'(a)') ' /'
    end if
    if(.not. present(unit)) then
      call TBKOSTER_flush(unit_rt)
      close(unit_rt)
    else
      call TBKOSTER_flush(unit)
    end if
    !deallocate(file_rt,property_rt)
  end subroutine write_txt_formatted
end module band_structure_mod
