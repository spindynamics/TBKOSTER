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
!  density_of_states.f90
!  TBKOSTER
module magnetic_anisotropy_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_mod
  use element_mod
  use energy_mod
  use hamiltonian_tb_mod
  use math_mod, only: i_unit, delta_function,fermi_function,theta_function,deg2rad,rad2deg
  use mesh_mod
#if defined(OpenMP_Fortran_FOUND)
  use omp_lib
#endif
  use precision_mod, only: rp
  use string_mod, only: TBKOSTER_flush, int2str, lower, real2str, sl, cmplx2str
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'na_mae', &
   'ia', &
   'Eref', &
   'angle' &
   ]

  type,public :: magnetic_anisotropy
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
    ! number of local MAE atomic site number
    integer :: na_mae
    ! Local MAE atomic site index
    integer,dimension(:),allocatable :: ia
    ! number of spin angles 
    integer :: nangle
     ! Local MAE atomic site index
    real(rp), dimension(2) :: angle1,angle2
    ! reference energy (default value is Ef)
    real(rp) :: Eref
    ! Total MAE
    real(rp) :: mae_tot
    ! Local MAE
    complex(rp),dimension(:),allocatable :: mae_s
    complex(rp),dimension(:),allocatable :: mae_p,mae_px,mae_py,mae_pz
    complex(rp),dimension(:),allocatable :: mae_d,mae_dxy,mae_dyz,mae_dzx, &
     mae_dx2y2,mae_dz2r2

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: add_mae_k
    procedure :: add_mae_local_k
    procedure :: initialize
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type magnetic_anisotropy

  ! Constructor
  interface magnetic_anisotropy
    procedure :: constructor
  end interface magnetic_anisotropy

contains
  function constructor(en) result(obj)
    class(energy),target,intent(in) :: en
    type(magnetic_anisotropy) :: obj

    obj%u => en%u
    obj%e => en%e
    obj%a => en%a
    obj%k => en%k
    obj%h => en%h
    obj%en => en
  end function constructor

  subroutine destructor(obj)
    type(magnetic_anisotropy) :: obj

    if(allocated(obj%ia))        deallocate(obj%ia)
    if(allocated(obj%mae_s))     deallocate(obj%mae_s)
    if(allocated(obj%mae_p))     deallocate(obj%mae_p)
    if(allocated(obj%mae_px))    deallocate(obj%mae_px)
    if(allocated(obj%mae_py))    deallocate(obj%mae_py)
    if(allocated(obj%mae_pz))    deallocate(obj%mae_pz)
    if(allocated(obj%mae_d))     deallocate(obj%mae_d)
    if(allocated(obj%mae_dxy))   deallocate(obj%mae_dxy)
    if(allocated(obj%mae_dyz))   deallocate(obj%mae_dyz)
    if(allocated(obj%mae_dzx))   deallocate(obj%mae_dzx)
    if(allocated(obj%mae_dx2y2)) deallocate(obj%mae_dx2y2)
    if(allocated(obj%mae_dz2r2)) deallocate(obj%mae_dz2r2)
  end subroutine destructor

  subroutine add_mae_k(obj,ik,isl,Eref)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    ! INPUT
    class(magnetic_anisotropy),intent(inout) :: obj
    integer,intent(in) :: ik,isl
    real(rp),intent(in) :: Eref
    ! LOCAL
    real(rp) :: en,ff
    integer :: ien,jmat,jmat2,jj,nn

    select case(obj%a%ns)
    case(1,2)
      write(*,*) 'for mae ns should be 4'
      stop
    case(4)
        do jmat=1,obj%h%nh
          nn = ik+(jmat-1)*obj%k%nx
          ff  = theta_function(-(obj%en%en_k(nn)-Eref)/obj%en%degauss,obj%en%smearing)&
               *(obj%en%en_k(nn)-Eref)
          obj%mae_tot = obj%mae_tot + ff*obj%k%w(ik)
        end do
    end select
  end subroutine add_mae_k

  subroutine add_mae_local_k(obj,ik,isl,v_k,Eref)
    ! INPUT
    class(magnetic_anisotropy),intent(inout) :: obj
    integer,intent(in) :: ik,isl
    real(rp),intent(in) :: Eref
    complex(rp),dimension(2,obj%h%nh,obj%h%nh),intent(in) :: v_k
    !	LOCAL
    real(rp) :: ff,en
    complex(rp) :: w_s,w_p,w_d
    complex(rp) :: w_px,w_py,w_pz
    complex(rp) :: w_dxy,w_dyz,w_dzx,w_dx2y2,w_dz2r2
    integer :: ia,ie,io,ia_mae,ispin,jspin,imat,jmat,jmat2,jj,nn
    integer :: ispin_up,ispin_dn,imat_ispin,imat_jspin

    if(obj%na_mae>0) then
      select case(obj%a%ns)
      case(1,2)
        write(*,*) 'for mae ns should be 4'
        stop
      case(4)

        ispin_up=1
        ispin_dn=2
        
        ! 1=s; 2=px; 3=py; 4=pz; 5=dxy; 6=dyz; 7=dzx; 8=dx^2-y^2; 9=d(3z^2-r^2)

        do ia_mae=1,obj%na_mae
          ia = obj%ia(ia_mae)
          ie = obj%a%ia2ie(ia)
          
            do jmat=1,obj%h%nh

              nn=ik+(jmat-1)*obj%k%nx

              ff=theta_function(-(obj%en%en_k(nn)-Eref)/obj%en%degauss,obj%en%smearing)&
               *(obj%en%en_k(nn)-Eref)

              do io=1,obj%e%no(ie)
                do ispin=1,2

                  imat_ispin=obj%h%iaos2ih(ia,io,ispin)

                  !do jspin=1,2
                   jspin=ispin !diagonal elements 
                    imat_jspin=obj%h%iaos2ih(ia,io,jspin)

                    if(obj%e%o(ie,io)==1) then
                      w_s=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                       +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                      obj%mae_s(ia_mae)=obj%mae_s(ia_mae)+w_s*ff*obj%k%w(ik)

                    elseif(obj%e%o(ie,io)>=2 .and. obj%e%o(ie,io)<=4) then
                      w_p=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                       +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                      obj%mae_p(ia_mae)=obj%mae_p(ia_mae)+w_p*ff*obj%k%w(ik)

                      if(obj%e%o(ie,io)==2) then
                        w_px=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%mae_px(ia_mae)=obj%mae_px(ia_mae)+w_px*ff*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==3) then
                        w_py=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%mae_py(ia_mae)=obj%mae_py(ia_mae)+w_py*ff*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==4) then
                        w_pz=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                        +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%mae_pz(ia_mae)=obj%mae_pz(ia_mae)+w_pz*ff*obj%k%w(ik)
                      end if

                    elseif(obj%e%o(ie,io)>=5 .and. obj%e%o(ie,io)<=9) then
                      w_d=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                       +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                      obj%mae_d(ia_mae)=obj%mae_d(ia_mae)+w_d*ff*obj%k%w(ik)

                      if(obj%e%o(ie,io)==5) then
                        w_dxy=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                        +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%mae_dxy(ia_mae)=obj%mae_dxy(ia_mae)+w_dxy*ff*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==6) then
                        w_dyz=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%mae_dyz(ia_mae)=obj%mae_dyz(ia_mae)+w_dyz*ff*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==7) then
                        w_dzx=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%mae_dzx(ia_mae)=obj%mae_dzx(ia_mae)+w_dzx*ff*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==8) then
                        w_dx2y2=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%mae_dx2y2(ia_mae)=obj%mae_dx2y2(ia_mae)+w_dx2y2*ff*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==9) then
                        w_dz2r2=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%mae_dz2r2(ia_mae)=obj%mae_dz2r2(ia_mae)+w_dz2r2*ff*obj%k%w(ik)
                      end if

                    end if

                  !end do ! fin de la boucle sur jspin
                end do ! fin de la boucle sur ispin
              end do ! fin de la boucle sur io
            end do ! fin de la boucle sur jmat
        end do  ! fin de la boucle sur ia_dos

      end select

    end if
  end subroutine add_mae_local_k

  subroutine initialize(obj)
    class(magnetic_anisotropy),intent(inout) :: obj

    obj%mae_tot = 0.0_rp
    ! Local MAE
    if(obj%na_mae>0) then
      if(allocated(obj%mae_s)) deallocate(obj%mae_s)
      allocate(obj%mae_s(obj%na_mae))
      if(allocated(obj%mae_p)) deallocate(obj%mae_p)
      allocate(obj%mae_p(obj%na_mae))
      if(allocated(obj%mae_px)) deallocate(obj%mae_px)
      allocate(obj%mae_px(obj%na_mae))
      if(allocated(obj%mae_py)) deallocate(obj%mae_py)
      allocate(obj%mae_py(obj%na_mae))
      if(allocated(obj%mae_pz)) deallocate(obj%mae_pz)
      allocate(obj%mae_pz(obj%na_mae))
      if(allocated(obj%mae_d)) deallocate(obj%mae_d)
      allocate(obj%mae_d(obj%na_mae))
      if(allocated(obj%mae_dxy)) deallocate(obj%mae_dxy)
      allocate(obj%mae_dxy(obj%na_mae))
      if(allocated(obj%mae_dyz)) deallocate(obj%mae_dyz)
      allocate(obj%mae_dyz(obj%na_mae))
      if(allocated(obj%mae_dzx)) deallocate(obj%mae_dzx)
      allocate(obj%mae_dzx(obj%na_mae))
      if(allocated(obj%mae_dx2y2)) deallocate(obj%mae_dx2y2)
      allocate(obj%mae_dx2y2(obj%na_mae))
      if(allocated(obj%mae_dz2r2)) deallocate(obj%mae_dz2r2)
      allocate(obj%mae_dz2r2(obj%na_mae))

      obj%mae_s = 0.0_rp
      obj%mae_p  = 0.0_rp
      obj%mae_px = 0.0_rp
      obj%mae_py = 0.0_rp
      obj%mae_pz = 0.0_rp
      obj%mae_d     = 0.0_rp
      obj%mae_dxy   = 0.0_rp
      obj%mae_dyz   = 0.0_rp
      obj%mae_dzx   = 0.0_rp
      obj%mae_dx2y2 = 0.0_rp
      obj%mae_dz2r2 = 0.0_rp
    end if
  end subroutine initialize

  !> Read object in text format from file (default: 'in_dos.txt')
  subroutine read_txt(obj,file)
    class(magnetic_anisotropy),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    integer  :: na_mae, nangle
    integer  :: iia,iangle
    integer,dimension(:),allocatable :: ia
    real(rp) :: Eref   
    real(rp),dimension(2) :: angle1,angle2
    real(rp) :: en_f,en   
    ! Namelist
    namelist /mae/na_mae,nangle,ia,Eref,angle1,angle2
    namelist /energy/en_f,en 

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_mae.txt'
    end if
    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'mae%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'mae%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    na_mae = 0
    open(unit=11,file='out_energy.txt',action='read')
    read(unit=11,nml=energy)
    Eref=en_f
    close(unit=11)
    allocate(ia(0))
    read(10,nml=mae,iostat=iostatus)
    deallocate(ia)
    allocate(ia(na_mae))
    if(na_mae > 0.and.na_mae<obj%a%na) then
      rewind(10)
      read(10,nml=mae)
    elseif(na_mae==obj%a%na) then
      do iia=1,obj%a%na
          ia(iia)=iia
      end do
    end if
    
    obj%na_mae = na_mae
    obj%nangle=nangle
    obj%angle1=angle1*deg2rad
    obj%angle2=angle2*deg2rad
    obj%Eref=Eref*obj%u%convert_energy('to','hau')
    call move_alloc(ia,obj%ia)
    call obj%initialize()
    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_dos.txt')
  subroutine write_txt(obj,file,unit)
    class(magnetic_anisotropy),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
     ! Namelist variables
    integer  :: na_mae, nangle
    integer  :: iia,iangle
    integer,dimension(:),allocatable :: ia
    real(rp) :: Eref   
    real(rp),dimension(2) :: angle1,angle2
    ! Namelist
    namelist /mae/na_mae,nangle, ia, Eref,angle1,angle2

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_mae.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    na_mae = obj%na_mae
    ia = obj%ia
    nangle=obj%nangle
    Eref=obj%Eref
    write(*,*) 'Eref= ',Eref
    angle1=obj%angle1*rad2deg 
    angle2=obj%angle2*rad2deg
    

    write(unit_rt,nml=mae)
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
    class(magnetic_anisotropy),intent(in) :: obj
     character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Local variables
    integer :: ia_mae, iangle,ip
    real(rp) :: mae_sum

   if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_mae.txt'
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
      write(unit_rt,'(a)') '&mae_out'
    end if

 do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('Eref')
         write(unit_rt,'(a)') ' Eref = ' // real2str(obj%Eref)
      case('angle')
        write(unit_rt,'(a)') ' nangle = ' // int2str(obj%nangle)
        write(unit_rt,'(a)') ' angle1 = ' // real2str(obj%angle1(1)*rad2deg) // ' , '&
                                          // real2str(obj%angle1(2)*rad2deg)
        write(unit_rt,'(a)') ' angle2 = ' // real2str(obj%angle2(1)*rad2deg) // ' , ' &
                                          // real2str(obj%angle2(2)*rad2deg)
        case('na_mae')
        write(unit_rt,'(a)') ' na_mae = ' // int2str(obj%na_mae)
      case('ia')
        do ia_mae=1,obj%na_mae
          write(unit_rt,'(a)') ' ia(' // int2str(ia_mae) // ') = ' &
           // int2str(obj%ia(ia_mae))
        end do
      case('mae')
        mae_sum= 0.0_rp

        do ia_mae=1,obj%na_mae
         
              write(unit_rt,'(a)') ' mae_s(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_s(ia_mae)*obj%u%convert_energy('from','hau'))
               
              write(unit_rt,'(a)') ' mae_p(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_p(ia_mae)*obj%u%convert_energy('from','hau'))
        
              write(unit_rt,'(a)') ' mae_px(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_px(ia_mae)*obj%u%convert_energy('from','hau'))
   
              write(unit_rt,'(a)') ' mae_py(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_py(ia_mae)*obj%u%convert_energy('from','hau'))
    
              write(unit_rt,'(a)') ' mae_pz(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_pz(ia_mae)*obj%u%convert_energy('from','hau'))
   
              write(unit_rt,'(a)') ' mae_d(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_d(ia_mae)*obj%u%convert_energy('from','hau'))
    
              write(unit_rt,'(a)') ' mae_dxy(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_dxy(ia_mae)*obj%u%convert_energy('from','hau'))
   
              write(unit_rt,'(a)') ' mae_dyz(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_dyz(ia_mae)*obj%u%convert_energy('from','hau'))
 
              write(unit_rt,'(a)') ' mae_dzx(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_dzx(ia_mae)*obj%u%convert_energy('from','hau'))
    
              write(unit_rt,'(a)') ' mae_dx2y2(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_dx2y2(ia_mae)*obj%u%convert_energy('from','hau'))
   
              write(unit_rt,'(a)') ' mae_dz2r2(' // int2str(ia_mae) // ') = ' &
               // cmplx2str(obj%mae_dz2r2(ia_mae)*obj%u%convert_energy('from','hau'))

               mae_sum=mae_sum+obj%mae_s(ia_mae)+obj%mae_p(ia_mae)+obj%mae_d(ia_mae)
        end do
        write(unit_rt,'(a)') ' mae_sum= ' // real2str(mae_sum*obj%u%convert_energy('from','hau'))
        write(unit_rt,'(a)') ' mae_tot= ' // real2str(obj%mae_tot*obj%u%convert_energy('from','hau'))
       end select
    end do   
 
    if(tag_rt) then
      write(unit_rt,'(a)') ' /'
    end if
    call TBKOSTER_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt,property_rt)
  end subroutine write_txt_formatted
end module magnetic_anisotropy_mod
