!
! Copyright (C) 2017
! Cyrille Barreteau <mailto:cyrille.barreteau@cea.fr>,
! Mathieu Cesar <mailto:mathieu.cesar@cea.fr>,
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
!  density_of_states.f90
!  DyNaMol
module density_of_states_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_mod
  use element_mod
  use energy_mod
  use hamiltonian_tb_mod
  use math_mod, only: i_unit, delta_function
  use mesh_mod
#if defined(OpenMP_Fortran_FOUND)
  use omp_lib
#endif
  use precision_mod, only: rp
  use string_mod, only: dynamol_flush, int2str, lower, real2str, sl, cmplx2str
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

  type,public :: density_of_states
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
    ! Number of energy points
    integer :: nen
    ! Energy range
    real(rp) :: en_min,en_max
    ! Energy delta
    real(rp) :: den
    ! Local DOS atomic site number
    integer :: na_dos
    ! Local DOS atomic site index
    integer,dimension(:),allocatable :: ia

    ! Total DOS
    real(rp),dimension(:,:),allocatable :: dos
    ! Local DOS
    complex(rp),dimension(:,:,:),allocatable :: dos_s
    complex(rp),dimension(:,:,:),allocatable :: dos_p,dos_px,dos_py,dos_pz
    complex(rp),dimension(:,:,:),allocatable :: dos_d,dos_dxy,dos_dyz,dos_dzx, &
     dos_dx2y2,dos_dz2r2

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: add_dos_k
    procedure :: add_dos_local_k
    procedure :: initialize
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type density_of_states

  ! Constructor
  interface density_of_states
    procedure :: constructor
  end interface density_of_states

contains
  function constructor(en) result(obj)
    class(energy),pointer,intent(in) :: en
    type(density_of_states) :: obj

    obj%u => en%u
    obj%e => en%e
    obj%a => en%a
    obj%k => en%k
    obj%h => en%h
    obj%en => en
  end function constructor

  subroutine destructor(obj)
    type(density_of_states) :: obj

    if(allocated(obj%ia))        deallocate(obj%ia)
    if(allocated(obj%dos))       deallocate(obj%dos)
    if(allocated(obj%dos_s))     deallocate(obj%dos_s)
    if(allocated(obj%dos_p))     deallocate(obj%dos_p)
    if(allocated(obj%dos_px))    deallocate(obj%dos_px)
    if(allocated(obj%dos_py))    deallocate(obj%dos_py)
    if(allocated(obj%dos_pz))    deallocate(obj%dos_pz)
    if(allocated(obj%dos_d))     deallocate(obj%dos_d)
    if(allocated(obj%dos_dxy))   deallocate(obj%dos_dxy)
    if(allocated(obj%dos_dyz))   deallocate(obj%dos_dyz)
    if(allocated(obj%dos_dzx))   deallocate(obj%dos_dzx)
    if(allocated(obj%dos_dx2y2)) deallocate(obj%dos_dx2y2)
    if(allocated(obj%dos_dz2r2)) deallocate(obj%dos_dz2r2)
  end subroutine destructor

  subroutine add_dos_k(obj,ik,isl)
    use, intrinsic :: iso_fortran_env, only: output_unit
    implicit none
    ! INPUT
    class(density_of_states),intent(inout) :: obj
    integer,intent(in) :: ik,isl
    ! LOCAL
    real(rp) :: en,ffd
    integer :: ien,jmat,jmat2,jj,nn

    select case(obj%a%ns)
    case(1,2)
      do ien=1,obj%nen
!        en = obj%en%en_f + obj%en_min + (ien-1)*obj%den
        en = obj%en_min + (ien-1)*obj%den
        do jmat = 1,obj%h%nh
          jmat2 = isl+(jmat-1)*obj%a%ns
          jj = ik + (jmat-1)*obj%k%nx
          nn = ik + (jmat2-1)*obj%k%nx
          ffd = delta_function(-(obj%en%en_k(nn)-en)/obj%en%degauss,obj%en%smearing) &
           /obj%en%degauss
          obj%dos(ien,isl) = obj%dos(ien,isl) + ffd*obj%k%w(ik)*obj%a%g_s
        end do
      end do ! end of ien loop
    case(4)
      do ien=1,obj%nen
!        en = obj%en%en_f + obj%en_min + (ien-1)*obj%den
         en = obj%en_min + (ien-1)*obj%den
        do jmat=1,obj%h%nh
          nn = ik+(jmat-1)*obj%k%nx
          ffd = delta_function(-(obj%en%en_k(nn)-en)/obj%en%degauss,obj%en%smearing) &
           /obj%en%degauss
          obj%dos(ien,1) = obj%dos(ien,1) + ffd*obj%k%w(ik)
        end do
      end do ! end of ien loop
    end select
  end subroutine add_dos_k

  subroutine add_dos_local_k(obj,ik,isl,v_k)
    ! INPUT
    class(density_of_states),intent(inout) :: obj
    integer,intent(in) :: ik,isl
    complex(rp),dimension(2,obj%h%nh,obj%h%nh),intent(in) :: v_k
    !	LOCAL
    real(rp) :: ffd,en
    complex(rp) :: w_s,w_p,w_d
    complex(rp) :: w_px,w_py,w_pz
    complex(rp) :: w_dxy,w_dyz,w_dzx,w_dx2y2,w_dz2r2
    integer :: ia,ie,io,ia_dos,ien,ispin,jspin,imat,jmat,jmat2,jj,nn
    integer :: ispin_up,ispin_dn,imat_ispin,imat_jspin

    if(obj%na_dos>0) then
      select case(obj%a%ns)
      case(1,2)

        do ia_dos=1,obj%na_dos
          ia = obj%ia(ia_dos)
          ie = obj%a%ia2ie(ia)
          do ien=1,obj%nen

!            en = obj%en%en_f + obj%en_min + (ien-1)*obj%den
          en = obj%en_min + (ien-1)*obj%den
            do jmat=1,obj%h%nh

              jmat2=isl+(jmat-1)*obj%a%ns
              jj=ik+(jmat-1)*obj%k%nx
              nn=ik+(jmat2-1)*obj%k%nx

              ffd = delta_function(-(obj%en%en_k(nn)-en)/obj%en%degauss, &
               obj%en%smearing)/obj%en%degauss

              do io=1,obj%e%no(ie)
                imat=obj%h%iaos2ih(ia,io,1)

                if(obj%e%o(ie,io)==1) then

                  w_s = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                  obj%dos_s(ia_dos,ien,isl) = obj%dos_s(ia_dos,ien,isl) &
                   + w_s*ffd*obj%k%w(ik)*obj%a%g_s

                elseif(obj%e%o(ie,io)>=2 .and. obj%e%o(ie,io)<=4) then
                  w_p = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                  obj%dos_p(ia_dos,ien,isl) = obj%dos_p(ia_dos,ien,isl) &
                   + w_p*ffd*obj%k%w(ik)*obj%a%g_s

                  if(obj%e%o(ie,io)==2) then

                    w_px = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                    obj%dos_px(ia_dos,ien,isl) = obj%dos_px(ia_dos,ien,isl) &
                     + w_px*ffd*obj%k%w(ik)*obj%a%g_s

                  elseif(obj%e%o(ie,io)==3) then

                    w_py = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                    obj%dos_py(ia_dos,ien,isl) = obj%dos_py(ia_dos,ien,isl) &
                     + w_py*ffd*obj%k%w(ik)*obj%a%g_s

                  elseif(obj%e%o(ie,io)==4) then

                    w_pz = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                    obj%dos_pz(ia_dos,ien,isl) = obj%dos_pz(ia_dos,ien,isl) &
                     + w_pz*ffd*obj%k%w(ik)*obj%a%g_s
                  end if

                elseif(obj%e%o(ie,io)>=5 .and. obj%e%o(ie,io)<=9) then
                  w_d = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                  obj%dos_d(ia_dos,ien,isl) = obj%dos_d(ia_dos,ien,isl) &
                   + w_d*ffd*obj%k%w(ik)*obj%a%g_s


                  if(obj%e%o(ie,io)==5) then

                    w_dxy = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                    obj%dos_dxy(ia_dos,ien,isl) = obj%dos_dxy(ia_dos,ien,isl) &
                     + w_dxy*ffd*obj%k%w(ik)*obj%a%g_s

                  elseif(obj%e%o(ie,io)==6) then

                    w_dyz = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                    obj%dos_dyz(ia_dos,ien,isl) = obj%dos_dyz(ia_dos,ien,isl) &
                     + w_dyz*ffd*obj%k%w(ik)*obj%a%g_s

                  elseif(obj%e%o(ie,io)==7) then

                    w_dzx = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                    obj%dos_dzx(ia_dos,ien,isl) = obj%dos_dzx(ia_dos,ien,isl) &
                     + w_dzx*ffd*obj%k%w(ik)*obj%a%g_s

                  elseif(obj%e%o(ie,io)==8) then

                    w_dx2y2 = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                    obj%dos_dx2y2(ia_dos,ien,isl) = obj%dos_dx2y2(ia_dos,ien,isl) &
                     + w_dx2y2*ffd*obj%k%w(ik)*obj%a%g_s

                  elseif(obj%e%o(ie,io)==9) then

                    w_dz2r2 = v_k(1,imat,jmat)*conjg(v_k(2,imat,jmat))
                    obj%dos_dz2r2(ia_dos,ien,isl) = obj%dos_dz2r2(ia_dos,ien,isl) &
                     + w_dz2r2*ffd*obj%k%w(ik)*obj%a%g_s
                  end if

                end if
              end do  !end of io loop

            end do  !end of jmat loop

          end do  !end of ien loop
        end do  !end of ia_dos loop

      case(4)

        ispin_up=1
        ispin_dn=2

        ! 1=s; 2=px; 3=py; 4=pz; 5=dxy; 6=dyz; 7=dzx; 8=dx^2-y^2; 9=d(3z^2-r^2)

        do ia_dos=1,obj%na_dos
          ia = obj%ia(ia_dos)
          ie = obj%a%ia2ie(ia)
          do ien=1,obj%nen

!            en= obj%en%en_f +obj%en_min+(ien-1)*obj%den
           en= obj%en_min+(ien-1)*obj%den
            do jmat=1,obj%h%nh

              nn=ik+(jmat-1)*obj%k%nx

              ffd=delta_function(-(obj%en%en_k(nn)-en)/obj%en%degauss,&
              obj%en%smearing)/obj%en%degauss

              do io=1,obj%e%no(ie)
                do ispin=1,2

                  imat_ispin=obj%h%iaos2ih(ia,io,ispin)

                  do jspin=1,2

                    imat_jspin=obj%h%iaos2ih(ia,io,jspin)

                    if(obj%e%o(ie,io)==1) then
                      w_s=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                       +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                      obj%dos_s(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                       =obj%dos_s(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                       +w_s*ffd*obj%k%w(ik)

                    elseif(obj%e%o(ie,io)>=2 .and. obj%e%o(ie,io)<=4) then
                      w_p=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                       +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                      obj%dos_p(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                       =obj%dos_p(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                       +w_p*ffd*obj%k%w(ik)

                      if(obj%e%o(ie,io)==2) then
                        w_px=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%dos_px(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         =obj%dos_px(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         +w_px*ffd*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==3) then
                        w_py=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%dos_py(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         =obj%dos_py(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         +w_py*ffd*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==3) then
                        w_pz=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                        +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%dos_pz(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         =obj%dos_pz(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         +w_pz*ffd*obj%k%w(ik)
                      end if

                    elseif(obj%e%o(ie,io)>=5 .and. obj%e%o(ie,io)<=9) then
                      w_d=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                       +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                      obj%dos_d(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                       =obj%dos_d(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                       +w_d*ffd*obj%k%w(ik)

                      if(obj%e%o(ie,io)==5) then
                        w_dxy=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                        +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%dos_dxy(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         =obj%dos_dxy(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         +w_dxy*ffd*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==6) then
                        w_dyz=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%dos_dyz(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         =obj%dos_dyz(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         +w_dyz*ffd*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==7) then
                        w_dzx=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%dos_dzx(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         =obj%dos_dzx(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         +w_dzx*ffd*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==8) then
                        w_dx2y2=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%dos_dx2y2(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         =obj%dos_dx2y2(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         +w_dx2y2*ffd*obj%k%w(ik)
                      elseif(obj%e%o(ie,io)==9) then
                        w_dz2r2=(v_k(1,imat_ispin,jmat)*conjg(v_k(2,imat_jspin,jmat)) &
                         +v_k(2,imat_jspin,jmat)*conjg(v_k(1,imat_ispin,jmat)))/2
                        obj%dos_dz2r2(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         =obj%dos_dz2r2(ia_dos,ien,obj%a%iss2is(ispin,jspin)) &
                         +w_dz2r2*ffd*obj%k%w(ik)
                      end if

                    end if

                  end do ! fin de la boucle sur ispin
                end do ! fin de la boucle sur jspin
              end do ! fin de la boucle sur io
            end do ! fin de la boucle sur jmat

          end do  ! fin de la boucle sur ien
        end do  ! fin de la boucle sur ia_dos

      end select

    end if
  end subroutine add_dos_local_k

  subroutine initialize(obj)
    class(density_of_states),intent(inout) :: obj

    ! Total DOS
    if(allocated(obj%dos)) deallocate(obj%dos)
    allocate(obj%dos(obj%nen,obj%a%nsl))

    obj%dos = 0.0_rp

    ! Local DOS
    if(obj%na_dos>0) then
      if(allocated(obj%dos_s)) deallocate(obj%dos_s)
      allocate(obj%dos_s(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_p)) deallocate(obj%dos_p)
      allocate(obj%dos_p(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_px)) deallocate(obj%dos_px)
      allocate(obj%dos_px(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_py)) deallocate(obj%dos_py)
      allocate(obj%dos_py(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_pz)) deallocate(obj%dos_pz)
      allocate(obj%dos_pz(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_d)) deallocate(obj%dos_d)
      allocate(obj%dos_d(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_dxy)) deallocate(obj%dos_dxy)
      allocate(obj%dos_dxy(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_dyz)) deallocate(obj%dos_dyz)
      allocate(obj%dos_dyz(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_dzx)) deallocate(obj%dos_dzx)
      allocate(obj%dos_dzx(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_dx2y2)) deallocate(obj%dos_dx2y2)
      allocate(obj%dos_dx2y2(obj%na_dos,obj%nen,obj%a%ns))
      if(allocated(obj%dos_dz2r2)) deallocate(obj%dos_dz2r2)
      allocate(obj%dos_dz2r2(obj%na_dos,obj%nen,obj%a%ns))

      obj%dos_s = 0.0_rp
      obj%dos_p  = 0.0_rp
      obj%dos_px = 0.0_rp
      obj%dos_py = 0.0_rp
      obj%dos_pz = 0.0_rp
      obj%dos_d     = 0.0_rp
      obj%dos_dxy   = 0.0_rp
      obj%dos_dyz   = 0.0_rp
      obj%dos_dzx   = 0.0_rp
      obj%dos_dx2y2 = 0.0_rp
      obj%dos_dz2r2 = 0.0_rp
    end if
  end subroutine initialize

  subroutine initialize_en(nen,en_min,en_max)
    integer,intent(out) :: nen
    real(rp),intent(out) :: en_min,en_max

    nen = 0
    en_min = 0.0_rp
    en_max = 0.0_rp
  end subroutine initialize_en

  !> Read object in text format from file (default: 'in_dos.txt')
  subroutine read_txt(obj,file)
    class(density_of_states),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    integer  :: nen
    real(rp) :: en_min, en_max
    integer  :: na_dos
    integer,dimension(:),allocatable :: ia
    ! Namelist
    namelist /dos/ nen, en_min, en_max, na_dos, ia

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_dos.txt'
    end if
    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'dos%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'dos%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    call initialize_en(nen,en_min,en_max)
    na_dos = 0
    allocate(ia(0))
    read(10,nml=dos,iostat=iostatus)
    deallocate(ia)
    allocate(ia(na_dos))
    if(na_dos > 0) then
      rewind(10)
      read(10,nml=dos)
    end if

    obj%nen = nen
    obj%en_min = en_min * obj%u%convert_energy('to','hau')
    obj%en_max = en_max * obj%u%convert_energy('to','hau')
    obj%na_dos = na_dos
    call move_alloc(ia,obj%ia)

    obj%den = (en_max-en_min)*obj%u%convert_energy('to','hau')/(nen-1)

    call obj%initialize()

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_dos.txt')
  subroutine write_txt(obj,file,unit)
    class(density_of_states),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    integer  :: nen
    real(rp) :: en_min, en_max
    integer  :: na_dos
    integer,dimension(obj%na_dos) :: ia
    ! Namelist
    namelist /dos/ nen, en_min, en_max, na_dos, ia

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_dos.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    nen = obj%nen
    en_min = obj%en_min * obj%u%convert_energy('from','hau')
    en_max = obj%en_max * obj%u%convert_energy('from','hau')
    na_dos = obj%na_dos
    ia = obj%ia

    write(unit_rt,nml=dos)
    call dynamol_flush(unit_rt)

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
    class(density_of_states),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Local variables
    integer :: ia_dos, ien, ip, is, isl

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_dos.txt'
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
      write(unit_rt,'(a)') '&dos'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('nen')
        write(unit_rt,'(a)') ' nen = ' // int2str(obj%nen)
      case('en_min')
        write(unit_rt,'(a)') ' en_min = ' // real2str(obj%en_min &
         * obj%u%convert_energy('from','hau'))
      case('en_max')
        write(unit_rt,'(a)') ' en_max = ' // real2str(obj%en_max &
         * obj%u%convert_energy('from','hau'))
      case('den')
        write(unit_rt,'(a)') ' den = ' // real2str(obj%den &
         * obj%u%convert_energy('from','hau'))
      case('na_dos')
        write(unit_rt,'(a)') ' na_dos = ' // int2str(obj%na_dos)
      case('ia')
        do ia_dos=1,obj%na_dos
          write(unit_rt,'(a)') ' ia(' // int2str(ia_dos) // ') = ' &
           // int2str(obj%ia(ia_dos))
        end do
      case('dos')
        do ien=1,obj%nen
          do isl=1,obj%a%nsl
            write(unit_rt,'(a)') ' dos_tot(' // int2str(ien) // ',' &
             // int2str(isl) // ') = ' // real2str(obj%dos(ien,isl))
          end do
        end do
      case('dos_s')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_s(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_s(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_p')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_p(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_p(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_px')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_px(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_px(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_py')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_py(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_py(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_pz')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_pz(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_pz(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_d')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_d(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_d(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_dxy')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_dxy(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_dxy(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_dyz')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_dyz(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_dyz(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_dzx')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_dzx(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_dzx(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_dx2y2')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_dx2y2(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_dx2y2(ia_dos,ien,is))
            end do
          end do
        end do
      case('dos_dz2r2')
        do ia_dos=1,obj%na_dos
          do ien=1,obj%nen
            do is=1,obj%a%ns
              write(unit_rt,'(a)') ' dos_dz2r2(' // int2str(ia_dos) // ',' &
               // int2str(ien) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%dos_dz2r2(ia_dos,ien,is))
            end do
          end do
        end do
      end select
    end do

    if(tag_rt) then
      write(unit_rt,'(a)') ' /'
    end if
    call dynamol_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt,property_rt)
  end subroutine write_txt_formatted
end module density_of_states_mod
