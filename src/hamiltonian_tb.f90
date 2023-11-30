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
!  hamiltonian_tb.f90
!  TBKOSTER
module hamiltonian_tb_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_tb_mod
#if defined(BLAS95_FOUND)
  use blas95, only: hemm
#endif
  use charge_mod
  use element_tb_mod
#if defined(LAPACK95_FOUND)
  use lapack95, only: hegv
#endif
  use math_mod, only: epsilon, sigma_x, sigma_y, sigma_z, cart2sph, rho2nm, &
   sph2cart
  use mesh_mod
  use precision_mod, only: rp
  use string_mod, only: sl, cmplx2str, int2str, lower, real2str, TBKOSTER_flush
  use units_mod
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'e_e_interaction', &
   'm_penalization' &
   ]

  type,public :: hamiltonian_tb
    ! Units
    class(units),pointer :: u
    ! Elements TB
    class(element_tb),pointer :: e_tb
    ! Atoms TB
    class(atom_tb),pointer :: a_tb
    ! Charges
    class(charge),pointer :: q
    ! Reciprocal space mesh
    class(mesh),pointer :: k

    !> Electron-electron interaction ; options: 'stoner' (default), 'ujb'
    character(len=6) :: e_e_interaction
    !> Magnetic penalization type ; options: 'none' (default), 'r', 'r,theta',
    !> 'r,theta,phi', 'theta', 'theta,phi', 'phi'
    character(len=11) :: m_penalization

    ! Hamiltonian matrix dimension
    integer :: nh
    ! Atom/orbital/spin-to-Hamiltonian index
    integer,dimension(:,:,:),allocatable :: iaos2ih
    ! Intra-atomic energy (diagonal elements of the Hamiltonian)
    real(rp),dimension(:,:),allocatable :: en_intra
    ! Real space Hamiltonian
    real(rp),dimension(:,:,:,:),allocatable :: h_r
    ! Real space overlap
    real(rp),dimension(:,:,:,:),allocatable :: s_r
    ! Hamiltonian renormalizations: EEI, LCN, magnetic penalization
    complex(rp),dimension(:,:,:,:),allocatable :: delta_h_eei
    complex(rp),dimension(:,:,:)  ,allocatable :: delta_v_lcn
    complex(rp),dimension(:,:,:)  ,allocatable :: delta_v_pen
    ! Reciprocal space Hamiltonian
    !complex(rp),dimension(:,:),allocatable :: h_k
    ! Reciprocal space overlap
    !complex(rp),dimension(:,:),allocatable :: s_k

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: add_delta_h_eei
    procedure :: add_delta_h_ov
    procedure :: add_delta_h_so
    !procedure :: add_delta_subsystem
    procedure :: build_projection_k
    procedure :: build_d_projection_k
    procedure :: build_v_k
    procedure :: build_w_k
    procedure :: calculate_b_pen
    procedure :: calculate_delta_h_eei
    procedure :: calculate_delta_v_lcn
    procedure :: calculate_delta_v_pen
    procedure :: calculate_h_r
    procedure :: calculate_s_r
    procedure :: initialize
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type hamiltonian_tb

  ! Constructor
  interface hamiltonian_tb
    procedure :: constructor
  end interface hamiltonian_tb

contains
  function constructor(a_tb,q,k) result(obj)
    class(atom_tb),target,intent(in) :: a_tb
    class(charge),target,intent(in) :: q
    class(mesh),target,intent(in) :: k
    type(hamiltonian_tb) :: obj

    obj%u => a_tb%u
    obj%e_tb => a_tb%e_tb
    obj%a_tb => a_tb
    obj%q => q
    obj%k => k

    call obj%initialize()
  end function constructor

  subroutine destructor(obj)
    type(hamiltonian_tb) :: obj

    if(allocated(obj%delta_v_pen)) deallocate(obj%delta_v_pen)
    if(allocated(obj%delta_v_lcn)) deallocate(obj%delta_v_lcn)
    if(allocated(obj%delta_h_eei)) deallocate(obj%delta_h_eei)
    if(allocated(obj%s_r))         deallocate(obj%s_r)
    if(allocated(obj%h_r))         deallocate(obj%h_r)
    if(allocated(obj%en_intra))    deallocate(obj%en_intra)
    if(allocated(obj%iaos2ih))     deallocate(obj%iaos2ih)
  end subroutine destructor

  subroutine add_delta_h_eei(obj,isl,h_k)
    ! INPUT
    class(hamiltonian_tb),intent(in) :: obj
    integer,intent(in) :: isl
    complex(rp),dimension(obj%nh,obj%nh),intent(inout) :: h_k
    ! LOCAL
    integer :: ia,ie,io1,io2,ispin,jspin,imat,jmat
#if defined(DEBUG)
    write(output_unit,*) 'DEBUG == Entering add_delta_h_eei'
#endif

    select case(obj%a_tb%ns)
    case(1,2)
      do ia=1,obj%a_tb%na
        ie = obj%a_tb%ia2ie(ia)
        do io1=1,obj%e_tb%no(ie)
          imat = obj%iaos2ih(ia,io1,isl)
          do io2=1,obj%e_tb%no(ie)
            jmat = obj%iaos2ih(ia,io2,isl)
            h_k(imat,jmat) = h_k(imat,jmat) + obj%delta_h_eei(ia,io1,io2,isl)
          end do
        end do
      end do
    case(4)
      do ia=1,obj%a_tb%na
        ie = obj%a_tb%ia2ie(ia)
        do io1=1,obj%e_tb%no(ie)
          do io2=1,obj%e_tb%no(ie)
            do ispin=1,2
              imat = obj%iaos2ih(ia,io1,ispin)
              do jspin=1,2
                jmat = obj%iaos2ih(ia,io2,jspin)
                h_k(imat,jmat) = h_k(imat,jmat) &
                 + obj%delta_h_eei(ia,io1,io2,obj%a_tb%iss2is(ispin,jspin))
              end do
            end do
          end do
        end do
      end do
    end select
  end subroutine add_delta_h_eei

  subroutine add_delta_h_ov(obj,isl,h_k,s_k)
    ! INPUT
    class(hamiltonian_tb),intent(in) :: obj
    integer,intent(in) :: isl
    complex(rp),dimension(obj%nh,obj%nh),intent(in) :: s_k
    complex(rp),dimension(obj%nh,obj%nh),intent(inout) :: h_k
    ! LOCAL
    integer :: ia1,ia2,ie1,ie2,io1,io2,l1,l2,ispin,jspin
    integer :: imat,jmat,imat_ispin,imat_jspin,jmat_ispin,jmat_jspin
    complex(rp),dimension(:,:,:), allocatable :: delta_v_ov
    complex(rp),dimension(:,:), allocatable  :: delta_h_ov
#if defined(DEBUG)
    write(output_unit,*) 'DEBUG == Entering add_delta_h_ov'
#endif

    if (allocated(delta_v_ov)) deallocate(delta_v_ov)
    if (allocated(delta_h_ov)) deallocate(delta_h_ov)
    allocate(delta_v_ov(obj%a_tb%na,3,obj%a_tb%ns),delta_h_ov(obj%nh,obj%nh))

    delta_v_ov = obj%delta_v_lcn + obj%delta_v_pen
    delta_h_ov = cmplx(0.0_rp,0.0_rp,kind=rp)

    select case(obj%a_tb%ns)
    case(1,2)
      do ia1=1,obj%a_tb%na
        ie1 = obj%a_tb%ia2ie(ia1)
        do io1=1,obj%e_tb%no(ie1)
          imat = obj%iaos2ih(ia1,io1,1)
          l1 = obj%e_tb%o2l(obj%e_tb%o(ie1,io1))
          do ia2=1,obj%a_tb%na
            ie2 = obj%a_tb%ia2ie(ia2)
            do io2=1,obj%e_tb%no(ie2)
              jmat = obj%iaos2ih(ia2,io2,1)
              l2 = obj%e_tb%o2l(obj%e_tb%o(ie2,io2))

              delta_h_ov(imat,jmat) = 0.5_rp*(delta_v_ov(ia1,l1,isl) &
               + delta_v_ov(ia2,l2,isl))*s_k(imat,jmat)
            end do
          end do
        end do
      end do
    case(4)
      do ispin=1,2
        do jspin=1,2
          do ia1=1,obj%a_tb%na
            ie1 = obj%a_tb%ia2ie(ia1)
            do io1=1,obj%e_tb%no(ie1)
              imat_ispin = obj%iaos2ih(ia1,io1,ispin)
              imat_jspin = obj%iaos2ih(ia1,io1,jspin)
              l1 = obj%e_tb%o2l(obj%e_tb%o(ie1,io1))

              do ia2=1,obj%a_tb%na
                ie2 = obj%a_tb%ia2ie(ia2)
                do io2=1,obj%e_tb%no(ie2)
                  jmat_ispin = obj%iaos2ih(ia2,io2,ispin)
                  jmat_jspin = obj%iaos2ih(ia2,io2,jspin)
                  l2 = obj%e_tb%o2l(obj%e_tb%o(ie2,io2))

                  delta_h_ov(imat_ispin,jmat_jspin) = 0.5_rp &
                   *(delta_v_ov(ia1,l1,obj%a_tb%iss2is(ispin,jspin)) &
                   + delta_v_ov(ia2,l2,obj%a_tb%iss2is(ispin,jspin))) &
                   *s_k(imat_ispin,jmat_ispin)
                end do
              end do
            end do
          end do
        end do
      end do
    end select

    h_k = h_k + delta_h_ov
    if (allocated(delta_v_ov)) deallocate(delta_v_ov)
    if (allocated(delta_h_ov)) deallocate(delta_h_ov)
  end subroutine add_delta_h_ov

  subroutine add_delta_h_so(obj,h_k)
    ! INPUT
    class(hamiltonian_tb),intent(in) :: obj
    complex(rp),dimension(obj%nh,obj%nh),intent(inout) :: h_k
    ! LOCAL
    integer :: ia,ie,io1,io2,l1,l2,is1,is2,ih1,ih2,ils1,ils2

    do ia=1,obj%a_tb%na
      ie = obj%a_tb%ia2ie(ia)
      do io1=1,obj%e_tb%no(ie)
        l1 = obj%e_tb%o2l(obj%e_tb%o(ie,io1))
        do io2=1,obj%e_tb%no(ie)
          l2 = obj%e_tb%o2l(obj%e_tb%o(ie,io2))
          if(l1==l2) then
            do is1=1,2
              ih1 = obj%iaos2ih(ia,io1,is1)
              ils1 = is1+(io1-1)*2
              do is2=1,2
                ih2 = obj%iaos2ih(ia,io2,is2)
                ils2 = is2+(io2-1)*2
                h_k(ih1,ih2) = h_k(ih1,ih2) + obj%e_tb%l_dot_s(ie,l1,ils1,ils2)
              end do
            end do
          end if
        end do
      end do
    end do
  end subroutine add_delta_h_so

  !> Calculate the k-space projection M(k) of the matrix M(r)
  function build_projection_k(obj,m_r,c_k) result(m_k)
    ! INPUT
    class(hamiltonian_tb),intent(in) :: obj
    real(rp),dimension(obj%a_tb%na,0:obj%a_tb%nn_max,obj%e_tb%no_max,obj%e_tb%no_max), intent(in) :: m_r
    complex(rp),dimension(obj%a_tb%na,obj%a_tb%nn_max,obj%a_tb%nsp), intent(in) :: c_k
    ! OUTPUT
    complex(rp),dimension(obj%nh,obj%nh) :: m_k
    ! LOCAL
    integer :: ispin,ia1,ia2,in,ie1,ie2,io1,io2,imat,jmat
#if defined(DEBUG)
    write(output_unit,*) 'DEBUG == Entering build_projection_k'
#endif

    m_k = cmplx(0.0_rp,0.0_rp,kind=rp)

    select case(obj%a_tb%ns)
    case(1,2)
      do ia1=1,obj%a_tb%na
        ie1 = obj%a_tb%ia2ie(ia1)
        do ia2=1,obj%a_tb%na
          ie2 = obj%a_tb%ia2ie(ia2)
          do io1=1,obj%e_tb%no(ie1)
            imat = obj%iaos2ih(ia1,io1,1)
            do io2=1,obj%e_tb%no(ie2)
              jmat = obj%iaos2ih(ia2,io2,1)
              do in=1,obj%a_tb%nn(ia1)
                if(obj%a_tb%ian2ia(ia1,in)==ia2) then
                  m_k(imat,jmat) = m_k(imat,jmat) + m_r(ia1,in,io1,io2)*c_k(ia1,in,1)
                end if
              end do
            end do
          end do
        end do
      end do

      do ia1=1,obj%a_tb%na
        ie1 = obj%a_tb%ia2ie(ia1)
        do io1=1,obj%e_tb%no(ie1)
          imat = obj%iaos2ih(ia1,io1,1)
          do io2=1,obj%e_tb%no(ie1)
            jmat = obj%iaos2ih(ia1,io2,1)
            m_k(imat,jmat) = m_k(imat,jmat) + m_r(ia1,0,io1,io2)
          end do
        end do
      end do

    case(4)
      do ispin=1,2

        do ia1=1,obj%a_tb%na
          ie1 = obj%a_tb%ia2ie(ia1)
          do ia2=1,obj%a_tb%na
            ie2 = obj%a_tb%ia2ie(ia2)
            do io1=1,obj%e_tb%no(ie1)
              imat = obj%iaos2ih(ia1,io1,ispin)
              do io2=1,obj%e_tb%no(ie2)
                jmat = obj%iaos2ih(ia2,io2,ispin)
                do in=1,obj%a_tb%nn(ia1)
                  if(obj%a_tb%ian2ia(ia1,in)==ia2) then
                    m_k(imat,jmat) = m_k(imat,jmat) &
                     + m_r(ia1,in,io1,io2)*c_k(ia1,in,ispin)
                  end if
                end do
              end do
            end do
          end do
        end do

        do ia1=1,obj%a_tb%na
          ie1 = obj%a_tb%ia2ie(ia1)
          do io1=1,obj%e_tb%no(ie1)
            imat = obj%iaos2ih(ia1,io1,ispin)
            do io2=1,obj%e_tb%no(ie1)
              jmat = obj%iaos2ih(ia1,io2,ispin)
              m_k(imat,jmat) = m_k(imat,jmat) + m_r(ia1,0,io1,io2)
            end do
          end do
        end do

      end do
    end select
  end function build_projection_k

  !> Calculate the k-space projection M'(k) of the matrix M'(r)
  function build_d_projection_k(obj,d_m_r,c_k) result(d_m_k)
    ! INPUT
    class(hamiltonian_tb),intent(in) :: obj
    real(rp),dimension(obj%a_tb%na,0:obj%a_tb%nn_max,obj%e_tb%no_max, &
     obj%e_tb%no_max,3),intent(in) :: d_m_r
    complex(rp),dimension(obj%a_tb%na,obj%a_tb%nn_max,obj%a_tb%nsp), &
     intent(in) :: c_k
    ! OUTPUT
    complex(rp),dimension(obj%nh,obj%nh,3,obj%a_tb%na) :: d_m_k
    ! LOCAL
    integer :: ispin,ia1,ia2,in,ie1,ie2,io1,io2,imat,jmat,ix

    d_m_k = cmplx(0.0_rp,0.0_rp,kind=rp)

    select case(obj%a_tb%ns)
    case(1,2)
    do ix=1,3
      do ia1=1,obj%a_tb%na
        ie1 = obj%a_tb%ia2ie(ia1)
        do ia2=1,obj%a_tb%na
          ie2 = obj%a_tb%ia2ie(ia2)
          do io1=1,obj%e_tb%no(ie1)
            imat = obj%iaos2ih(ia1,io1,1)
            do io2=1,obj%e_tb%no(ie2)
              jmat = obj%iaos2ih(ia2,io2,1)
              d_m_k(imat,jmat,ix,ia1) = cmplx(0.0_rp,0.0_rp,kind=rp) 
              d_m_k(jmat,imat,ix,ia1) = cmplx(0.0_rp,0.0_rp,kind=rp) 
              do in=1,obj%a_tb%nn(ia1)
                if(obj%a_tb%ian2ia(ia1,in)==ia2) then
                  d_m_k(imat,jmat,ix,ia1) = d_m_k(imat,jmat,ix,ia1) &
                  + d_m_r(ia1,in,io1,io2,ix)*c_k(ia1,in,1)
                end if
              end do
              d_m_k(jmat,imat,ix,ia1) = conjg(d_m_k(imat,jmat,ix,ia1))
            end do
          end do
        end do
      end do
    end do
     
    !     do ix=1,3
    !      do ia1=1,obj%a_tb%na
    !        ie1 = obj%a_tb%ia2ie(ia1)
    !        do io1=1,obj%e_tb%no(ie1)
    !          imat = obj%iaos2ih(ia1,io1,1)
    !          do io2=1,obj%e_tb%no(ie1)
    !            jmat = obj%iaos2ih(ia1,io2,1)
    !            d_m_k(imat,jmat,ix,ia1) = d_m_k(imat,jmat,ix,ia1) + d_m_r(ia1,0,io1,io2,ix)
    !          end do
    !        end do
    !      end do
    !     end do

    case(4)
    do ix=1,3
      do ispin=1,2
 
        do ia1=1,obj%a_tb%na
          ie1 = obj%a_tb%ia2ie(ia1)
          do ia2=1,obj%a_tb%na
            ie2 = obj%a_tb%ia2ie(ia2)
            do io1=1,obj%e_tb%no(ie1)
              imat = obj%iaos2ih(ia1,io1,ispin)
              do io2=1,obj%e_tb%no(ie2)
                jmat = obj%iaos2ih(ia2,io2,ispin)
                d_m_k(imat,jmat,ix,ia1) = cmplx(0.0_rp,0.0_rp,kind=rp)
                d_m_k(jmat,imat,ix,ia1) = cmplx(0.0_rp,0.0_rp,kind=rp)
                do in=1,obj%a_tb%nn(ia1)
                  if(obj%a_tb%ian2ia(ia1,in)==ia2) then
                    d_m_k(imat,jmat,ix,ia1) = d_m_k(imat,jmat,ix,ia1) &
                     + d_m_r(ia1,in,io1,io2,ix)*c_k(ia1,in,ispin)
                  end if
                end do
                d_m_k(jmat,imat,ix,ia1) = conjg(d_m_k(imat,jmat,ix,ia1))
              end do
            end do
          end do
        end do
      end do
    end do
!        do ia1=1,obj%a_tb%na
!          ie1 = obj%a_tb%ia2ie(ia1)
!          do io1=1,obj%e_tb%no(ie1)
!            imat = obj%iaos2ih(ia1,io1,ispin)
!            do io2=1,obj%e_tb%no(ie1)
!              jmat = obj%iaos2ih(ia1,io2,ispin)
!              d_m_k(imat,jmat,ix,ia1) = d_m_k(imat,jmat,ix,ia1) + d_m_r(ia1,0,io1,io2,ix)
!            end do
!          end do
!        end do
!      end do
!    end do
    end select
  end function build_d_projection_k

  function build_v_k(obj,ik,isl) result(v_k)
    class(hamiltonian_tb),intent(inout) :: obj
    ! INPUT
    integer,intent(in) :: ik,isl
    ! OUTPUT
    complex(rp),dimension(2,obj%nh,obj%nh) :: v_k
    ! LOCAL
    complex(rp),dimension(:,:,:), allocatable :: c_k
    complex(rp),dimension(:,:), allocatable   :: s_k, s_k_work
    real(rp),dimension(:), allocatable        :: w_k
    complex(rp),dimension(:,:), allocatable   :: v_k1,v_k2
    real(rp),dimension(3) :: k_point ! a k-point

#if !defined(LAPACK95_FOUND)
    complex(rp),dimension(:), allocatable :: work
    real(rp),dimension(:), allocatable :: rwork
    integer :: info
    integer :: lwork
    lwork = max(1,2*obj%nh-1)
    allocate(work(lwork),rwork(max(1,3*obj%nh-2)))
#endif
#if defined(DEBUG)
    write(output_unit,*) 'DEBUG == Entering build_v_k'
    write(output_unit,'(I5,1X,F10.7,1X,F10.7,1X,F10.7)') ik,obj%k%x(ik,:)
#endif
    if (.not.allocated(c_k))      allocate(c_k(obj%a_tb%na,obj%a_tb%nn_max,obj%a_tb%nsp))
    if (.not.allocated(s_k))      allocate(s_k(obj%nh,obj%nh))
    if (.not.allocated(s_k_work)) allocate(s_k_work(obj%nh,obj%nh))
    if (.not.allocated(w_k))      allocate(w_k(obj%nh))
    if (.not.allocated(v_k1))     allocate(v_k1(obj%nh,obj%nh))
    if (.not.allocated(v_k2))     allocate(v_k2(obj%nh,obj%nh))
#if defined(DEBUG)
    write(output_unit,*) "DEBUG ==> 1"
#endif
    call TBKOSTER_flush(output_unit)

    ! Build reciprocal space projections
    k_point(:) = obj%k%x(ik,:)
    c_k = obj%a_tb%build_c_k(k_point)
    v_k1 = obj%build_projection_k(obj%h_r,c_k)
    s_k =  obj%build_projection_k(obj%s_r,c_k)
    s_k_work = s_k

#if defined(DEBUG)
    write(output_unit,*) "DEBUG ==> 2"
    call TBKOSTER_flush(output_unit)
#endif    

    !	Add renormalization
    call obj%add_delta_h_eei(isl,v_k1)
    call obj%add_delta_h_ov(isl,v_k1,s_k)
    if(obj%a_tb%ns==4) then
      call obj%add_delta_h_so(v_k1)
    end if
    !call obj%add_delta_h_subsystem()

    !	Build the eigenvectors
#if defined(LAPACK95_FOUND)
    call hegv(v_k1,s_k_work,w_k,1,'V','U')
#else
    call zhegv(1,'V','U',obj%nh,v_k1,obj%nh,s_k_work,obj%nh,w_k,work, &
     lwork,rwork,info)
    deallocate(rwork,work)
#endif

    ! Build the eigenvectors tilde
#if defined(BLAS95_FOUND)
    call hemm(s_k,v_k1,v_k2)
#else
    call zhemm('L','U',obj%nh,obj%nh,cmplx(1.0_rp,0.0_rp,kind=rp),s_k,obj%nh, &
     v_k1,obj%nh,cmplx(0.0_rp,0.0_rp,kind=rp),v_k2,obj%nh)
#endif
    v_k(1,:,:) = v_k1(:,:)
    v_k(2,:,:) = v_k2(:,:)
  
    deallocate(v_k2,v_k1,w_k,s_k_work,s_k,c_k)
#if defined(DEBUG)  
    write(output_unit,*) 'DEBUG == Leaving build_v_k'
    call TBKOSTER_flush(output_unit)
#endif

  end function build_v_k

  function build_w_k(obj,ik,isl) result(w_k)
    ! INPUT
    class(hamiltonian_tb),intent(in) :: obj
    integer,intent(in) :: ik,isl
    ! OUTPUT
    real(rp),dimension(obj%nh) :: w_k
    ! LOCAL
    complex(rp),dimension(:,:,:), allocatable :: c_k
    complex(rp),dimension(:,:), allocatable :: h_k, s_k
    real(rp), dimension(3) :: k_point ! a k-point

#if !defined(LAPACK95_FOUND)
    complex(rp),dimension(:), allocatable:: work
    real(rp),dimension(:), allocatable :: rwork
    integer :: lwork
    integer :: info
    lwork = max(1,2*obj%nh-1)
    allocate(work(lwork),rwork(max(1,3*obj%nh-2)))
#endif
#if defined(DEBUG)
    write(output_unit,*) 'DEBUG == Entering build_w_k'
    write(output_unit,'(I5,1X,F10.7,1X,F10.7,1X,F10.7)') ik,obj%k%x(ik,:)
#endif
    if (.not.allocated(c_k)) allocate(c_k(obj%a_tb%na,obj%a_tb%nn_max,obj%a_tb%nsp))
    if (.not.allocated(h_k)) allocate(h_k(obj%nh,obj%nh))
    if (.not.allocated(s_k)) allocate(s_k(obj%nh,obj%nh))
    ! Build reciprocal space projections
    k_point(:)=obj%k%x(ik,:)
    c_k = obj%a_tb%build_c_k(k_point)
    h_k = obj%build_projection_k(obj%h_r,c_k)
    s_k = obj%build_projection_k(obj%s_r,c_k)

    !	Add renormalization
    call obj%add_delta_h_ov(isl,h_k,s_k)
    call obj%add_delta_h_eei(isl,h_k)
    if(obj%a_tb%ns==4) then
      call obj%add_delta_h_so(h_k)
    end if
    !call obj%add_delta_h_subsystem()

    !	Build the eigenvalues
#if defined(LAPACK95_FOUND)
    call hegv(h_k,s_k,w_k,1,'N','U')
#else
    call zhegv(1,'N','U',obj%nh,h_k,obj%nh,s_k,obj%nh,w_k,work,lwork,rwork,info)
    deallocate(work,rwork)
#endif
  deallocate(s_k,h_k,c_k)
  
  end function build_w_k

  subroutine calculate_b_pen(obj)
    ! INPUT
    class(hamiltonian_tb),intent(inout) :: obj
    ! LOCAL
    integer :: ia,l,ispin,jspin,ispin_rev,sigma
    real(rp) :: m_z, m_r
    real(rp),dimension(3) :: m_cart, m_sph, b_pen_1, b_pen_2

   !write(output_unit,*) "====> Entering calculate_b_pen"

    obj%a_tb%b_pen = 0.0_rp

    select case(obj%a_tb%ns)
    case(2)
      select case(obj%m_penalization)
      case('none')
        ! Do nothing
      case('r')
        do ia=1,obj%a_tb%na
          m_z = 0.0_rp
          do l=1,3
            m_z = m_z + obj%q%q_mul_in(ia,l,0) - obj%q%q_mul_in(ia,l,1)
          end do

          obj%a_tb%b_pen(ia,3) = -2.0_rp*obj%a_tb%lambda_pen(ia) &
           *(m_z-obj%a_tb%m_pen(ia,1))
        end do

      end select
    case(4)
      select case(obj%m_penalization)
      case('none')
        ! Do nothing
      case('r')
        do ia=1,obj%a_tb%na
          m_cart = sum(obj%q%q_mul_in(ia,:,1:3),1)
          m_r = norm2(m_cart)

          if(m_r>epsilon.and.obj%a_tb%lambda_pen(ia)>epsilon) then
            obj%a_tb%b_pen(ia,:) = -2.0_rp*obj%a_tb%lambda_pen(ia)*m_cart &
             *(m_r-obj%a_tb%m_pen(ia,1))/m_r
          else
            obj%a_tb%b_pen(ia,:) = 0.0_rp
          end if
        end do
      case('r,theta')
        do ia=1,obj%a_tb%na
          m_cart = sum(obj%q%q_mul_in(ia,:,1:3),1)
          m_r = norm2(m_cart)
          m_sph = cart2sph(m_cart)

          if(m_r>epsilon.and.obj%a_tb%lambda_pen(ia)>epsilon) then
            b_pen_1 = -2.0_rp*obj%a_tb%lambda_pen(ia)*m_cart &
             *(m_r-obj%a_tb%m_pen(ia,1))/m_r
          else
            b_pen_1 = 0.0_rp
          end if

          if(m_r>epsilon .and. abs(m_cart(3)/m_r)<1-epsilon.and.obj%a_tb%lambda_pen(ia)>epsilon) then
            b_pen_2(1) = -m_cart(1)*m_cart(3)
            b_pen_2(2) = -m_cart(2)*m_cart(3)
            b_pen_2(3) = (m_cart(1)**2+m_cart(2)**2)
            b_pen_2 = 2.0_rp*obj%a_tb%lambda_pen(ia) &
             *(m_sph(2)-obj%a_tb%m_pen(ia,2)) &
             /(sqrt(1.0_rp-(m_cart(3)/m_r)**2)*m_r**3)*b_pen_2
          else
            b_pen_2 = 0.0_rp
          end if

          obj%a_tb%b_pen(ia,:) = b_pen_1 + b_pen_2
        end do
      case('r,theta,phi')
        do ia=1,obj%a_tb%na
          m_cart = sum(obj%q%q_mul_in(ia,:,1:3),1)

          obj%a_tb%b_pen(ia,:) = -2.0_rp*obj%a_tb%lambda_pen(ia) &
           *(m_cart - sph2cart(obj%a_tb%m_pen(ia,:)))
        end do
      case('theta')
        do ia=1,obj%a_tb%na
          m_cart = sum(obj%q%q_mul_in(ia,:,1:3),1)
          m_r = norm2(m_cart)
          m_sph = cart2sph(m_cart)

          if(m_r>epsilon .and. abs(m_cart(3)/m_r)<1-epsilon.and.obj%a_tb%lambda_pen(ia)>epsilon) then
            obj%a_tb%b_pen(ia,1) = -m_cart(1)*m_cart(3)
            obj%a_tb%b_pen(ia,2) = -m_cart(2)*m_cart(3)
            obj%a_tb%b_pen(ia,3) =  m_cart(1)**2+m_cart(2)**2
            obj%a_tb%b_pen(ia,:) = 2.0_rp*obj%a_tb%lambda_pen(ia) &
             *(m_sph(2)-obj%a_tb%m_pen(ia,2)) &
             /(sqrt(1.0_rp-(m_cart(3)/m_r)**2)*m_r**3)*obj%a_tb%b_pen(ia,:)
          else
            obj%a_tb%b_pen(ia,:) = 0.0_rp
          end if
        end do
      case('theta,phi')
        do ia=1,obj%a_tb%na
          m_cart = sum(obj%q%q_mul_in(ia,:,1:3),1)
          m_r = norm2(m_cart)

         if(m_r>epsilon.and.obj%a_tb%lambda_pen(ia)>epsilon) then
          obj%a_tb%b_pen(ia,:) = -obj%a_tb%lambda_pen(ia) &
           *(m_cart/m_r - sph2cart(obj%a_tb%m_pen(ia,:))/obj%a_tb%m_pen(ia,1))
	  else
           obj%a_tb%b_pen(ia,:) = 0.0_rp	  
	  endif
        end do
      case('phi')
        do ia=1,obj%a_tb%na
          m_cart = sum(obj%q%q_mul_in(ia,:,1:3),1)
          m_r = sqrt(m_cart(1)**2+m_cart(2)**2)
          m_sph = cart2sph(m_cart)

          if(m_r>epsilon .and. abs(m_cart(1)/m_r)<1.0_rp-epsilon.and.obj%a_tb%lambda_pen(ia)>epsilon) then
            obj%a_tb%b_pen(ia,1) = -m_cart(2)*m_cart(2)
            obj%a_tb%b_pen(ia,2) =  m_cart(1)*m_cart(2)
            obj%a_tb%b_pen(ia,3) =  0.0_rp
            obj%a_tb%b_pen(ia,:) = 2.0_rp*obj%a_tb%lambda_pen(ia) &
             *(m_sph(3)-obj%a_tb%m_pen(ia,3)) &
             /(sqrt(1.0_rp-(m_cart(1)/m_r)**2)*m_r**3)*obj%a_tb%b_pen(ia,:)
          else
            obj%a_tb%b_pen(ia,:) = 0.0_rp
          end if
          if(m_cart(2)<0.0_rp) then
            obj%a_tb%b_pen(ia,:) = -obj%a_tb%b_pen(ia,:)
          end if
        end do
      end select
    end select
  end subroutine calculate_b_pen

  subroutine calculate_delta_h_eei(obj)
    ! INPUT
    class(hamiltonian_tb),intent(inout) :: obj
    ! LOCAL
    complex(rp),dimension(2,2) :: rho
    real(rp) :: n
    real(rp) :: m_z
    real(rp),dimension(3) :: m_cart
    real(rp) :: i_stoner_d
    real(rp),dimension(3) :: i_stoner
    real(rp) :: test, n_d
    integer :: ia,ie,io,l,ispin,jspin,ispin_rev
    integer :: io1,io2,io3,io4,l1,l2,l3,l4,o1,o2,o3,o4

    obj%delta_h_eei = cmplx(0.0_rp,0.0_rp,kind=rp)

    select case(trim(obj%e_e_interaction))
    case('stoner')
      select case(obj%a_tb%ns)
      case(2)
        do ispin=1,2
          ispin_rev = 3-ispin
          do ia=1,obj%a_tb%na
            ie = obj%a_tb%ia2ie(ia)
            i_stoner(1) = obj%e_tb%i_stoner_d(ie)/10
            i_stoner(2) = obj%e_tb%i_stoner_d(ie)/10
            i_stoner(3) = obj%e_tb%i_stoner_d(ie)

            m_z = 0.0_rp
            do io=1,obj%e_tb%no(ie)
              if(obj%e_tb%o2l(obj%e_tb%o(ie,io))==3) then
                m_z = m_z + real(obj%q%rho_net_in(ia,io,io,ispin) &
                 - obj%q%rho_net_in(ia,io,io,ispin_rev))
              end if
            end do

            do io=1,obj%e_tb%no(ie)
              l = obj%e_tb%o2l(obj%e_tb%o(ie,io))
              obj%delta_h_eei(ia,io,io,ispin) = obj%delta_h_eei(ia,io,io,ispin) &
               - i_stoner(l)/2*m_z
            end do
          end do
        end do
      case(4)
        do ia=1,obj%a_tb%na  ! loop on the atoms of the system
          ie = obj%a_tb%ia2ie(ia)
          i_stoner(1) = obj%e_tb%i_stoner_d(ie)/10
          i_stoner(2) = obj%e_tb%i_stoner_d(ie)/10
          i_stoner(3) = obj%e_tb%i_stoner_d(ie)

          rho = 0.0_rp
          do ispin=1,2
            do jspin=1,2
              do io=1,obj%e_tb%no(ie)
                if(obj%e_tb%o2l(obj%e_tb%o(ie,io))==3) then
                  rho(ispin,jspin) = rho(ispin,jspin) &
                   + obj%q%rho_net_in(ia,io,io,obj%a_tb%iss2is(ispin,jspin))
                end if
              end do
            end do
          end do
          call rho2nm(rho,n,m_cart)
          do ispin=1,2
            do jspin=1,2
              do io=1,obj%e_tb%no(ie) ! loop on the orbitals
                obj%delta_h_eei(ia,io,io,obj%a_tb%iss2is(ispin,jspin))&
                 = obj%delta_h_eei(ia,io,io,obj%a_tb%iss2is(ispin,jspin))&
                 - 0.5_rp*i_stoner(obj%e_tb%o2l(obj%e_tb%o(ie,io)))&
                 *(m_cart(1)*sigma_x(ispin,jspin)&
                 + m_cart(2)*sigma_y(ispin,jspin)&
                 + m_cart(3)*sigma_z(ispin,jspin))
              end do
            end do
          end do
        end do
      end select

    case('ujb')
      select case(obj%a_tb%ns)
      case(1,2)

        do ispin=1,obj%a_tb%ns
          ispin_rev = obj%a_tb%ns+1-ispin

          do ia=1,obj%a_tb%na
            ie = obj%a_tb%ia2ie(ia)
            do io1=1,obj%e_tb%no(ie)
              do io2=1,obj%e_tb%no(ie)
                l1 = obj%e_tb%o2l(obj%e_tb%o(ie,io1))
                l2 = obj%e_tb%o2l(obj%e_tb%o(ie,io2))

                if(l1==3 .and. l2==3) then
                  o1 = obj%e_tb%o(ie,io1)-4
                  o2 = obj%e_tb%o(ie,io2)-4
                  do io3=1,obj%e_tb%no(ie)
                    l3 = obj%e_tb%o2l(obj%e_tb%o(ie,io3))
                    do io4=1,obj%e_tb%no(ie)
                      l4 = obj%e_tb%o2l(obj%e_tb%o(ie,io4))
                      if(l3==3 .and. l4==3) then
                        o3 = obj%e_tb%o(ie,io3)-4
                        o4 = obj%e_tb%o(ie,io4)-4
                        obj%delta_h_eei(ia,io1,io2,ispin) &
                         = obj%delta_h_eei(ia,io1,io2,ispin)&
                         + (obj%e_tb%v_osc(ie,o3,o1,o4,o2)&
                         - obj%e_tb%v_osc(ie,o3,o1,o2,o4))&
                         *obj%q%rho_net_in(ia,io4,io3,ispin)&
                         + obj%e_tb%v_osc(ie,o3,o1,o4,o2)&
                         *obj%q%rho_net_in(ia,io4,io3,ispin_rev)
                      end if
                    end do
                  end do

                end if

              end do
            end do
          end do

          do ia=1,obj%a_tb%na
            ie = obj%a_tb%ia2ie(ia)

            n_d = 0.0_rp
            do io1=1,obj%e_tb%no(ie)
              if(obj%e_tb%o2l(obj%e_tb%o(ie,io1))==3) then
                n_d = n_d + real(obj%q%rho_net_in(ia,io1,io1,ispin) &
                 + obj%q%rho_net_in(ia,io1,io1,ispin_rev))
              end if
            end do

            do io1=1,obj%e_tb%no(ie)
              if(obj%e_tb%o2l(obj%e_tb%o(ie,io1))==3) then
                obj%delta_h_eei(ia,io1,io1,ispin) = obj%delta_h_eei(ia,io1,io1,ispin) &
                 - (9*obj%e_tb%u_dd(ie) - 2*obj%e_tb%j_dd(ie))/10*n_d
              end if
            end do

          end do

          do ia=1,obj%a_tb%na  !
            ie = obj%a_tb%ia2ie(ia)
            i_stoner_d = (obj%e_tb%u_dd(ie)+6*obj%e_tb%j_dd(ie))/5.0D0
            i_stoner(1) = i_stoner_d/10
            i_stoner(2) = i_stoner_d/10
            i_stoner(3) = i_stoner_d

            m_z = 0.0_rp
            do io1=1,obj%e_tb%no(ie)
              if(obj%e_tb%o2l(obj%e_tb%o(ie,io1))==3) then
                m_z = m_z + real(obj%q%rho_net_in(ia,io1,io1,ispin) &
                 - obj%q%rho_net_in(ia,io1,io1,ispin_rev))
              end if
            end do

            do io1=1,obj%e_tb%no(ie)
              if(obj%e_tb%o2l(obj%e_tb%o(ie,io1))==1 &
               .or. obj%e_tb%o2l(obj%e_tb%o(ie,io1))==2) then
                obj%delta_h_eei(ia,io1,io1,ispin) = obj%delta_h_eei(ia,io1,io1,ispin)&
                 - i_stoner(obj%e_tb%o2l(obj%e_tb%o(ie,io1)))/2*m_z
              end if
            end do
          end do

        end do

      case(4)
        !
        ! upup and dndn
        !
        do ia=1,obj%a_tb%na
          ie = obj%a_tb%ia2ie(ia)
          do io1=1,obj%e_tb%no(ie)
            l1 = obj%e_tb%o2l(obj%e_tb%o(ie,io1))
            do io2=1,obj%e_tb%no(ie)
              do ispin=1,2
                ispin_rev = 3-ispin
                l2 = obj%e_tb%o2l(obj%e_tb%o(ie,io2))
                if(l1==3 .and. l2==3) then
                  o1 = obj%e_tb%o(ie,io1)-4
                  o2 = obj%e_tb%o(ie,io2)-4
                  do io3=1,obj%e_tb%no(ie)
                    l3 = obj%e_tb%o2l(obj%e_tb%o(ie,io3))
                    do io4=1,obj%e_tb%no(ie)
                      l4 = obj%e_tb%o2l(obj%e_tb%o(ie,io4))
                      if(l3==3 .and. l4==3) then
                        o3 = obj%e_tb%o(ie,io3)-4
                        o4 = obj%e_tb%o(ie,io4)-4
                        obj%delta_h_eei(ia,io1,io2,obj%a_tb%iss2is(ispin,ispin))&
                         = obj%delta_h_eei(ia,io1,io2,obj%a_tb%iss2is(ispin,ispin))&
                         + (obj%e_tb%v_osc(ie,o3,o1,o4,o2)&
                         - obj%e_tb%v_osc(ie,o3,o1,o2,o4))&
                         *obj%q%rho_net_in(ia,io4,io3,obj%a_tb%iss2is(ispin,ispin))&
                         + obj%e_tb%v_osc(ie,o3,o1,o4,o2)&
                         *obj%q%rho_net_in(ia,io4,io3,obj%a_tb%iss2is(ispin_rev,ispin_rev))
                      end if
                    end do
                  end do

                end if

              end do
            end do
          end do
        end do

        do ia=1,obj%a_tb%na
          ie = obj%a_tb%ia2ie(ia)

          rho = 0.0_rp
          do ispin=1,2
            do jspin=1,2
              do io1=1,obj%e_tb%no(ie)
                if(obj%e_tb%o2l(obj%e_tb%o(ie,io1))==3) then
                  rho(ispin,jspin) = rho(ispin,jspin) &
                   + obj%q%rho_net_in(ia,io1,io1,obj%a_tb%iss2is(ispin,jspin))
                end if
              end do
            end do
          end do
          call rho2nm(rho,n_d,m_cart)
          do io1=1,obj%e_tb%no(ie)
            do ispin=1,2
              if(obj%e_tb%o2l(obj%e_tb%o(ie,io1))==3) then
                obj%delta_h_eei(ia,io1,io1,obj%a_tb%iss2is(ispin,ispin)) &
                 = obj%delta_h_eei(ia,io1,io1,obj%a_tb%iss2is(ispin,ispin)) &
                 - (9*obj%e_tb%u_dd(ie)-2*obj%e_tb%j_dd(ie))/10*n_d
              end if
            end do
          end do

        end do

        !
        ! updn and dnup
        !
        do ia=1,obj%a_tb%na
          ie = obj%a_tb%ia2ie(ia)
          do io1=1,obj%e_tb%no(ie)
            l1 = obj%e_tb%o2l(obj%e_tb%o(ie,io1))
            do io2=1,obj%e_tb%no(ie)
              do ispin=1,2
                ispin_rev = 3-ispin
                l2 = obj%e_tb%o2l(obj%e_tb%o(ie,io2))
                if(l1==3 .and. l2==3) then
                  o1 = obj%e_tb%o(ie,io1)-4
                  o2 = obj%e_tb%o(ie,io2)-4
                  do io3=1,obj%e_tb%no(ie)
                    l3 = obj%e_tb%o2l(obj%e_tb%o(ie,io3))
                    do io4=1,obj%e_tb%no(ie)
                      l4 = obj%e_tb%o2l(obj%e_tb%o(ie,io4))
                      if(l3==3 .and. l4==3) then
                        o3 = obj%e_tb%o(ie,io3)-4
                        o4 = obj%e_tb%o(ie,io4)-4

                        obj%delta_h_eei(ia,io1,io2,obj%a_tb%iss2is(ispin,ispin_rev))&
                         = obj%delta_h_eei(ia,io1,io2,obj%a_tb%iss2is(ispin,ispin_rev))&
                         - obj%e_tb%v_osc(ie,o3,o1,o2,o4)&
                         *obj%q%rho_net_in(ia,io4,io3,obj%a_tb%iss2is(ispin,ispin_rev))
                      end if
                    end do
                  end do

                end if

              end do
            end do
          end do
        end do

        !
        ! take care of s and p
        !
        do ia=1,obj%a_tb%na  ! loop on the atoms of the system
          i_stoner_d = (obj%e_tb%u_dd(ie)+6*obj%e_tb%j_dd(ie))/5.0D0
          i_stoner(1) = i_stoner_d/10
          i_stoner(2) = i_stoner_d/10
          i_stoner(3) = i_stoner_d

          rho = 0.0_rp
          do ispin=1,2
            do jspin=1,2
              do io1=1,obj%e_tb%no(ie)
                if(obj%e_tb%o2l(obj%e_tb%o(ie,io1))==3) then
                  rho(ispin,jspin) = rho(ispin,jspin)&
                   + obj%q%rho_net_in(ia,io1,io1,obj%a_tb%iss2is(ispin,jspin))
                end if
              end do
            end do
          end do
          call rho2nm(rho,n_d,m_cart)
          do ispin=1,2
            do jspin=1,2
              do io1=1,obj%e_tb%no(ie)
                l1 = obj%e_tb%o2l(obj%e_tb%o(ie,io1))
                if(l1==1 .or. l1==2) then
                  obj%delta_h_eei(ia,io1,io1,obj%a_tb%iss2is(ispin,jspin))&
                   = obj%delta_h_eei(ia,io1,io1,obj%a_tb%iss2is(ispin,jspin))&
                   - 0.5*i_stoner(obj%e_tb%o2l(obj%e_tb%o(ie,io1)))&
                   *(m_cart(1)*sigma_x(ispin,jspin)&
                   + m_cart(2)*sigma_y(ispin,jspin)&
                   + m_cart(3)*sigma_z(ispin,jspin))
                end if

              end do
            end do
          end do

        end do

        do ia=1,obj%a_tb%na
          ie = obj%a_tb%ia2ie(ia)
          do io1=1,obj%e_tb%no(ie)
            do io2=1,obj%e_tb%no(ie)
              do ispin=1,2
                do jspin=1,2
                  test = abs(obj%delta_h_eei(ia,io1,io2,obj%a_tb%iss2is(ispin,jspin))&
                   -conjg(obj%delta_h_eei(ia,io2,io1,obj%a_tb%iss2is(jspin,ispin))))
                  if(test>epsilon) then
                    write(*,*) 'Non Hermitian obj%delta_h_eei matrix',io1,io2,ispin,jspin
                  end if

                  test = abs(obj%q%rho_net_in(ia,io1,io2,obj%a_tb%iss2is(ispin,jspin))&
                   -conjg(obj%q%rho_net_in(ia,io2,io1,obj%a_tb%iss2is(jspin,ispin))))
                  if(test>epsilon) then
                    write(*,*) 'Non Hermitian rho_net matrix',io1,io2,ispin,jspin
                  end if

                end do
              end do
            end do
          end do
        end do

      end select
    end select
  end subroutine calculate_delta_h_eei

  subroutine calculate_delta_v_lcn(obj)
    ! INPUT
    class(hamiltonian_tb),intent(inout) :: obj
    ! LOCAL
    real(rp) :: dn,dn_d
    integer :: ia,ie,l,ispin,ispin_rev

    if (.not.allocated(obj%delta_v_lcn)) allocate(obj%delta_v_lcn(obj%a_tb%na,3,obj%a_tb%ns))

    obj%delta_v_lcn = cmplx(0.0_rp,0.0_rp,kind=rp)

    select case(obj%a_tb%ns)
    case(1)
      ispin=1
      do ia=1,obj%a_tb%na
        ie = obj%a_tb%ia2ie(ia)

        dn = 0.0_rp
        do l=1,3
          dn = dn + obj%q%q_mul_in(ia,l,ispin-1)
        end do
        dn = dn - obj%e_tb%q(ie)
        dn_d = obj%q%q_mul_in(ia,3,ispin-1) - obj%e_tb%q_d(ie)

        do l=1,3
          obj%delta_v_lcn(ia,l,ispin) = obj%delta_v_lcn(ia,l,ispin) &
           + obj%e_tb%u_lcn(ie)*dn
        end do
        obj%delta_v_lcn(ia,3,ispin) = obj%delta_v_lcn(ia,3,ispin) &
         + obj%e_tb%u_lcn_d(ie)*dn_d
      end do
    case(2)
      do ispin=1,2
        ispin_rev = 3-ispin
        do ia=1,obj%a_tb%na
          ie = obj%a_tb%ia2ie(ia)

          dn = 0.0_rp
          do l=1,3
            dn = dn + obj%q%q_mul_in(ia,l,ispin-1) &
             + obj%q%q_mul_in(ia,l,ispin_rev-1)
          end do
          dn = dn - obj%e_tb%q(ie)
          dn_d = obj%q%q_mul_in(ia,3,ispin-1) &
           + obj%q%q_mul_in(ia,3,ispin_rev-1) - obj%e_tb%q_d(ie)

          do l=1,3
            obj%delta_v_lcn(ia,l,ispin) = obj%delta_v_lcn(ia,l,ispin) &
             + obj%e_tb%u_lcn(ie)*dn
          end do
          obj%delta_v_lcn(ia,3,ispin) = obj%delta_v_lcn(ia,3,ispin) &
           + obj%e_tb%u_lcn_d(ie)*dn_d
        end do
      end do
    case(4)
      do ia=1,obj%a_tb%na  ! loop on the atoms of the system
        ie = obj%a_tb%ia2ie(ia)
        dn = 0.0_rp
        do l=1,3
          dn = dn + obj%q%q_mul_in(ia,l,0)
        end do
        dn = dn - obj%e_tb%q(ie)
        dn_d = obj%q%q_mul_in(ia,3,0) - obj%e_tb%q_d(ie)
        ! U delta n Id
        do l=1,3
          obj%delta_v_lcn(ia,l,1) = obj%delta_v_lcn(ia,l,1) + obj%e_tb%u_lcn(ie)*dn
          obj%delta_v_lcn(ia,l,2) = obj%delta_v_lcn(ia,l,2) + obj%e_tb%u_lcn(ie)*dn
        end do
        obj%delta_v_lcn(ia,3,1) = obj%delta_v_lcn(ia,3,1) + obj%e_tb%u_lcn_d(ie)*dn_d
        obj%delta_v_lcn(ia,3,2) = obj%delta_v_lcn(ia,3,2) + obj%e_tb%u_lcn_d(ie)*dn_d
      end do
    end select
  end subroutine calculate_delta_v_lcn

  subroutine calculate_delta_v_pen(obj)
    ! INPUT
    class(hamiltonian_tb),intent(inout) :: obj
    ! LOCAL
    integer :: ia,l,ispin,jspin,sigma

    if (.not.allocated(obj%delta_v_pen)) allocate(obj%delta_v_pen(obj%a_tb%na,3,obj%a_tb%ns))
    obj%delta_v_pen = cmplx(0.0_rp,0.0_rp,kind=rp)

    if(obj%m_penalization/= 'none') then
      call obj%calculate_b_pen()

      select case(obj%a_tb%ns)
      case(2)
        do ia=1,obj%a_tb%na ! atom loop
          do l=1,3 ! subshell loop
            do ispin=1,2
              sigma = 3-2*ispin ! sigma = +-1
              !obj%delta_v_pen(ia,l,ispin) &
              ! = obj%delta_v_pen(ia,l,ispin) &
              ! - obj%a_tb%b_pen(ia,3)
              obj%delta_v_pen(ia,l,ispin) &
               = obj%delta_v_pen(ia,l,ispin) &
               - obj%a_tb%b_pen(ia,3)*sigma
            end do
          end do ! subshell loop
        end do ! atom loop
      case(4)
        do ia=1,obj%a_tb%na ! atom loop
          do l=1,3 ! subshell loop (s,p,d)
            do ispin=1,2
              do jspin=1,2
                obj%delta_v_pen(ia,l,obj%a_tb%iss2is(ispin,jspin))&
                 = obj%delta_v_pen(ia,l,obj%a_tb%iss2is(ispin,jspin))&
                 - (obj%a_tb%b_pen(ia,1)*sigma_x(ispin,jspin)&
                  + obj%a_tb%b_pen(ia,2)*sigma_y(ispin,jspin)&
                  + obj%a_tb%b_pen(ia,3)*sigma_z(ispin,jspin))
              end do
            end do
          end do ! subshell loop
        end do ! atom loop
      end select
    end if

  end subroutine calculate_delta_v_pen

  subroutine calculate_h_r(obj)
    ! INPUT
    class(hamiltonian_tb),intent(inout) :: obj
    ! LOCAL
    integer :: ia,ie,io
#if defined(DEBUG)
    write(output_unit,*) 'DEBUG == Entering calculate_h_r'
    call TBKOSTER_flush(output_unit)
#endif

    select case(obj%e_tb%tb_type)
    case('nrl')
      obj%en_intra = obj%a_tb%build_en_intra()
      obj%h_r = obj%a_tb%build_b_r()
      obj%h_r(:,0,:,:) = 0.0_rp

      do ia=1,obj%a_tb%na
        ie=obj%a_tb%ia2ie(ia)
        do io=1,obj%e_tb%no(ie)
          obj%h_r(ia,0,io,io) = obj%en_intra(ia,io)
        end do
      end do
    case('wan','mod')
      obj%h_r = obj%a_tb%build_b_r()
    end select

  end subroutine calculate_h_r
  
  subroutine calculate_s_r(obj)
    ! INPUT
    class(hamiltonian_tb),intent(inout) :: obj
    ! LOCAL
    integer :: ia,ie,io
    obj%s_r = obj%a_tb%build_s_r()
    obj%s_r(:,0,:,:) = 0.0_rp
    do ia=1,obj%a_tb%na
      ie=obj%a_tb%ia2ie(ia)
      do io=1,obj%e_tb%no(ie)
        obj%s_r(ia,0,io,io) = 1.0_rp
      end do
    end do
  end subroutine calculate_s_r

  !> Check the validity of the electron-electron interaction
  subroutine check_e_e_interaction(e_e_interaction)
    character(len=*),intent(in) :: e_e_interaction

    if(e_e_interaction /= 'stoner' .and. e_e_interaction /= 'ujb') then
      write(error_unit,*) 'element%check_e_e_interaction(): &
       &element%e_e_interaction must be one of: ''stoner'', ''ujb'''
      error stop
    end if
  end subroutine check_e_e_interaction

  !> Check the validity of the magnetic penalization
  subroutine check_m_penalization(m_penalization)
    character(len=*) :: m_penalization

    if(m_penalization /= 'none' &
     .and. m_penalization /= 'r' &
     .and. m_penalization /= 'r,theta' &
     .and. m_penalization /= 'r,theta,phi' &
     .and. m_penalization /= 'theta' &
     .and. m_penalization /= 'theta,phi' &
     .and. m_penalization /= 'phi') then
      write(error_unit,*) 'atom%check_m_penalization(): atom%m_penalization must be &
       &one of: ''none'', ''r'', ''r,theta'', ''r,theta,phi'', ''theta'', &
       &''theta,phi'', ''phi'''
      error stop
    end if
  end subroutine check_m_penalization

  subroutine initialize(obj)
    class(hamiltonian_tb),intent(inout) :: obj
    integer :: ih,ia,ie,io,is

    if(allocated(obj%iaos2ih)) deallocate(obj%iaos2ih)

    select case(obj%a_tb%ns)
    case(1,2)
      
      allocate(obj%iaos2ih(obj%a_tb%na,obj%e_tb%no_max,obj%a_tb%ns))

      ih = 0
      do ia=1,obj%a_tb%na
        ie = obj%a_tb%ia2ie(ia)
        do io=1,obj%e_tb%no(ie)
          ih = ih+1
          do is=1,obj%a_tb%ns
            obj%iaos2ih(ia,io,is) = ih
          end do
        end do
      end do
      obj%nh = ih
    case(4)
      
      allocate(obj%iaos2ih(obj%a_tb%na,obj%e_tb%no_max,2))

      ih = 0
      do ia=1,obj%a_tb%na
        ie = obj%a_tb%ia2ie(ia)
        do io=1,obj%e_tb%no(ie)
          do is=1,2
            ih = ih+1
            obj%iaos2ih(ia,io,is) = ih
          end do
        end do
      end do
      obj%nh = ih
    end select

    if(allocated(obj%en_intra)) deallocate(obj%en_intra)
    allocate(obj%en_intra(obj%a_tb%na,obj%e_tb%no_max))
    if(allocated(obj%h_r)) deallocate(obj%h_r)
    allocate(obj%h_r(obj%a_tb%na, 0:obj%a_tb%nn_max, obj%e_tb%no_max, obj%e_tb%no_max))
    if(allocated(obj%s_r)) deallocate(obj%s_r)
    allocate(obj%s_r(obj%a_tb%na, 0:obj%a_tb%nn_max, obj%e_tb%no_max, obj%e_tb%no_max))
    if(allocated(obj%delta_h_eei)) deallocate(obj%delta_h_eei)
    allocate(obj%delta_h_eei(obj%a_tb%na, obj%e_tb%no_max, obj%e_tb%no_max, obj%a_tb%ns))
    if(allocated(obj%delta_v_lcn)) deallocate(obj%delta_v_lcn)
    allocate(obj%delta_v_lcn(obj%a_tb%na,3,obj%a_tb%ns))
    if(allocated(obj%delta_v_pen)) deallocate(obj%delta_v_pen)
    allocate(obj%delta_v_pen(obj%a_tb%na,3,obj%a_tb%ns))

    obj%en_intra = 0.0_rp
    obj%h_r = 0.0_rp
    obj%s_r = 0.0_rp
    obj%delta_h_eei = cmplx(0.0_rp,0.0_rp,kind=rp)
    obj%delta_v_lcn = cmplx(0.0_rp,0.0_rp,kind=rp)
    obj%delta_v_pen = cmplx(0.0_rp,0.0_rp,kind=rp)
  end subroutine initialize

  !> Initialize electron-electron interaction
  subroutine initialize_eei(e_e_interaction)
    character(len=6),intent(out) :: e_e_interaction

    e_e_interaction = 'stoner'
  end subroutine initialize_eei

  !> Initialize magnetic penalization
  subroutine initialize_pen(m_penalization)
    character(len=11),intent(out) :: m_penalization

    m_penalization = 'none'
  end subroutine initialize_pen

  !> Read object in text format from file (default: 'in_hamiltonian_tb.txt')
  subroutine read_txt(obj,file)
    class(hamiltonian_tb),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    character(len=6) :: e_e_interaction
    character(len=11) :: m_penalization
    ! Namelist
    namelist /hamiltonian_tb/ e_e_interaction, m_penalization

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_hamiltonian_tb.txt'
    end if

    inquire(unit=10, opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'hamiltonian_tb%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'hamiltonian_tb%read_txt(): file ', file_rt, &
       ' not found'
      error stop
    end if

    call initialize_eei(e_e_interaction)
    call initialize_pen(m_penalization)

    read(10,nml=hamiltonian_tb)

    e_e_interaction = trim(lower(e_e_interaction))
    m_penalization = trim(lower(m_penalization))
    call check_e_e_interaction(e_e_interaction)
    call check_m_penalization(m_penalization)
    obj%e_e_interaction = e_e_interaction
    obj%m_penalization = m_penalization

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_hamiltonian_tb.txt')
  subroutine write_txt(obj,file,unit)
    class(hamiltonian_tb),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=len(obj%e_e_interaction)) :: e_e_interaction
    character(len=len(obj%m_penalization)) :: m_penalization
    ! Namelist
    namelist /hamiltonian_tb/ e_e_interaction, m_penalization

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_hamiltonian_tb.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    e_e_interaction = obj%e_e_interaction
    m_penalization = obj%m_penalization

    write(unit_rt,nml=hamiltonian_tb)
    call TBKOSTER_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_hamiltonian_tb.txt'), if tag (default: .true.) the namelist opening
  !> and closing tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(hamiltonian_tb),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Local variables
    integer :: ia, io, io1, io2, in, ip, is, isp, l

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_hamiltonian_tb.txt'
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
      write(unit_rt,'(a)') '&hamiltonian_tb'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('e_e_interaction')
        write(unit_rt,'(a)') ' e_e_interaction = ''' &
         // trim(obj%e_e_interaction) // ''''
      case('m_penalization')
        write(unit_rt,'(a)') ' m_penalization = ''' &
         // trim(obj%m_penalization) // ''''
      case('nh')
        write(unit_rt,'(a)') ' nh = ' // int2str(obj%nh)
      case('ns')
        write(unit_rt,'(a)') ' ns = ' // int2str(obj%a_tb%ns)
      case('iaos2ih')
        do ia=1,obj%a_tb%na
          do io=1,obj%e_tb%no(obj%a_tb%ia2ie(ia))
            do isp=1,obj%a_tb%nsp
              write(unit_rt,'(a)') ' iaos2ih(' // int2str(ia) // ',' &
               // int2str(io) // ',' // int2str(isp) // ') = ' &
               // int2str(obj%iaos2ih(ia,io,isp))
            end do
          end do
        end do
      case('en_intra')
        do ia=1,obj%a_tb%na
          do io=1,obj%e_tb%no_max
            write(unit_rt,'(a)') ' en_intra(' // int2str(ia) // ',' &
             // int2str(io)// ') = ' // real2str(obj%en_intra(ia,io))
          end do
        end do
      case('h_r')
        do ia=1,obj%a_tb%na
          do in=1,obj%a_tb%nn_max
            do io1=1,obj%e_tb%no_max
              do io2=1,obj%e_tb%no_max
                write(unit_rt,'(a)') ' h_r(' // int2str(ia) // ',' &
                 // int2str(in) // ',' // int2str(io1) // ',' // int2str(io2) &
                 // ') = ' // real2str(obj%h_r(ia,in,io1,io2))
              end do
            end do
          end do
        end do
      case('s_r')
        do ia=1,obj%a_tb%na
          do in=1,obj%a_tb%nn_max
            do io1=1,obj%e_tb%no_max
              do io2=1,obj%e_tb%no_max
                write(unit_rt,'(a)') ' s_r(' // int2str(ia) // ',' &
                 // int2str(in) // ',' // int2str(io1) // ',' // int2str(io2) &
                 // ') = ' // real2str(obj%s_r(ia,in,io1,io2))
              end do
            end do
          end do
        end do
      case('delta_h_eei')
        do ia=1,obj%a_tb%na
          do io1=1,obj%e_tb%no_max
            do io2=1,obj%e_tb%no_max
              do is=1,obj%a_tb%ns
                write(unit_rt,'(a)') ' delta_h_eei(' // int2str(ia) // ',' &
                 // int2str(io1) // ',' // int2str(io2) // ',' // int2str(is) &
                 // ') = ' // cmplx2str(obj%delta_h_eei(ia,io1,io2,is))
              end do
            end do
          end do
        end do
      case('delta_v_lcn')
        do ia=1,obj%a_tb%na
          do l=1,3
            do is=1,obj%a_tb%ns
              write(unit_rt,'(a)') ' delta_v_lcn(' // int2str(ia) // ',' &
               // int2str(l) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%delta_v_lcn(ia,l,is))
            end do
          end do
        end do
      case('delta_v_pen')
        do ia=1,obj%a_tb%na
          do l=1,3
            do is=1,obj%a_tb%ns
              write(unit_rt,'(a)') ' delta_v_pen(' // int2str(ia) // ',' &
               // int2str(l) // ',' // int2str(is) // ') = ' &
               // cmplx2str(obj%delta_v_pen(ia,l,is))
            end do
          end do
        end do
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
end module hamiltonian_tb_mod
