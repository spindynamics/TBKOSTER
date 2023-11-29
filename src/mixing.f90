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
!  mixing.f90
!  TBKOSTER
module mixing_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use atom_mod
  use charge_mod
  use element_mod
  use math_mod, only: i_unit,inverse,is_hermitian_r4,nm2rho,rho2nm
  use precision_mod, only: rp
  use string_mod, only: sl, int2str, lower, real2str, TBKOSTER_flush
  implicit none
  private

  !> Derived type properties for i/o methods
  character(len=sl),dimension(*),parameter :: property_list = &
   [character(len=sl) :: &
   'type', &
   'alpha', &
   'n_init', &
   'n_hist' &
   ]

  type,public :: mixing
    ! Elements
    class(element),pointer :: e
    ! Atoms
    class(atom),pointer :: a
    ! Charges
    class(charge),pointer :: q

    ! Mixing type ; options:
    !  'linear'
    !  'anderson' (not implemented)
    !  'broyden'
    character(len=7) :: type
    ! Linear mixing factor
    real(rp) :: alpha
    ! Linear mixing number of steps
    integer :: n_init
    ! Broyden mixing maximum number of saved iterations
    integer :: n_hist
    ! Broyden mixing only on the diagonal
    ! (e.g. TB+U: broyden sur les diagonal, linear sur les off-diagonal)
    logical :: broyden_diagonal

    ! Dimensions
    integer :: nals,naos,naoos
    ! Broyden mixing
    integer :: vlen
    real(rp),dimension(:,:),allocatable :: vector
    real(rp),dimension(:,:),allocatable :: u, vt
    real(rp),dimension(:)  ,allocatable :: f, df, vold, w

  contains
    ! Destructor
    final :: destructor
    ! Procedures
    procedure :: initialize_broyden
    procedure :: mix
    procedure :: mix_broyden
    procedure :: read_txt
    procedure :: write_txt
    procedure :: write_txt_formatted
  end type mixing

  ! Constructor
  interface mixing
    procedure :: constructor
  end interface mixing

contains
  function constructor(q) result(obj)
    class(charge),target,intent(in) :: q
    type(mixing) :: obj

    obj%e => q%e
    obj%a => q%a
    obj%q => q
  end function constructor

  subroutine destructor(obj)
    type(mixing) :: obj

    if(allocated(obj%vector)) deallocate(obj%vector)
    if(allocated(obj%u))      deallocate(obj%u)
    if(allocated(obj%vt))     deallocate(obj%vt)
    if(allocated(obj%f))      deallocate(obj%f)
    if(allocated(obj%df))     deallocate(obj%df)
    if(allocated(obj%vold))   deallocate(obj%vold)
    if(allocated(obj%w))      deallocate(obj%w)
  end subroutine destructor

  subroutine check_type(type)
    character(len=*),intent(in) :: type

    if(type /= 'linear' &
     .and. type /= 'anderson' &
     .and. type /= 'broyden') then
      write(error_unit,*) 'mixing%check_type(): mixing%type must be one of: &
       &''linear'', ''anderson'', ''broyden'''
      error stop
    end if
  end subroutine check_type

  subroutine initialize_broyden(obj)
    class(mixing),intent(inout) :: obj

    obj%nals = obj%a%na*3*obj%a%ns
    obj%naos = obj%a%na*obj%e%no_max*obj%a%ns
    obj%naoos = obj%a%na*obj%e%no_max*obj%e%no_max*obj%a%ns

    if(obj%broyden_diagonal) then
      obj%vlen = obj%nals + obj%naos
    else
      select case(obj%a%ns)
      case(1,2)
        obj%vlen = obj%nals + obj%naoos
      case(4)
        obj%vlen = obj%nals + 2*obj%naoos
      end select
    end if

    if(allocated(obj%vector)) deallocate(obj%vector)
    allocate(obj%vector(obj%vlen,2))
    if(allocated(obj%u)) deallocate(obj%u)
    allocate(obj%u(obj%vlen,obj%n_hist))
    if(allocated(obj%vt)) deallocate(obj%vt)
    allocate(obj%vt(obj%vlen,obj%n_hist))
    if(allocated(obj%f)) deallocate(obj%f)
    allocate(obj%f(obj%vlen))
    if(allocated(obj%df)) deallocate(obj%df)
    allocate(obj%df(obj%vlen))
    if(allocated(obj%vold)) deallocate(obj%vold)
    allocate(obj%vold(obj%vlen))
    if(allocated(obj%w)) deallocate(obj%w)
    allocate(obj%w(obj%n_hist))
  end subroutine initialize_broyden

  subroutine initialize_mixing(type,alpha,n_init,n_hist,broyden_diagonal)
    character(len=7),intent(out) :: type
    real(rp),intent(out) :: alpha
    integer,intent(out) :: n_init,n_hist
    logical,intent(out) :: broyden_diagonal

    type = 'broyden'
    alpha = 0.3_rp
    n_init = 1
    n_hist = 50
    !if(n_init==n_hist) type='linear'
    broyden_diagonal = .false.
  end subroutine initialize_mixing

  subroutine mix(obj,iter)
    ! INPUT
    class(mixing),intent(inout) :: obj
    integer,intent(in) :: iter
    ! LOCAL
    integer :: ia,il,io1,io2,is1,is2
    integer :: ispin,jspin
    integer :: nn
    ! Real and imaginary parts of rho
    complex(rp),dimension(2,2) :: rho_in_i,rho_out_i,rho_in_r,rho_out_r
    real(rp) :: m_cart_in_i(3),n_in_i,m_cart_out_i(3),n_out_i
    real(rp) :: m_cart_in_r(3),n_in_r,m_cart_out_r(3),n_out_r
    real(rp) :: rms

    select case(obj%a%ns)
    case(1,2)

      select case(obj%type)
      case('broyden')

        if(obj%broyden_diagonal) then

          do ia=1,obj%a%na
            do il=1,3
              do is1=1,obj%a%ns
                nn = (ia-1)*3*obj%a%ns + (il-1)*obj%a%ns + is1
                obj%vector(nn,1) = obj%q%q_mul_in(ia,il,is1-1)
                obj%vector(nn,2) = obj%q%q_mul_out(ia,il,is1-1)
              end do
            end do
          end do

          rms = sqrt(sum((obj%vector(:,1)-obj%vector(:,2))**2))

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do is2=1,obj%a%ns
                nn = (ia-1)*obj%e%no_max*obj%a%ns + (io1-1)*obj%a%ns &
                 + is2 + obj%nals
                obj%vector(nn,1) = real(obj%q%rho_net_in(ia,io1,io1,is2))
                obj%vector(nn,2) = real(obj%q%rho_net_out(ia,io1,io1,is2))
              end do
            end do
          end do

          rms = sqrt(sum((obj%vector(:,1)-obj%vector(:,2))**2))

          call obj%mix_broyden(iter,rms)

          do ia=1,obj%a%na
            do il=1,3
              do is1=1,obj%a%ns
                nn = (ia-1)*3*obj%a%ns + (il-1)*obj%a%ns + is1
                obj%q%q_mul_out(ia,il,is1-1) = obj%vector(nn,2)
              end do
            end do
          end do

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do is2=1,obj%a%ns
                nn = (ia-1)*obj%e%no_max*obj%a%ns + (io1-1)*obj%a%ns &
                 + is2 + obj%nals
                obj%q%rho_net_out(ia,io1,io1,is2) = obj%vector(nn,2)
              end do
            end do
          end do

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do io2=1,obj%e%no_max
                if(io1/=io2) then
                  do is2=1,obj%a%ns
                    obj%q%rho_net_out(ia,io1,io2,is2) &
                     = obj%alpha*obj%q%rho_net_out(ia,io1,io2,is2) &
                     + (1-obj%alpha)*obj%q%rho_net_in(ia,io1,io2,is2)
                  end do !end is2
                end if
              end do  !end io1
            end do  !end io2
          end do !end ia

        else !if(.not. obj%broyden_diagonal) then

          do ia=1,obj%a%na
            do il=1,3
              do is1=1,obj%a%ns
                nn = (ia-1)*3*obj%a%ns + (il-1)*obj%a%ns + is1
                obj%vector(nn,1) = obj%q%q_mul_in(ia,il,is1-1)
                obj%vector(nn,2) = obj%q%q_mul_out(ia,il,is1-1)
              end do
            end do
          end do

          rms = sqrt(sum((obj%vector(:,1)-obj%vector(:,2))**2))

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do io2=1,obj%e%no_max
                do is2=1,obj%a%ns
                  nn = (ia-1)*obj%e%no_max*obj%e%no_max*obj%a%ns &
                   + (io1-1)*obj%e%no_max*obj%a%ns &
                   + (io2-1)*obj%a%ns + is2 + obj%nals
                  obj%vector(nn,1) = real(obj%q%rho_net_in(ia,io1,io2,is2))
                  obj%vector(nn,2) = real(obj%q%rho_net_out(ia,io1,io2,is2))
                end do
              end do
            end do
          end do

          rms = sqrt(sum((obj%vector(:,1)-obj%vector(:,2))**2))

          call obj%mix_broyden(iter,rms)

          do ia=1,obj%a%na
            do il=1,3
              do is1=1,obj%a%ns
                nn = (ia-1)*3*obj%a%ns + (il-1)*obj%a%ns + is1
                obj%q%q_mul_out(ia,il,is1-1) = obj%vector(nn,2)
              end do
            end do
          end do

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do io2=1,obj%e%no_max
                do is2=1,obj%a%ns
                  nn = (ia-1)*obj%e%no_max*obj%e%no_max*obj%a%ns &
                   + (io1-1)*obj%e%no_max*obj%a%ns &
                   + (io2-1)*obj%a%ns + is2 + obj%nals
                  obj%q%rho_net_out(ia,io1,io2,is2) = obj%vector(nn,2)
                end do
              end do
            end do
          end do

        end if

      case('linear')
        obj%q%q_mul_out = obj%alpha*obj%q%q_mul_out &
         + (1-obj%alpha)*obj%q%q_mul_in
        obj%q%rho_net_out = obj%alpha*obj%q%rho_net_out &
         + (1-obj%alpha)*obj%q%rho_net_in
      end select

    case(4)

      if(.not. is_hermitian_r4(obj%q%rho_net_in)) then
        write(error_unit,'(a)') 'mixing%mix(): obj%q%rho_net_in not Hermitian'
      end if
      if(.not. is_hermitian_r4(obj%q%rho_net_out)) then
        write(error_unit,'(a)') 'mixing%mix(): obj%q%rho_net_out not Hermitian'
      end if

      select case(obj%type)
      case('broyden')

        if(obj%broyden_diagonal) then

          do ia=1,obj%a%na
            do il=1,3
              do is1=1,obj%a%ns
                nn = (ia-1)*3*obj%a%ns + (il-1)*obj%a%ns + is1
                obj%vector(nn,1) = obj%q%q_mul_in(ia,il,is1-1)
                obj%vector(nn,2) = obj%q%q_mul_out(ia,il,is1-1)
              end do   !end is1
            end do   !end il
          end do   !end ia

          rms = sqrt(sum((obj%vector(:,1)-obj%vector(:,2))**2))

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max

              do ispin=1,2
                do jspin=1,2
                  rho_in_r(ispin,jspin)  = obj%q%rho_net_in(ia,io1,io1,obj%a%iss2is(ispin,jspin))
                  rho_out_r(ispin,jspin) = obj%q%rho_net_out(ia,io1,io1,obj%a%iss2is(ispin,jspin))
                end do !end ispin
              end do !end jspin

              call rho2nm(rho_in_r,n_in_r,m_cart_in_r)
              call rho2nm(rho_out_r,n_out_r,m_cart_out_r)

              do is2=1,obj%a%ns
                nn = (ia-1)*obj%e%no_max*obj%a%ns + (io1-1)*obj%a%ns &
                 + is2 + obj%nals
                if(is2==1) then
                  obj%vector(nn,1) = n_in_r
                  obj%vector(nn,2) = n_out_r
                elseif(is2==2) then
                  obj%vector(nn,1) = m_cart_in_r(1)
                  obj%vector(nn,2) = m_cart_out_r(1)
                elseif(is2==3) then
                  obj%vector(nn,1) = m_cart_in_r(2)
                  obj%vector(nn,2) = m_cart_out_r(2)
                elseif(is2==4) then
                  obj%vector(nn,1) = m_cart_in_r(3)
                  obj%vector(nn,2) = m_cart_out_r(3)
                end if

              end do !end is2
            end do !end io1
          end do !end ia

          rms = sqrt(sum((obj%vector(:,1)-obj%vector(:,2))**2))

          call obj%mix_broyden(iter,rms)

          do ia=1,obj%a%na
            do il=1,3
              do is1=1,obj%a%ns
                nn = (ia-1)*3*obj%a%ns + (il-1)*obj%a%ns + is1
                obj%q%q_mul_out(ia,il,is1-1) = obj%vector(nn,2)
              end do !end il
            end do !end is1
          end do !end ia

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do is2=1,obj%a%ns
                nn = (ia-1)*obj%e%no_max*obj%a%ns + (io1-1)*obj%a%ns &
                 + is2 + obj%nals

                if(is2==1) then
                  n_out_r = obj%vector(nn,2)
                elseif(is2==2) then
                  m_cart_out_r(1) = obj%vector(nn,2)
                elseif(is2==3) then
                  m_cart_out_r(2) = obj%vector(nn,2)
                elseif(is2==4) then
                  m_cart_out_r(3) = obj%vector(nn,2)
                end if
              end do  !end is2

              call nm2rho(n_out_r,m_cart_out_r,rho_out_r)

              do ispin=1,2
                do jspin=1,2
                  obj%q%rho_net_out(ia,io1,io1,obj%a%iss2is(ispin,jspin)) &
                   = rho_out_r(ispin,jspin)
                end do   !end ispin
              end do   !end jspin

            end do  !end io1
          end do  !end ia

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do io2=1,obj%e%no_max
                if(io1/=io2) then
                  do ispin=1,2
                    do jspin=1,2
                      obj%q%rho_net_out(ia,io1,io2,obj%a%iss2is(ispin,jspin)) &
                       = obj%alpha*obj%q%rho_net_out(ia,io1,io2,obj%a%iss2is(ispin,jspin)) &
                       + (1-obj%alpha)*obj%q%rho_net_in(ia,io1,io2,obj%a%iss2is(ispin,jspin))
                    end do !end ispin
                  end do !end jspin
                end if
              end do  !end io1
            end do  !end io2
          end do !end ia

        else !if(.not. obj%broyden_diagonal) then

          do ia=1,obj%a%na
            do il=1,3
              do is1=1,obj%a%ns
                nn = (ia-1)*3*obj%a%ns + (il-1)*obj%a%ns + is1
                obj%vector(nn,1) = obj%q%q_mul_in(ia,il,is1-1)
                obj%vector(nn,2) = obj%q%q_mul_out(ia,il,is1-1)
              end do !end is1
            end do !end il
          end do !end ia

          rms = sqrt(sum((obj%vector(:,1)-obj%vector(:,2))**2))

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do io2=1,obj%e%no_max

                do ispin=1,2
                  do jspin=1,2
                    rho_in_r(ispin,jspin) = 0.5_rp &
                     * (obj%q%rho_net_in(ia,io1,io2,obj%a%iss2is(ispin,jspin)) &
                     + obj%q%rho_net_in(ia,io2,io1,obj%a%iss2is(ispin,jspin)))
                    rho_out_r(ispin,jspin) = 0.5_rp &
                     * (obj%q%rho_net_out(ia,io1,io2,obj%a%iss2is(ispin,jspin)) &
                     + obj%q%rho_net_out(ia,io2,io1,obj%a%iss2is(ispin,jspin)))
                    rho_in_i(ispin,jspin) = -i_unit*0.5_rp &
                     * (obj%q%rho_net_in(ia,io1,io2,obj%a%iss2is(ispin,jspin)) &
                     - obj%q%rho_net_in(ia,io2,io1,obj%a%iss2is(ispin,jspin)))
                    rho_out_i(ispin,jspin) = -i_unit*0.5_rp &
                     *(obj%q%rho_net_out(ia,io1,io2,obj%a%iss2is(ispin,jspin)) &
                     - obj%q%rho_net_out(ia,io2,io1,obj%a%iss2is(ispin,jspin)))
                  end do !end ispin
                end do !end jspin

                call rho2nm(rho_in_i,n_in_i,m_cart_in_i)
                call rho2nm(rho_out_i,n_out_i,m_cart_out_i)
                call rho2nm(rho_in_r,n_in_r,m_cart_in_r)
                call rho2nm(rho_out_r,n_out_r,m_cart_out_r)

                do is2=1,obj%a%ns
                  nn = (ia-1)*obj%e%no_max*obj%e%no_max*2*obj%a%ns &
                   + (io1-1)*obj%e%no_max*2*obj%a%ns &
                   + (io2-1)*2*obj%a%ns + is2 + obj%nals
                  if(is2==1) then
                    obj%vector(nn,1) = n_in_r
                    obj%vector(nn,2) = n_out_r
                  elseif(is2==2) then
                    obj%vector(nn,1) = m_cart_in_r(1)
                    obj%vector(nn,2) = m_cart_out_r(1)
                  elseif(is2==3) then
                    obj%vector(nn,1) = m_cart_in_r(2)
                    obj%vector(nn,2) = m_cart_out_r(2)
                  elseif(is2==4) then
                    obj%vector(nn,1) = m_cart_in_r(3)
                    obj%vector(nn,2) = m_cart_out_r(3)
                  end if

                end do  !end is2
                do is2=obj%a%ns+1,2*obj%a%ns
                  nn = (ia-1)*obj%e%no_max*obj%e%no_max*2*obj%a%ns &
                   + (io1-1)*obj%e%no_max*2*obj%a%ns &
                   + (io2-1)*2*obj%a%ns + is2 + obj%nals
                  if(is2==5) then
                    obj%vector(nn,1) = n_in_i
                    obj%vector(nn,2) = n_out_i
                  elseif(is2==6) then
                    obj%vector(nn,1) = m_cart_in_i(1)
                    obj%vector(nn,2) = m_cart_out_i(1)
                  elseif(is2==7) then
                    obj%vector(nn,1) = m_cart_in_i(2)
                    obj%vector(nn,2) = m_cart_out_i(2)
                  elseif(is2==8) then
                    obj%vector(nn,1) = m_cart_in_i(3)
                    obj%vector(nn,2) = m_cart_out_i(3)
                  end if

                end do !end is2

              end do !end io1
            end do !end io2
          end do !end ia

          rms = sqrt(sum((obj%vector(:,1)-obj%vector(:,2))**2))

          call obj%mix_broyden(iter,rms)

          do ia=1,obj%a%na
            do il=1,3
              do is1=1,obj%a%ns
                nn = (ia-1)*3*obj%a%ns + (il-1)*obj%a%ns + is1
                obj%q%q_mul_out(ia,il,is1-1) = obj%vector(nn,2)
              end do !end is1
            end do !end il
          end do !end ia

          do ia=1,obj%a%na
            do io1=1,obj%e%no_max
              do io2=1,obj%e%no_max

                do is2=1,obj%a%ns
                  nn = (ia-1)*obj%e%no_max*obj%e%no_max*2*obj%a%ns &
                   + (io1-1)*obj%e%no_max*2*obj%a%ns &
                   + (io2-1)*2*obj%a%ns + is2 + obj%nals

                  if(is2==1) then
                    n_out_r = obj%vector(nn,2)
                  elseif(is2==2) then
                    m_cart_out_r(1) = obj%vector(nn,2)
                  elseif(is2==3) then
                    m_cart_out_r(2) = obj%vector(nn,2)
                  elseif(is2==4) then
                    m_cart_out_r(3) = obj%vector(nn,2)
                  end if
                end do  !end is2

                do is2=obj%a%ns+1,2*obj%a%ns
                  nn = (ia-1)*obj%e%no_max*obj%e%no_max*2*obj%a%ns &
                   + (io1-1)*obj%e%no_max*2*obj%a%ns &
                   + (io2-1)*2*obj%a%ns + is2 + obj%nals

                  if(is2==5) then
                    n_out_i = obj%vector(nn,2)
                  elseif(is2==6) then
                    m_cart_out_i(1) = obj%vector(nn,2)
                  elseif(is2==7) then
                    m_cart_out_i(2) = obj%vector(nn,2)
                  elseif(is2==8) then
                    m_cart_out_i(3) = obj%vector(nn,2)
                  end if
                end do  !end is2

                call nm2rho(n_out_i,m_cart_out_i,rho_out_i)
                call nm2rho(n_out_r,m_cart_out_r,rho_out_r)

                do ispin=1,2
                  do jspin=1,2
                    obj%q%rho_net_out(ia,io1,io2,obj%a%iss2is(ispin,jspin)) &
                     = rho_out_r(ispin,jspin) + i_unit*rho_out_i(ispin,jspin)
                  end do !end ispin
                end do !end jspin

              end do !end io1
            end do !end io2
          end do !end ia

        end if

      case('linear')
        obj%q%q_mul_out = obj%alpha*obj%q%q_mul_out + (1-obj%alpha)*obj%q%q_mul_in
        obj%q%rho_net_out = obj%alpha*obj%q%rho_net_out + (1-obj%alpha)*obj%q%rho_net_in
      end select

      if(.not. is_hermitian_r4(obj%q%rho_net_in)) then
        write(error_unit,'(a)') 'mixing%mix(): obj%q%rho_net_in not Hermitian'
      end if
      if(.not. is_hermitian_r4(obj%q%rho_net_out)) then
        write(error_unit,'(a)') 'mixing%mix(): obj%q%rho_net_out not Hermitian'
      end if
    end select
  end subroutine mix

  subroutine mix_broyden(obj,iter,rms)
    !===========================================================================
    ! History: Original code written by D.D. Johnson (see PRB 38, 12807)
    ! Note: There are a few typos in that paper but the code is working!
    ! Rewritten by W. A. Shelton for LSMS code 6/21/94
    !   this version is easy to read (no goto!!!! more comments ...)
    !   and is setup for MPP machines (not tested)
    ! Rewritten by T. C. Schulthess, ORNL, March 97
    !   this version should work for any code (see comments below)
    ! Bug fixes: TCS, 8/5/97 see comments below
    !===========================================================================
    !
    !     further comments on how to use this subroutine:
    !     (if nothing useful stands here I had no time yet to update these
    !     comments, please consult usage in lkkr code version 0.6.3 or later,
    !     or call Thomas Schulthess (423) 5768067)
    !
    !     obj%vector(r,i) -> i=1: old vector (input), scratch (ouput)
    !                     -> i=2: new vector (input), mixed vector (output)
    !     obj%vlen    -> length of vectors
    !     obj%alpha   -> linear mixing factor
    !     rms         -> RMS difference between old and new vector
    !     iter        -> iteration number (if 1, linear mixing, broyden reset)
    !     obj%n_hist  -> number of iterations that are used for mixing (<=obj%n_hist)
    !                 -> maximum number of iterations that can be saved
    !     obj%u, obj%vt, obj%f, obj%df, obj%vold, and obj%w are working arrays
    !     a, b, d, cm, and ipiv are working arrays that need not be saved
    !
    !     See declaration for exact dimentsions and types
    !
    !     There are two options for matrix inversions, a Gaussian
    !     elimination routine called invert1 and calls to lapack routines
    !     with pivoting (see comments "using invert1" and "using lapack").
    !     Obviously only one can be used, comment out the other one.
    !
    !     When using this subroutine in a parallel code in which only parts
    !     of the vectors are known on every node, make sure that the calls
    !     to gldsum (global sum) are correct (LKKR and LSMS codes have
    !     different calls).
    !
    !     In a serial code, either comment out the calls to glbsum or
    !     provide a dummy subroutine
    !
    !===========================================================================
    ! INPUT
    class(mixing),intent(inout) :: obj
    integer :: iter
    real(rp) :: rms
    ! LOCAL
    real(rp),dimension(:,:), allocatable :: a, b, d
    real(rp),dimension(:), allocatable :: cm
    integer,dimension(:), allocatable :: ipiv
    !
    integer :: i,j,k,info
    real(rp) :: fac1,fac2,fnorm,dfnorm,w0,aij,gmi,cmj,wtmp
    !
    integer :: last_iter,nn

    last_iter = iter-1
    allocate(a(obj%n_hist,obj%n_hist),b(obj%n_hist,obj%n_hist),d(obj%n_hist,obj%n_hist))
    allocate(cm(obj%n_hist))
    allocate(ipiv(obj%n_hist))

    if(iter<=obj%n_init) then
      !=========================================================================
      ! first obj%n_init iterations: perform linear mixing, load f and vold, set
      !                  different pointers and variables
      !=========================================================================
      do k = 1,obj%vlen
        obj%f(k) = obj%vector(k,2) - obj%vector(k,1)
        obj%vold(k) = obj%vector(k,1)
      end do

      do k = 1,obj%vlen          ! this is the linear mixing
        obj%vector(k,2) = obj%vector(k,1) + obj%alpha * obj%f(k)
      end do
      !=========================================================================
    else
      !=========================================================================
      !            iter > obj%n_init: perform non-linear mixing
      !=========================================================================
      if(iter > obj%n_hist) then ! set current length of broyden cycle
        nn = obj%n_hist
      else
        nn = last_iter
      end if

      w0 = 1E-2_rp              ! set weighting factor for the zeroth iteration

      !---- find: f[i] := vector(2)[i] - vector(1)[i]

      do k = 1,obj%vlen
        obj%df(k) = obj%vector(k,2) - obj%vector(k,1) - obj%f(k)
      end do
      do k = 1,obj%vlen
        obj%f(k) = obj%vector(k,2) - obj%vector(k,1)
      end do

      !---- find: fnorm  := |f|

      dfnorm = 0.0_rp
      fnorm = 0.0_rp
      do k = 1,obj%vlen
        dfnorm = dfnorm + obj%df(k)*obj%df(k)
        fnorm  = fnorm  + obj%f(k)*obj%f(k)
      end do

      dfnorm = sqrt( dfnorm )
      fnorm  = sqrt( fnorm )

      !---- set: vector(2) := alpha*df/|df| + (vector(1) - vold)/|df|

      fac2 = 1.0_rp/dfnorm
      fac1 = obj%alpha*fac2
      do k = 1,obj%vlen
        obj%vector(k,2) = fac1*obj%df(k) + fac2*(obj%vector(k,1) - obj%vold(k))
        obj%vold(k) = obj%vector(k,1)
        obj%vector(k,1) = fac2*obj%df(k)
      end do

      !---- store vector(1) and vector(2) in the stacks u and vt restpectively

      !=========================================================================
      if(iter-1<=obj%n_hist) then
        !=======================================================================
        !   Load the first obj%n_hist iterations in increasing iteration count
        !=======================================================================
        do i=1,obj%vlen
          obj%u(i,iter-1) = obj%vector(i,2)
        end do
        !=======================================================================
        do i=1,obj%vlen
          obj%vt(i,iter-1) = obj%vector(i,1)
        end do
        !=======================================================================
      else
        !=======================================================================
        !     Re-load so that the ordering is in increasing iteration count
        !=======================================================================
        do j=1,obj%n_hist-1
          !          write(output_unit,'('' IN BROY_SAV: j,j+1 '',2is2)') j,j+1
          do i=1,obj%vlen
            obj%u(i,j) = obj%u(i,j+1)
          end do
          !=====================================================================
          do i=1,obj%vlen
            obj%vt(i,j) = obj%vt(i,j+1)
          end do
          !=====================================================================
        end do
        !=======================================================================
        !     Load current charge densities in the last storage location
        !=======================================================================
        do i=1,obj%vlen
          obj%u(i,obj%n_hist) = obj%vector(i,2)
        end do
        !=======================================================================
        do i=1,obj%vlen
          obj%vt(i,obj%n_hist) = obj%vector(i,1)
        end do
        !=======================================================================
      end if
      !=========================================================================

      !---- calculate coefficient matrices, a(i,j), and sum cm(i) for corrections:

      do j=1,nn - 1          ! off diagonal elements of a(i,j)
        do i = j+1,nn
          aij = 0.0_rp
          do k = 1,obj%vlen
            aij = aij + obj%vt(k,j)*obj%vt(k,i)
          end do
          a(i,j) = aij
          a(j,i) = aij
        end do
      end do

      do i = 1,nn             ! diagonal elements a(i,i) and cm(i)
        aij = 0.0_rp
        cmj = 0.0_rp
        do k=1,obj%vlen
          cmj = cmj + obj%vt(k,i)*obj%f(k)
          aij = aij + obj%vt(k,i)*obj%vt(k,i)
        end do
        a(i,i) = aij
        cm(i) = cmj
      end do

      !---- shift down weights in stack

      if(iter-1 > obj%n_hist)then
        do i=1,obj%n_hist-1
          obj%w(i) = obj%w(i+1)
        end do
      end if
      wtmp = 0.0_rp
      if(rms > 1.0E-09_rp) wtmp = 2.0_rp*sqrt(0.01_rp/rms)
      if(wtmp < 1.0_rp)    wtmp = 1.0_rp
      if(iter > obj%n_hist) then
        obj%w(obj%n_hist) = wtmp
      else
        obj%w(last_iter) = wtmp
      end if

      !---- now calculate the b-matrix:
      !>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> uses lapack
      do i=1,nn
        do j=1,nn
          b(j,i) = a(j,i)*obj%w(j)*obj%w(i)
        end do
        b(i,i) = w0**2 + a(i,i)*obj%w(i)*obj%w(i)
      end do
      call dgetrf(nn,nn,b,obj%n_hist,ipiv,info)
      call dgetri(nn,b,obj%n_hist,ipiv,d,nn,info)
      !b = inverse(b)
      !  write(output_unit,*) ' optimum lwork', d(1,1)
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< uses lapack
      !---- mix vectors:

      do k=1,obj%vlen
        obj%vector(k,2) = obj%vold(k) + obj%alpha * obj%f(k)
      end do
      do i=1,nn
        gmi = 0.0_rp
        do j=1,nn
          gmi = gmi + cm(j)*b(j,i)*obj%w(j)
        end do
        do k=1,obj%vlen
          obj%vector(k,2) = obj%vector(k,2) - gmi*obj%u(k,i)*obj%w(i)
        end do
      end do

    end if
    deallocate(ipiv,cm,d,b,a)
  end subroutine mix_broyden

  !> Read object in text format from file (default: 'in_scf.txt')
  subroutine read_txt(obj,file)
    class(mixing),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt
    integer :: iostatus
    logical :: isopen
    ! Namelist variables
    character(len=7) :: type
    real(rp) :: alpha
    integer :: n_init, n_hist
    logical :: broyden_diagonal
    ! Namelist
    namelist /mixing/ type, alpha, n_init, n_hist

    if(present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_mixing.txt'
    end if

    inquire(unit=10,opened=isopen)
    if (isopen) then
      write(error_unit,'(a)') 'mixing%read_txt() : Unit 10 is already open'
      error stop
    else
      open(unit=10,file=file_rt,action='read',iostat=iostatus,status='old')
    end if
    if(iostatus /= 0) then
      write(error_unit,*) 'mixing%read_txt(): file ', file_rt, ' not found'
      error stop
    end if

    call initialize_mixing(type,alpha,n_init,n_hist,broyden_diagonal)
    read(10,nml=mixing)
    type = lower(type)
    call check_type(trim(type))

    obj%type = type
    obj%alpha = alpha
    obj%n_init = n_init
    obj%n_hist  = n_hist
    obj%broyden_diagonal = broyden_diagonal

    call obj%initialize_broyden()

    close(unit=10)
    !deallocate(file_rt)
  end subroutine read_txt

  !> Write object in text format to unit (default: 10), if it's a file
  !> its name is set to file (default: 'out_mixing.txt')
  subroutine write_txt(obj,file,unit)
    class(mixing),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Namelist variables
    character(len=7) :: type
    real(rp) :: alpha
    integer :: n_init, n_hist
    ! Namelist
    namelist /mixing/ type, alpha, n_init, n_hist

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_mixing.txt'
    end if
    if(present(unit)) then
      unit_rt = unit
    else
      unit_rt = 10
    end if

    if(.not. present(unit)) then
      open(unit=unit_rt,file=file_rt,action='write')
    end if

    type = obj%type
    alpha = obj%alpha
    n_init = obj%n_init
    n_hist = obj%n_hist

    write(unit_rt,nml=mixing)
    call TBKOSTER_flush(unit_rt)

    if(.not. present(unit)) then
      close(unit_rt)
    end if
    !deallocate(file_rt)
  end subroutine write_txt

  !> Write property (default: property_list) in text format to unit
  !> (default: 10), if it's a file its name is set to file (default:
  !> 'out_mixing.txt'), if tag (default: .true.) the namelist opening and
  !> closing tags are written
  subroutine write_txt_formatted(obj,file,property,tag,unit)
    class(mixing),intent(in) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable         :: file_rt
    character(len=*),dimension(:),intent(in),optional :: property
    character(len=:),dimension(:),allocatable         :: property_rt
    logical,intent(in),optional :: tag
    logical                     :: tag_rt
    integer,intent(in),optional :: unit
    integer                     :: unit_rt
    ! Local variables
    integer :: ip

    if(present(file)) then
      file_rt = file
    else
      file_rt = 'out_mixing.txt'
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
      write(unit_rt,'(a)') '&mixing'
    end if

    do ip=1,size(property_rt)
      select case(lower(trim(property_rt(ip))))
      case('type')
        write(unit_rt,'(a)') ' type = ''' // trim(obj%type) // ''''
      case('alpha')
        write(unit_rt,'(a)') ' alpha = ' // real2str(obj%alpha)
      case('n_init')
        write(unit_rt,'(a)') ' n_init = ' // int2str(obj%n_init)
      case('n_hist')
        write(unit_rt,'(a)') ' n_hist = ' // int2str(obj%n_hist)
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
end module mixing_mod
