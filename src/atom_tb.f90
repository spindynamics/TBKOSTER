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
!  atom_tb.f90
!  TBKOSTER
module atom_tb_mod
  use, intrinsic :: iso_fortran_env, only: error_unit, output_unit
  use constant_mod, only :e_ry
  use atom_mod
  use element_tb_mod
  use lattice_mod
  use math_mod, only: sqrt_three, fermi_function, one_third, two_third, fermi_function_derivative
  use precision_mod, only: rp
  use string_mod, only : TBKOSTER_flush
  implicit none
  private

  type,public,extends(atom) :: atom_tb
    !> Elements TB
    class(element_tb),pointer :: e_tb

  contains
    ! Procedures
    procedure :: build_b_r
    procedure :: build_d_b_r
    procedure :: build_en_intra
    procedure :: build_d_en_intra
    procedure :: build_s_r
    procedure :: build_d_s_r
    procedure :: read_txt
  end type atom_tb

  ! Constructor
  interface atom_tb
    procedure :: constructor
  end interface atom_tb

contains
  function constructor(e_tb,l_r,l_k) result(obj)
    class(element_tb),target,intent(in) :: e_tb
    class(lattice),target,intent(in) :: l_r
    class(lattice),target,intent(in) :: l_k
    type(atom_tb) :: obj
    
    ! Parent type constructor
    obj%atom = atom(e_tb,l_r,l_k)

    ! Derived type constructor
    obj%e_tb => e_tb
  end function constructor

  function build_b_r(obj) result(B)
    ! INPUT
    class(atom_tb),intent(in) :: obj
    ! OUTPUT
    real(rp), dimension(:,:,:,:), allocatable :: B
    ! LOCAL
    real(rp) :: RR(3),R,f_cut
    real(rp) :: temp(9,9)
    real(rp) :: Overlap(10),I_Overlap(10) ! SS,SP,PP(2),SD,PD(2),DD(3)
    integer, dimension(:), allocatable ::  weight   
    integer  :: ia1,ia2,in,ie1,ie2,io1,io2,lbeta,step,i1,i2,i3,dummy1,dummy2,icase,ncase,norb
    real(rp) :: Btemp
    logical :: file_existence, isopen
    write(output_unit,*) 'DEBUG == Entering atom_tb & build_b_r'
    call TBKOSTER_flush(output_unit)

    if (.not.allocated(B)) allocate(B(obj%na,0:obj%nn_max,obj%e%no_max,obj%e%no_max))
    B = 0.0_rp
    select case(obj%e_tb%tb_type)
      case('nrl')
        do ia1=1,obj%na
          ie1 = obj%ia2ie(ia1)
          do in=1,obj%nn(ia1)
            ia2 = obj%ian2ia(ia1,in)
            ie2 = obj%ia2ie(ia2)

            RR(:) = obj%rn(ia1,in,:)
            R = norm2(RR)

            select case(obj%e_tb%nrl_type(ie1))
              case('old','new','cyr')
                do lbeta=1,10
                  step = 4*(lbeta-1)
                  Overlap(lbeta)=integral_parametrization_nrl(&
                  obj%e_tb%p(ie1,18+step),obj%e_tb%p(ie1,19+step),&
                  obj%e_tb%p(ie1,20+step),obj%e_tb%p(ie1,21+step),R)
                end do
              case('pow')
              do lbeta=1,10
                step = 4*(lbeta-1)
                Overlap(lbeta)=integral_parametrization_pow(&
                obj%e_tb%p(ie1,18+step),obj%e_tb%p(ie1,19+step),&
                obj%e_tb%p(ie1,21+step),R)
              end do
            end select

            f_cut = fermi_function(R-obj%e_tb%r_0(ie1), 1/obj%e_tb%r_l(ie1))
            Overlap(:)=f_cut*Overlap(:)

            if(ie1/=ie2) then
              
              select case(obj%e_tb%nrl_type(ie2))
                case('old','new','cyr')
                  do lbeta=1,10
                    step = 4*(lbeta-1)
                    I_Overlap(lbeta)=integral_parametrization_nrl(&
                    obj%e_tb%p(ie2,18+step),obj%e_tb%p(ie2,19+step),&
                    obj%e_tb%p(ie2,20+step),obj%e_tb%p(ie2,21+step),R)
                  end do
                case('pow')
                do lbeta=1,10
                  step = 4*(lbeta-1)
                  I_Overlap(lbeta)=integral_parametrization_pow(&
                  obj%e_tb%p(ie2,18+step),obj%e_tb%p(ie2,19+step),&
                  obj%e_tb%p(ie2,21+step),R)
                end do
              end select

              f_cut = fermi_function(R-obj%e_tb%r_0(ie2), 1/obj%e_tb%r_l(ie2))
              I_Overlap(:)=f_cut*I_Overlap(:)

              do lbeta=1,10
                Overlap(lbeta)=0.5_rp*(Overlap(lbeta)+I_Overlap(lbeta))
              end do
            end if
      
            temp = slater_koster(RR(:)/R,Overlap)

            do io1=1,obj%e%no(ie1)
              do io2=1,obj%e%no(ie2)
                B(ia1,in,io1,io2) = temp(obj%e%o(ie1,io1),obj%e%o(ie2,io2))
              end do
            end do
          end do
        end do
      case('wan')
        inquire(file='hr.dat',exist=file_existence)
        if (.not.file_existence) then
          write(error_unit,*) 'file hr.dat not present'
          error stop
        endif
        inquire(unit=10,opened=isopen)
        if (isopen) then
          write(error_unit,'(a)') 'build_b_r : Unit 10 is already open'
          error stop
        endif
        open(unit=10,file='hr.dat',action='read')
        read(10,*)
        read(10,*) norb
        read(10,*) ncase
        if (.not.allocated(weight)) allocate(weight(ncase))
        read(10,*) (weight(icase),icase=1,ncase)

        do icase=1,ncase
          do ia2=1,obj%na 
            ie2 = obj%ia2ie(ia2)
            do io2=1,obj%e%no(ie2)
                do ia1=1,obj%na
                  ie1 = obj%ia2ie(ia1)
                  do io1=1,obj%e%no(ie1)
                    read(10,*) i1,i2,i3,dummy1,dummy2,Btemp
                    in=obj%iapbc2in(ia1,ia2,i1+obj%pbc(1)+1,i2+obj%pbc(2)+1,i3+obj%pbc(3)+1)
                    ! Unit conversion from eV atomic units to Hartree atomic units       
                    B(ia1,in,io1,io2)=Btemp/weight(icase)*0.5_rp/e_ry
                  end do
                end do
            end do 
          end do  
        end do  

        if (allocated(weight)) deallocate(weight)
        close(unit=10)

        !   do ip1=1,2*obj%pbc(1)+1
        !    do ip2=1,2*obj%pbc(2)+1
        !     do ip3=1,2*obj%pbc(3)+1
        !        do ia2=1,obj%na 
        !           ie2 = obj%ia2ie(ia2)
        !           do io2=1,obj%e%no(ie2)
        !            do ia1=1,obj%na
        !               ie1 = obj%ia2ie(ia1)
        !               do io1=1,obj%e%no(ie1)
        !                read(10,*) i1,i2,i3,imat1,imat2,Btemp
        !                in=obj%iapbc2in(ia1,ia2,ip1,ip2,ip3)
        ! Unit conversion from eV atomic units to Hartree atomic units             
        !                B(ia1,in,io1,io2)=Btemp*0.5_rp/e_ry
        !               end do 
        !             end do
        !           end do
        !         end do
        !      end do
        !     end do
        !   end do
      
      case('mod')
        inquire(file='mod.dat',exist=file_existence)
        if (.not.file_existence) then
          write(error_unit,*) 'file mod.dat not present'
          error stop
        endif
        inquire(unit=10,opened=isopen)
        if (isopen) then
          write(error_unit,'(a)') 'build_b_r : Unit 10 is already open'
          error stop
        endif
        open(unit=10,file='mod.dat',action='read')
        read(10,*) ncase
        do icase=1,ncase
          read(10,*) i1,i2,i3,ia1,io1,ia2,io2,Btemp
          in=obj%iapbc2in(ia1,ia2,i1+obj%pbc(1)+1,i2+obj%pbc(2)+1,i3+obj%pbc(3)+1)
          ! Unit conversion from eV atomic units to Hartree atomic units             
          B(ia1,in,io1,io2)=Btemp*0.5_rp/e_ry
        end do
        close(unit=10)
    end select
    write(output_unit,*) "DEBUG B(1,1,1,1) = ", B(1,1,1,1)
    write(output_unit,*) 'DEBUG == Exiting atom_tb & build_b_r'
    call TBKOSTER_flush(output_unit)
  end function build_b_r

  ! Routine to calculate the derivative of the hopping matrix (d_B) 
  function build_d_b_r(obj) result(d_B)
    ! INPUT
    class(atom_tb),intent(in) :: obj
    ! OUTPUT
    real(rp), dimension(:,:,:,:,:), allocatable :: d_B
    ! LOCAL
    real(rp) :: AA(3),RR(3),R,f_cut,d_f_cut
    real(rp) :: d_Btemp(9,9,3) 
    real(rp) :: Overlap(10),d_Overlap(10,3),d_Overlaptemp(10) !SS,SP,PP(2),SD,PD(2),DD(3)
    real(rp) :: I_Overlap(10),d_I_Overlap(10,3),d_I_Overlaptemp(10) ! SS,SP,PP(2),SD,PD(2),DD(3)

    integer :: ia1,ia2,in,ie1,ie2,io1,io2,lbeta,ix,step

    if (.not.allocated(d_B)) allocate(d_B(obj%na,0:obj%nn_max,obj%e%no_max,obj%e%no_max,3))

    d_Btemp = 0.0_rp ; Overlap = 0.0_rp ; d_Overlap = 0.0_rp ; d_Overlaptemp = 0.0_rp
    d_B = 0.0_rp ; I_Overlap = 0.0_rp ; d_I_Overlap = 0.0_rp ; d_I_Overlaptemp = 0.0_rp

    ! For the hopping matrix elements
    do ia1=1,obj%na
      ie1 = obj%ia2ie(ia1)
      do in=1,obj%nn(ia1)
        ia2 = obj%ian2ia(ia1,in)
        ie2 = obj%ia2ie(ia2)

        RR(:) = obj%rn(ia1,in,:)
        R = norm2(RR)
        aa(:) = RR(:)/R

        f_cut = fermi_function(R-obj%e_tb%r_0(ie1), 1/obj%e_tb%r_l(ie1))
        d_f_cut = fermi_function_derivative(R-obj%e_tb%r_0(ie1),1/obj%e_tb%r_l(ie1))

        !Getting the hopping NRL parameters in the source file and computing the
        !value of the function. For the derivative, storing part of if on
        !d_Overlaptemp to later compute it at d_Overlap
        select case(obj%e_tb%nrl_type(ie1))
        case('old','new','cyr')
              
          do lbeta=1,10
            step = 4*(lbeta-1)

            !Calculating f(r)
            Overlap(lbeta)=integral_parametrization_nrl(&
            obj%e_tb%p(ie1,18+step),obj%e_tb%p(ie1,19+step),&
            obj%e_tb%p(ie1,20+step),obj%e_tb%p(ie1,21+step),R)

            !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r))
            d_Overlaptemp(lbeta)=d_integral_parametrization_nrl(&
            obj%e_tb%p(ie1,18+step),obj%e_tb%p(ie1,19+step),&
            obj%e_tb%p(ie1,20+step),obj%e_tb%p(ie1,21+step),R)*f_cut+&
            Overlap(lbeta)*d_f_cut
          end do
          
        case('pow')
                  
          do lbeta=1,10
            step = 4*(lbeta-1)

            !Calculating f(r)
            Overlap(lbeta)=integral_parametrization_pow(&
            obj%e_tb%p(ie1,18+step),obj%e_tb%p(ie1,19+step),&
            obj%e_tb%p(ie1,21+step),R)

            !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r))
            d_Overlaptemp(lbeta)=d_integral_parametrization_pow(&
            obj%e_tb%p(ie1,18+step),obj%e_tb%p(ie1,19+step),&
            obj%e_tb%p(ie1,21+step),R)*f_cut+Overlap(lbeta)*d_f_cut
          end do
        end select
        Overlap(:)=f_cut*Overlap(:)
        ! Calculating the derivative of the hopping NRL parameters for the
        ! general case
        do ix=1,3
          d_Overlap(:,ix)=-aa(ix)*d_Overlaptemp(:)
        end do

        !In case of alloys
        if(ie1/=ie2) then
          f_cut = fermi_function(R-obj%e_tb%r_0(ie1), 1/obj%e_tb%r_l(ie1))
          d_f_cut = fermi_function_derivative(R-obj%e_tb%r_0(ie1),1/obj%e_tb%r_l(ie1))

          select case(obj%e_tb%nrl_type(ie2))
          case('old','new','cyr')
                 
          do lbeta=1,10
            step = 4*(lbeta-1)

            !Calculating f(r) 
            I_Overlap(lbeta)=integral_parametrization_nrl(&
            obj%e_tb%p(ie2,18+step),obj%e_tb%p(ie2,19+step),&
            obj%e_tb%p(ie2,20+step),obj%e_tb%p(ie2,21+step),R)

            !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r)) 
            d_I_Overlaptemp(lbeta)=d_integral_parametrization_nrl(&
            obj%e_tb%p(ie2,18+step),obj%e_tb%p(ie2,19+step),&
            obj%e_tb%p(ie2,20+step),obj%e_tb%p(ie2,21+step),R)*f_cut+&
            I_Overlap(lbeta)*d_f_cut
          end do
        
          case('pow')
                  
          do lbeta=1,10
            step = 4*(lbeta-1)

            !Calculating f(r)
            I_Overlap(lbeta)=integral_parametrization_pow(&
            obj%e_tb%p(ie2,18+step),obj%e_tb%p(ie2,19+step),&
            obj%e_tb%p(ie2,21+step),R)

            !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r)) 
            d_I_Overlaptemp(lbeta)=d_integral_parametrization_pow(&
            obj%e_tb%p(ie2,18+step),obj%e_tb%p(ie2,19+step),&
            obj%e_tb%p(ie2,21+step),R)*f_cut+&
            I_Overlap(lbeta)*d_f_cut
          end do
          
          end select

          I_Overlap(:)=f_cut*I_Overlap(:) 
         
          Overlap(:)=0.5_rp*(Overlap(:)+I_Overlap(:))

          ! Calculating the derivative of the hopping NRL parameters, in case of
          ! alloys
          do ix=1,3
            d_Overlap(:,ix)=-aa(ix)*(d_Overlaptemp(:)+d_I_Overlaptemp(:))/2._rp
          end do
        end if

        !Lastly, computing the derivative of the hopping matrix d_B
        do ix=1,3
          call d_slater_koster(R,aa,Overlap,d_Overlap,d_Btemp,ix)
          do io1=1,obj%e%no(ie1)
            do io2=1,obj%e%no(ie2)
              d_B(ia1,in,io1,io2,ix) = d_Btemp(obj%e%o(ie1,io1),obj%e%o(ie2,io2),ix)
            end do
          end do
        end do
      end do
    end do
  
  end function build_d_b_r

  function build_s_r(obj) result(S)
    ! INPUT
    class(atom_tb),intent(in) :: obj
    ! OUTPUT
    real(rp), dimension(:,:,:,:), allocatable :: S
    ! LOCAL
    real(rp) :: RR(3),R,f_cut
    real(rp) :: temp(9,9)
    real(rp) :: Overlap(10),I_Overlap(10)
    integer :: ia1,ia2,in,ie1,ie2,io1,io2,lbeta,step

    if (.not.allocated(S)) allocate(S(obj%na,0:obj%nn_max,obj%e%no_max,obj%e%no_max))
    select case(obj%e_tb%tb_type)
    case('nrl')

    do ia1=1,obj%na
      ie1 = obj%ia2ie(ia1)
      do in=1,obj%nn(ia1)
        ia2 = obj%ian2ia(ia1,in)
        ie2 = obj%ia2ie(ia2)

        RR = obj%rn(ia1,in,:)
        R = norm2(RR(:))

        select case(obj%e_tb%nrl_type(ie1))
          case('old','cyr')
            do lbeta=1,10
              step = 4*(lbeta-1)
              Overlap(lbeta)=integral_parametrization_nrl(&
              obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
              obj%e_tb%p(ie1,60+step),obj%e_tb%p(ie1,61+step),R)
            end do
          case('new')
            do lbeta=1,10
              step = 4*(lbeta-1)
              Overlap(lbeta)=integral_parametrization_nrl_new1(&
              obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
              obj%e_tb%p(ie1,60+step),obj%e_tb%p(ie1,61+step),R)
            end do
          case('pow')
          do lbeta=1,10
            step = 4*(lbeta-1)
            Overlap(lbeta)=integral_parametrization_pow(&
            obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
            obj%e_tb%p(ie1,61+step),R)
          end do
        end select

        f_cut = fermi_function(R-obj%e_tb%r_0(ie1), 1/obj%e_tb%r_l(ie1))
        Overlap(:)=f_cut*Overlap(:)

        if(ie1/=ie2) then
          select case(obj%e_tb%nrl_type(ie2))
            case('old','cyr')
              do lbeta=1,10
                step = 4*(lbeta-1)
                I_Overlap(lbeta)=integral_parametrization_nrl(&
                obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
                obj%e_tb%p(ie2,60+step),obj%e_tb%p(ie2,61+step),R)
              end do
              
            case('new')
              do lbeta=1,10
                step = 4*(lbeta-1)
                I_Overlap(lbeta)=integral_parametrization_nrl_new1(&
                obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
                obj%e_tb%p(ie2,60+step),obj%e_tb%p(ie2,61+step),R)
              end do 
            case('pow')
            do lbeta=1,10
              step = 4*(lbeta-1)
              I_Overlap(lbeta)=integral_parametrization_pow(&
              obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
              obj%e_tb%p(ie2,61+step),R)
            end do
          end select

          f_cut = fermi_function(R-obj%e_tb%r_0(ie2), 1/obj%e_tb%r_l(ie2))
          I_Overlap(:)=f_cut*I_Overlap(:)

          do lbeta=1,10
            Overlap(lbeta)=0.5_rp*(Overlap(lbeta)+I_Overlap(lbeta))
          end do
        end if

        temp = slater_koster(RR(:)/R,Overlap)

        do io1=1,obj%e%no(ie1)
          do io2=1,obj%e%no(ie2)
            S(ia1,in,io1,io2) = temp(obj%e%o(ie1,io1),obj%e%o(ie2,io2))
          end do
        end do
      end do
    end do

    case('mod','wan')
      S=0.0_rp
    end select
  end function build_s_r

  function build_d_s_r(obj) result(d_S)
    ! INPUT
    class(atom_tb),intent(in) :: obj
    ! OUTPUT
    real(rp), dimension(:,:,:,:,:), allocatable :: d_S
    ! LOCAL
    real(rp) :: AA(3),RR(3),R,f_cut,d_f_cut,newparamtype
    real(rp) :: d_Stemp(9,9,3)
    real(rp) :: Overlap(10),d_Overlap(10,3),d_Overlaptemp(10) !SS,SP,PP(2),SD,PD(2),DD(3)
    real(rp) :: I_Overlap(10),d_I_Overlap(10,3),d_I_Overlaptemp(10) !SS,SP,PP(2),SD,PD(2),DD(3)

    integer :: ia1,ia2,in,ie1,ie2,io1,io2,lbeta,ix,step

    if(.not.allocated(d_S)) allocate(d_S(obj%na,0:obj%nn_max,obj%e%no_max,obj%e%no_max,3))

    d_Stemp = 0.0_rp ; Overlap = 0.0_rp ; d_Overlap = 0.0_rp ; d_Overlaptemp = 0.0_rp
    d_S = 0.0_rp ; I_Overlap = 0.0_rp ; d_I_Overlap = 0.0_rp ; d_I_Overlaptemp = 0.0_rp  
  
    do ia1=1,obj%na
      ie1 = obj%ia2ie(ia1)
      do in=1,obj%nn(ia1)
        ia2 = obj%ian2ia(ia1,in)
        ie2 = obj%ia2ie(ia2)

        RR = obj%rn(ia1,in,:)
        R = norm2(RR(:))
        aa(:) = RR(:)/R

        f_cut = fermi_function(R-obj%e_tb%r_0(ie1), 1/obj%e_tb%r_l(ie1))
        d_f_cut = fermi_function_derivative(R-obj%e_tb%r_0(ie1),1/obj%e_tb%r_l(ie1))

        !Getting the hopping NRL parameters in the source file and computing the
        !value of the function. For the derivative, storing part of if on
        !d_Overlaptemp to later compute it at d_Overlap
        select case(obj%e_tb%nrl_type(ie1))
        case('old','cyr')
                    
          do lbeta=1,10
            step = 4*(lbeta-1)

            !Calculating f(r)
            Overlap(lbeta)=integral_parametrization_nrl(&
            obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
            obj%e_tb%p(ie1,60+step),obj%e_tb%p(ie1,61+step),R)

            !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r))
            d_Overlaptemp(lbeta)=d_integral_parametrization_nrl(&
            obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
            obj%e_tb%p(ie1,60+step),obj%e_tb%p(ie1,61+step),R)*f_cut+&
            Overlap(lbeta)*d_f_cut
          end do
          
        case('new')
          !Calculating f(r)          
          do lbeta=1,10
            step = 4*(lbeta-1)
           if(lbeta==2.or.lbeta==5.or.lbeta==6.or.lbeta==7)then
            Overlap(lbeta)=integral_parametrization_nrl_new2(&
            obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
            obj%e_tb%p(ie1,60+step),obj%e_tb%p(ie1,61+step),R)
           else  
            Overlap(lbeta)=integral_parametrization_nrl_new1(&
            obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
            obj%e_tb%p(ie1,60+step),obj%e_tb%p(ie1,61+step),R)
           end if
          end do
          !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r)) 
          do lbeta=1,10
            step = 4*(lbeta-1)
            newparamtype=1._rp
           if(lbeta==2.or.lbeta==5.or.lbeta==6.or.lbeta==7)then
            newparamtype=0._rp
           end if
            d_Overlaptemp(lbeta)=d_integral_parametrization_nrl_new(&
            obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
            obj%e_tb%p(ie1,60+step),obj%e_tb%p(ie1,61+step),R,newparamtype)+&
            Overlap(lbeta)*d_f_cut
          end do
        case('pow')
          !Calculating f(r)          
          do lbeta=1,10
            step = 4*(lbeta-1)
            Overlap(lbeta)=integral_parametrization_pow(&
            obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
            obj%e_tb%p(ie1,61+step),R)
          end do
          !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r)) 
          do lbeta=1,10
            step = 4*(lbeta-1)
            d_Overlaptemp(lbeta)=d_integral_parametrization_pow(&
            obj%e_tb%p(ie1,58+step),obj%e_tb%p(ie1,59+step),&
            obj%e_tb%p(ie1,61+step),R)*f_cut+&
            Overlap(lbeta)*d_f_cut
          end do
        end select

        Overlap(:)=f_cut*Overlap(:)
        ! Calculating the derivative of the hopping NRL parameters for the
        ! general case
        do ix=1,3
         d_Overlap(:,ix)=-aa(ix)*d_Overlaptemp(:)
        end do

        !In case of alloys
        if(ie1/=ie2) then

          f_cut = fermi_function(R-obj%e_tb%r_0(ie1), 1/obj%e_tb%r_l(ie1))
        d_f_cut = fermi_function_derivative(R-obj%e_tb%r_0(ie1),1/obj%e_tb%r_l(ie1))

          select case(obj%e_tb%nrl_type(ie2))
          case('old','cyr')
          !Calculating f(r)          
          do lbeta=1,10
            step = 4*(lbeta-1)
            I_Overlap(lbeta)=integral_parametrization_nrl(&
            obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
            obj%e_tb%p(ie2,60+step),obj%e_tb%p(ie2,61+step),R)
          end do
          !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r)) 
          do lbeta=1,10
            step = 4*(lbeta-1)
            d_I_Overlaptemp(lbeta)=d_integral_parametrization_nrl(&
            obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
            obj%e_tb%p(ie2,60+step),obj%e_tb%p(ie2,61+step),R)*f_cut+&
            I_Overlap(lbeta)*d_f_cut
          end do
          case('new')
          !Calculating f(r)          
          do lbeta=1,10
            step = 4*(lbeta-1)
           if(lbeta==2.or.lbeta==5.or.lbeta==6.or.lbeta==7)then
            I_Overlap(lbeta)=integral_parametrization_nrl_new2(&
            obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
            obj%e_tb%p(ie2,60+step),obj%e_tb%p(ie2,61+step),R)
           else
            I_Overlap(lbeta)=integral_parametrization_nrl_new1(&
            obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
            obj%e_tb%p(ie2,60+step),obj%e_tb%p(ie2,61+step),R)
           end if
          end do
          !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r)) 
          do lbeta=1,10
            step = 4*(lbeta-1)
            newparamtype=1._rp
           if(lbeta==2.or.lbeta==5.or.lbeta==6.or.lbeta==7)then
            newparamtype=0._rp
           end if
            d_I_Overlaptemp(lbeta)=d_integral_parametrization_nrl_new(&
            obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
            obj%e_tb%p(ie2,60+step),obj%e_tb%p(ie2,61+step),R,newparamtype)+&
            Overlap(lbeta)*d_f_cut
          end do
          case('pow')
          !Calculating f(r)          
          do lbeta=1,10
            step = 4*(lbeta-1)
            I_Overlap(lbeta)=integral_parametrization_pow(&
            obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
            obj%e_tb%p(ie2,61+step),R)
          end do
          !The derivative as (f'(r)*fcut(r)+f(r)*fcut'(r)) 
          do lbeta=1,10
            step = 4*(lbeta-1)
            d_I_Overlaptemp(lbeta)=d_integral_parametrization_pow(&
            obj%e_tb%p(ie2,58+step),obj%e_tb%p(ie2,59+step),&
            obj%e_tb%p(ie2,61+step),R)*f_cut+&
            I_Overlap(lbeta)*d_f_cut
          end do
          end select

          I_Overlap(:)=f_cut*I_Overlap(:)

          Overlap(:)=0.5_rp*(Overlap(:)+I_Overlap(:))

         ! Calculating the derivative of the hopping NRL parameters, in case of
         ! alloys
         do ix=1,3
          d_Overlap(:,ix)=-aa(ix)*(d_Overlaptemp(:)+d_I_Overlaptemp(:))/2._rp
         end do
        end if

       !Lastly, computing the derivative of the hopping matrix d_S
       do ix=1,3
        call d_slater_koster(R,aa,Overlap,d_Overlap,d_Stemp,ix)
        do io1=1,obj%e%no(ie1)
          do io2=1,obj%e%no(ie2)
            d_S(ia1,in,io1,io2,ix) = d_Stemp(obj%e%o(ie1,io1),obj%e%o(ie2,io2),ix)
          end do
        end do
       end do
      end do
    end do
  end function build_d_s_r

  function build_en_intra(obj) result(en_intra)
    ! INPUT
    class(atom_tb),intent(in) :: obj
    ! OUTPUT
    real(rp), dimension(:,:), allocatable :: en_intra
    ! LOCAL
    real(rp) :: rho,r,r_0,r_l,f_cut,en_intra_s,en_intra_p,en_intra_d
    integer  :: ia1,ia2,in,ie1,ie2,io

    write(output_unit,*) 'DEBUG == Entering atom_tb & build_en_intra'
    call TBKOSTER_flush(output_unit)

    if (.not.allocated(en_intra)) allocate(en_intra(obj%na,obj%e%no_max))
    en_intra = 0.0_rp

    do ia1=1,obj%na
      ie1 = obj%ia2ie(ia1)
      rho = 0.0_rp
      do in=1,obj%nn(ia1)
        ia2 = obj%ian2ia(ia1,in)
        ie2 = obj%ia2ie(ia2)
        r = norm2(obj%rn(ia1,in,:))
        r_0 = (obj%e_tb%r_0(ie2) + obj%e_tb%r_0(ie1))/2._rp
        r_l = (obj%e_tb%r_l(ie2) + obj%e_tb%r_l(ie1))/2._rp
        f_cut = fermi_function(r-r_0,1._rp/r_l)
        rho = rho + f_cut*exp(-obj%e_tb%p(ie2,1)**2*r)
        !rho = rho + f_cut*exp(-0.25*(obj%e_tb%p(ie2,1)+obj%e_tb%p(ie1,1))**2*r)
      end do
      en_intra_s = onsite(obj%e_tb%p(ie1,2), obj%e_tb%p(ie1,3), &
       obj%e_tb%p(ie1,4), obj%e_tb%p(ie1,5), obj%e_tb%p(ie1,6), rho)
      en_intra_p = onsite(obj%e_tb%p(ie1,7), obj%e_tb%p(ie1,8), &
       obj%e_tb%p(ie1,9), obj%e_tb%p(ie1,10), obj%e_tb%p(ie1,11), rho)
      en_intra_d = onsite(obj%e_tb%p(ie1,12), obj%e_tb%p(ie1,13), &
       obj%e_tb%p(ie1,14), obj%e_tb%p(ie1,15), obj%e_tb%p(ie1,16), rho)

      do io=1,obj%e%no(ie1)
        if(obj%e%o(ie1,io)==1) then
          en_intra(ia1,io) = en_intra_s
        elseif(obj%e%o(ie1,io)>=2 .and. obj%e%o(ie1,io)<=4) then
          en_intra(ia1,io) = en_intra_p
        elseif(obj%e%o(ie1,io)>=5 .and. obj%e%o(ie1,io)<=9) then
          en_intra(ia1,io) = en_intra_d
        end if
      end do
    end do
    write(output_unit,*) 'DEBUG == onsite sum(en_intra(1,:))=',sum(en_intra(1,:))
    write(output_unit,*) 'DEBUG == Exiting atom_tb & build_en_intra'
    call TBKOSTER_flush(output_unit)
    
  end function build_en_intra

  ! Function build_d_en_intra to be used on the molecular dynamics
  function build_d_en_intra(obj) result(d_en_intra)                  
    ! INPUT
    class(atom_tb),intent(in) :: obj
    ! OUTPUT
    real(rp),dimension(obj%na,obj%e%no_max,3) :: d_en_intra 
    ! LOCAL
    real(rp) :: rho,r,r_0,r_l,f_cut,d_en_intra_s,d_en_intra_p,d_en_intra_d,t1,t2
    real(rp),dimension(3) :: rr
    integer  :: ia1,ia2,in,ie1,ie2,io,ix

    d_en_intra = 0.0_rp
    
    do ia1=1,obj%na
      ie1 = obj%ia2ie(ia1)
      rho = 0.0_rp

      do in=1,obj%nn(ia1)
        ia2 = obj%ian2ia(ia1,in)
        ie2 = obj%ia2ie(ia2)
        rr(:) = obj%rn(ia1,in,:)
        r = norm2(rr(:))
        r_0 = (obj%e_tb%r_0(ie2) + obj%e_tb%r_0(ie1))/2._rp
        r_l = (obj%e_tb%r_l(ie2) + obj%e_tb%r_l(ie1))/2._rp
        f_cut = fermi_function(r-r_0,1._rp/r_l)

        t1 = (obj%e_tb%p(ie2,1))**2
        rho = rho + f_cut*exp(-t1*r)
        t2 = f_cut*exp(-t1*r)*(t1+(1._rp-f_cut)/r_l)
      end do
     
      do ix=1,3
        d_en_intra_s = d_onsite(obj%e_tb%p(ie1,2), obj%e_tb%p(ie1,3), &
        obj%e_tb%p(ie1,4), obj%e_tb%p(ie1,5), obj%e_tb%p(ie1,6), rho)
        d_en_intra_p = d_onsite(obj%e_tb%p(ie1,7), obj%e_tb%p(ie1,8), &
        obj%e_tb%p(ie1,9), obj%e_tb%p(ie1,10), obj%e_tb%p(ie1,11), rho)
        d_en_intra_d = d_onsite(obj%e_tb%p(ie1,12), obj%e_tb%p(ie1,13), &
        obj%e_tb%p(ie1,14), obj%e_tb%p(ie1,15), obj%e_tb%p(ie1,16), rho)
        do io=1,obj%e%no(ie1)
          if(obj%e%o(ie1,io)==1) then
            d_en_intra(ia1,io,ix) = d_en_intra_s
          elseif(obj%e%o(ie1,io)>=2 .and. obj%e%o(ie1,io)<=4) then
            d_en_intra(ia1,io,ix) = d_en_intra_p
          elseif(obj%e%o(ie1,io)>=5 .and. obj%e%o(ie1,io)<=9) then
            d_en_intra(ia1,io,ix) = d_en_intra_d
          end if
        end do
      end do
    end do
    !write(*,*)'d_onesite',sum(d_en_intra(1,:,1)),sum(d_en_intra(1,:,2)),sum(d_en_intra(1,:,3))
  end function build_d_en_intra

  function integral_parametrization_nrl(e,f,fbar,g,R) result(p)
    real(rp),intent(in) :: e,f,fbar,g,R
    real(rp) :: p
    p = (e+f*R+fbar*R*R)*exp(-g*g*R)
  end function integral_parametrization_nrl

  !Function d_integral_parametrization_nrl to be used in the molecular dynamics
  function d_integral_parametrization_nrl(e,f,fbar,g,R) result(p)
    real(rp),intent(in) :: e,f,fbar,g,R          
    real(rp) :: p

    p = exp(-g*g*R)*((f+2._rp*fbar*R)-(e+f*R+fbar*R*R)*(g*g))
  end function d_integral_parametrization_nrl

  function integral_parametrization_nrl_new1(e,f,fbar,g,R) result(p)
    real(rp),intent(in) :: e,f,fbar,g,R
    real(rp) :: p, R2, R3

    R2=R*R
    R3=R*R2

    p = (1._rp+e*R+f*R2+fbar*R3)*exp(-g*g*R)
  end function integral_parametrization_nrl_new1

  function integral_parametrization_nrl_new2(e,f,fbar,g,R) result(p)
    real(rp),intent(in) :: e,f,fbar,g,R
    real(rp) :: p, R2, R3

    R2=R*R
    R3=R*R2

    p = (e*R+f*R2+fbar*R3)*exp(-g*g*R)
  end function integral_parametrization_nrl_new2

  !Function d_integral_parametrization_nrl_new to be used in the molecular
  !dynamics
  function d_integral_parametrization_nrl_new(e,f,fbar,g,R,newparam) result(p)
    real(rp),intent(in) :: e,f,fbar,g,R,newparam
    real(rp) :: p, R2, R3, g2

    g2=g*g
    R2=R*R
    R3=R*R2

    p=exp(-g2*R)*((e+2._rp*f*R+3._rp*fbar*R2)-(newparam+e*R+f*R2+fbar*R3)*g2) 
  end function d_integral_parametrization_nrl_new

  function integral_parametrization_pow(e,f,g,R) result(p)
    real(rp),intent(in) :: e,f,g,R
    real(rp) :: p

    if(e==0.0_rp) then
      p = 0.0_rp
    else
      p = e*(f/R)**(g)
    endif
  end function integral_parametrization_pow

  ! Function d_integral_parametrization_pow to be used in the molecular dynamics  
  function d_integral_parametrization_pow(e,f,g,R) result(p)
    real(rp),intent(in) :: e,f,g,R
    real(rp) :: p

    if(e==0.0_rp) then
      p = 0.0_rp
    else
      p = e*(f/R)**(g)*(-g/R)
    endif
  end function d_integral_parametrization_pow

  function onsite(c,b1,b2,b3,b4,rho) result(p)
    real(rp), intent(in) :: c,b1,b2,b3,b4,rho
    real(rp) :: p, t1, t2, t3

    t1 = rho**(one_third)
    t2 = t1*t1            ! rho**(2/3)
    t3 = t2*t2            ! rho**(4/3)

    p = c + b1*t1 + b2*t2 + b3*t3 + b4*rho*rho
  end function onsite

  ! Function d_onsite to be used in the molecular dynamics
  ! dp/dx = dp/drho * drho/dx
  function d_onsite(c,b1,b2,b3,b4,rho) result(p)
    real(rp), intent(in) :: c,b1,b2,b3,b4,rho
    real(rp) :: p, t1, t2

    t1 = rho**(one_third)
    t2 = t1*t1           ! rho**(2/3)

    p = (one_third*b1/t2 + two_third*b2/t1 + (4._rp/3._rp)*b3*t1 + 2._rp*b4*rho)
  end function d_onsite

  ! Build the Slater Koster transformation
  function slater_koster(R,Overlap) result(H)
    real(rp), intent(in) :: R(3) ! X,Y,Z Direction cosines
    real(rp), intent(in) :: Overlap(10) ! SS,SP,PP(2),SD,PD(2),DD(3)
    real(rp) :: H(9,9)
    real(rp) :: X,Y,Z,XX,YY,ZZ
    real(rp) :: SS,SP,SD,PP(2),PD(2),DD(3)
    integer :: k,km,l

    X=R(1); Y=R(2); Z=R(3) 

    ! List of bond integrals
    SS=Overlap(1)    ! V_{ss\sigma}
    SP=Overlap(2)    ! V_{sp\sigma}
    PP(1)=Overlap(3) ! V_{pp\sigma}
    PP(2)=Overlap(4) ! V_{pp\pi}
    SD=Overlap(5)    ! V_{sd\sigma}
    PD(1)=Overlap(6) ! V_{pd\sigma}
    PD(2)=Overlap(7) ! V_{pd\pi}
    DD(1)=Overlap(8) ! V_{dd\sigma}
    DD(2)=Overlap(9) ! V_{dd\pi}
    DD(3)=Overlap(10)! V_{dd\delta}

    XX=X*X; YY=Y*Y; ZZ=Z*Z

    !
    ! 1=s; 2=px; 3=py; 4=pz; 5=dxy; 6=dyz; 7=dzx; 8=dx^2-y^2; 9=d(3z^2-r^2)
    ! See Wikipedia Tight Binding

    ! ss coupling
    H(1,1)=SS          ! E_{s,s}

    ! sp coupling
    H(1,2)=X*SP        ! E_{s,x}
    H(1,3)=Y*SP        ! E_{s,y}
    H(1,4)=Z*SP        ! E_{s,z}

    H(2,1)=-H(1,2)
    H(3,1)=-H(1,3)
    H(4,1)=-H(1,4)

    !  sd coupling
    H(1,5)=sqrt_three*X*Y*SD             ! E_{s,xy}
    H(1,6)=sqrt_three*Y*Z*SD             ! E_{s,yz}
    H(1,7)=sqrt_three*Z*X*SD             ! E_{s,zx}
    H(1,8)=0.5_rp*sqrt_three*(XX-YY)*SD  ! E_{s,x^2-y^2}
    H(1,9)=(ZZ-0.5_rp*(XX+YY))*SD        ! E_{s,3z^2-r^2}

    H(5,1)=H(1,5)
    H(6,1)=H(1,6)
    H(7,1)=H(1,7)
    H(8,1)=H(1,8)
    H(9,1)=H(1,9)

    ! pp coupling
    H(2,2)=XX*PP(1)+(1._rp-XX)*PP(2)     ! E_{x,x}
    H(2,3)=X*Y*(PP(1)-PP(2))             ! E_{x,y}
    H(2,4)=X*Z*(PP(1)-PP(2))             ! E_{x,z}
    H(3,3)=YY*PP(1)+(1._rp-YY)*PP(2)     ! E_{y,y}
    H(3,4)=Y*Z*(PP(1)-PP(2))             ! E_{y,z}
    H(4,4)=ZZ*PP(1)+(1._rp-ZZ)*PP(2)     ! E_{z,z}

    do k=2,3
      km=k-1
      do l=1,km
        H(k+1,l+1)=H(l+1,k+1)
      end do
    end do

    ! pd coupling
    H(2,5)=Y*(sqrt_three*XX*PD(1)+(1._rp-2._rp*XX)*PD(2))         ! E_{x,xy}
    H(2,6)=X*Y*Z*(sqrt_three*PD(1)-2._rp*PD(2))                   ! E_{x,yz}
    H(2,7)=Z*(sqrt_three*XX*PD(1)+(1._rp-2._rp*XX)*PD(2))         ! E_{x,zx}
    H(2,8)=X*(0.5_rp*sqrt_three*(XX-YY)*PD(1)+(1._rp-XX+YY)*PD(2))! E_{x,x^2-y^2)
    H(2,9)=X*((ZZ-0.5_rp*(XX+YY))*PD(1)-sqrt_three*ZZ*PD(2))    ! E_{x,3z^2-r^2}
    H(3,5)=X*(sqrt_three*YY*PD(1)+(1._rp-2._rp*YY)*PD(2))       ! E_{y,xy}
    H(3,6)=Z*(sqrt_three*YY*PD(1)+(1._rp-2._rp*YY)*PD(2))       ! E_{y,yz}
    H(3,7)=H(2,6)      ! E_{y,zx}=E_{x,yz}=E_{z,xy} (circular permutation)
    H(3,8)=Y*(0.5_rp*sqrt_three*(XX-YY)*PD(1)-(1._rp+XX-YY)*PD(2))!E_{y,x^2-y^2}
    H(3,9)=Y*((ZZ-0.5_rp*(XX+YY))*PD(1)-sqrt_three*ZZ*PD(2)) ! E_{y,3z^2-r^2}
    H(4,5)=H(2,6)
    H(4,6)=Y*(sqrt_three*ZZ*PD(1)+(1._rp-2._rp*ZZ)*PD(2)) ! E_{z,yz}
    H(4,7)=X*(sqrt_three*ZZ*PD(1)+(1._rp-2._rp*ZZ)*PD(2)) ! E_{z,xz}
    H(4,8)=Z*(XX-YY)*(0.5_rp*sqrt_three*PD(1)-PD(2))      ! E_{z,x^2-y^2}
    H(4,9)=Z*((ZZ-0.5_rp*(XX+YY))*PD(1)+sqrt_three*(XX+YY)*PD(2))!E_{z,3z^2-r^2}

    do k=1,3
      do l=1,5
        H(l+4,k+1)=-H(k+1,l+4)
      end do
    end do

    ! dd coupling
    ! E_{xy,xy}
    H(5,5)=3._rp*XX*YY*DD(1)+(XX+YY-4._rp*XX*YY)*DD(2)+(ZZ+XX*YY)*DD(3)
    ! E_{xy,yz}
    H(5,6)=X*Z*(3._rp*YY*DD(1)+(1._rp-4._rp*YY)*DD(2)+(YY-1._rp)*DD(3))
    ! E_{xy,zx}
    H(5,7)=Y*Z*(3._rp*XX*DD(1)+(1._rp-4._rp*XX)*DD(2)+(XX-1._rp)*DD(3))
    ! E_{xy,x^2-y^2}
    H(5,8)=(1.5_rp*DD(1)-2._rp*DD(2)+0.5_rp*DD(3))*(XX-YY)*X*Y
    ! E_{xy,3z^2-r^2}
    H(5,9)=X*Y*sqrt_three*((ZZ-(XX+YY)/2._rp)*DD(1)-2._rp*ZZ*DD(2)+0.5_rp*(1._rp+ZZ)*DD(3))
    ! E_{yz,yz}
    H(6,6)=3._rp*YY*ZZ*DD(1)+(YY+ZZ-4._rp*YY*ZZ)*DD(2)+(XX+YY*ZZ)*DD(3)
    ! E_{yz,zx}
    H(6,7)=X*Y*(3._rp*ZZ*DD(1)+(1._rp-4._rp*ZZ)*DD(2)+(ZZ-1._rp)*DD(3))
    ! E_{yz,x^2-y^2}
    H(6,8)=Y*Z*(1.5_rp*(XX-YY)*DD(1)-(1._rp+2._rp*(XX-YY))*DD(2)+(1._rp+(XX-YY)/2._rp)*DD(3))
    ! E_{yz,3z^2-r^2}
    H(6,9)=Y*Z*sqrt_three*((ZZ-(XX+YY)/2._rp)*DD(1)+(XX+YY-ZZ)*DD(2)-0.5_rp*(XX+YY)*DD(3))
    ! E_{zx,zx}
    H(7,7)=3._rp*XX*ZZ*DD(1)+(XX+ZZ-4._rp*XX*ZZ)*DD(2)+(YY+ZZ*XX)*DD(3)
    ! E_{zx,x^2-y^2}
    H(7,8)=Z*X*(1.5_rp*(XX-YY)*DD(1)+(1._rp-2._rp*(XX-YY))*DD(2)-(1._rp-(XX-YY)/2._rp)*DD(3))
    ! E_{zx,3z^2-r^2}
    H(7,9)=Z*X*sqrt_three*((ZZ-(XX+YY)/2._rp)*DD(1)+(XX+YY-ZZ)*DD(2)-0.5_rp*(XX+YY)*DD(3))
    ! E_{x^2-y^2,x^2-y^2}
    H(8,8)=3._rp/4._rp*(XX-YY)**2*DD(1)+(XX+YY-(XX-YY)**2)*DD(2)+((XX-YY)**2/4._rp+ZZ)*DD(3)
    ! E_{x^2-y^2,3z^2-r^2}
    H(8,9)=sqrt_three*(XX-YY)*(0.5_rp*(ZZ-(XX+YY)/2._rp)*DD(1)-ZZ*DD(2)+0.25_rp*(1._rp+ZZ)*DD(3))
    ! E_{3z^2-r^2,3z^2-r^2}
    H(9,9)=(ZZ-(XX+YY)/2._rp)**2*DD(1)+3._rp*ZZ*(XX+YY)*DD(2)+3._rp/4._rp*(XX+YY)**2*DD(3)

    do k=2,5
      km=k-1
      do l=1,km
        H(k+4,l+4)=H(l+4,k+4)
      end do
    end do
  end function slater_koster

  ! Build the derivative of the Slater Koster transformation
  subroutine d_slater_koster(r,AA,Overlap,d_Overlap,h,ix)

    real(rp), intent(out) :: h(9,9,3)
    real(rp), intent(in)  :: AA(3),Overlap(10),d_Overlap(10,3),r
    integer, intent(in)   :: ix

    real(rp) :: SS,SP,SD,PP(2),PD(2),DD(3)
    real(rp) :: d_SS(3),d_SP(3),d_SD(3),d_PP(2,3),d_PD(2,3),d_DD(3,3)
    real(rp) :: cx,cy,cz,cx2,cy2,cz2,cz3,d_cx,d_cy,d_cz
    real(rp) :: X, Y, Z
    integer :: k,km,l

    h(:,:,ix)=0.0_rp
    
    !
    !      1=s; 2=px; 3=py; 4=pz; 5=dxy; 6=dyz; 7=dzx; 8=dx^2-y^2; 9=d(3z^2-r^2)
    !

    X=AA(1); Y=AA(2); Z=AA(3)

    SS=Overlap(1)
    SP=Overlap(2)
    PP(1)=Overlap(3)
    PP(2)=Overlap(4)
    SD=Overlap(5)
    PD(1)=Overlap(6)
    PD(2)=Overlap(7)
    DD(1)=Overlap(8)
    DD(2)=Overlap(9)
    DD(3)=Overlap(10)

    d_SS(:)=d_Overlap(1,:)
    d_SP(:)=d_Overlap(2,:)
    d_PP(1,:)=d_Overlap(3,:)
    d_PP(2,:)=d_Overlap(4,:)
    d_SD=d_Overlap(5,:)
    d_PD(1,:)=d_Overlap(6,:)
    d_PD(2,:)=d_Overlap(7,:)
    d_DD(1,:)=d_Overlap(8,:)
    d_DD(2,:)=d_Overlap(9,:)
    d_DD(3,:)=d_Overlap(10,:)

    cx=x
    cy=y
    cz=z
    cx2=cx*cx
    cy2=cy*cy
    cz2=cz*cz
    cz3=cz*cz*cz

    select case (ix)
    case (1)
      d_cx = -(1._rp-cx2)/r
      d_cy = +(cx*cy)/r
      d_cz = +(cx*cz)/r
    case (2)
      d_cx = +(cx*cy)/r
      d_cy = -(1._rp-cy2)/r
      d_cz = +(cy*cz)/r
    case (3)
      d_cx = +(cx*cz)/r
      d_cy = +(cy*cz)/r
      d_cz = -(1._rp-cz2)/r
    end select

    ! ss coupling
    h(1,1,ix) = d_SS(ix)

    ! sp coupling
    h(1,2,ix) = d_cx*SP+cx*d_SP(ix)
    h(1,3,ix) = d_cy*SP+cy*d_SP(ix)
    h(1,4,ix) = d_cz*SP+cz*d_SP(ix)

    h(2,1,ix) = -h(1,2,ix)
    h(3,1,ix) = -h(1,3,ix)
    h(4,1,ix) = -h(1,4,ix)

    !  sd coupling

    h(1,5,ix) = sqrt_three*(d_cx*cy*SD+cx*d_cy*SD+cx*cy*d_SD(ix))
    h(1,6,ix) = sqrt_three*(d_cy*cz*SD+cy*d_cz*SD+cy*cz*d_SD(ix))
    h(1,7,ix) = sqrt_three*(d_cz*cx*SD+cz*d_cx*SD+cz*cx*d_SD(ix))
    h(1,8,ix) = 0.5_rp*sqrt_three*(2._rp*(cx*d_cx-cy*d_cy)*SD+(cx2-cy2)*&
    d_SD(ix))
    h(1,9,ix) = (2._rp*cz*d_cz-(cx*d_cx+cy*d_cy))*SD+(cz2-0.5_rp*(cx2+cy2))&
    *d_SD(ix)

    h(5,1,ix)=h(1,5,ix)
    h(6,1,ix)=h(1,6,ix)
    h(7,1,ix)=h(1,7,ix)
    h(8,1,ix)=H(1,8,ix)
    h(9,1,ix)=h(1,9,ix)
    
    ! pp coupling
    
    h(2,2,ix) = 2._rp*cx*d_cx*PP(1)+cx2*d_PP(1,ix)-2._rp*cx*d_cx*PP(2)+&
    (1._rp-cx2)*d_PP(2,ix)
    h(2,3,ix) = (d_cx*cy+cx*d_cy)*(PP(1)-PP(2))+cx*cy*(d_PP(1,ix)-d_PP(2,ix))
    h(2,4,ix) = d_cx*cz*(PP(1)-PP(2))+cx*d_cz*(PP(1)-PP(2))+cx*cz*(d_PP(1,ix)-&
    d_PP(2,ix))
    h(3,3,ix) = 2._rp*cy*d_cy*PP(1)+cy2*d_PP(1,ix)-2._rp*cy*d_cy*PP(2)+&
    (1._rp-cy2)*d_PP(2,ix)
    h(3,4,ix) = (d_cy*cz+cy*d_cz)*(PP(1)-PP(2))+cy*cz*(d_PP(1,ix)-d_PP(2,ix))
    h(4,4,ix) = 2._rp*cz*d_cz*PP(1)+cz2*d_PP(1,ix)-2._rp*cz*d_cz*PP(2)+&
    (1._rp-cz2)*d_PP(2,ix)

    do k=2,3
      km=k-1
      do l=1,km
        h(k+1,l+1,ix)=h(l+1,k+1,ix)
      end do
    end do

    ! pd coupling

    h(2,5,ix) = (2._rp*cx*d_cx*cy+cx2*d_cy)*(sqrt_three*PD(1)-&
    2._rp*PD(2))+cx2*cy*(sqrt_three*d_PD(1,ix)-2._rp*d_PD(2,ix))+&
    d_cy*PD(2)+cy*d_PD(2,ix)
    h(2,6,ix) = (d_cx*cy*cz+cx*d_cy*cz+cx*cy*d_cz)*(sqrt_three*PD(1)-2._rp* &
    PD(2))+cx*cy*cz*(sqrt_three*d_PD(1,ix)-2._rp*d_PD(2,ix))
    h(2,7,ix) = (2._rp*cx*d_cx*cz+cx2*d_cz)*(sqrt_three*PD(1)-2._rp* &
    PD(2))+cx2*cz*(sqrt_three*d_PD(1,ix)-2._rp*d_PD(2,ix))+d_cz*PD(2)+cz* &
    d_PD(2,ix)
    h(2,8,ix) = (d_cx*(cx2-cy2)+2._rp*cx*(cx*d_cx-cy*d_cy))*(0.5_rp*sqrt_three* &
    PD(1)-PD(2))+cx*(cx2-cy2)*(0.5_rp*sqrt_three*d_PD(1,ix)-d_PD(2,ix))+d_cx* &
    PD(2)+cx*d_PD(2,ix)
    h(2,9,ix) = (d_cx*cz2+2._rp*cx*cz*d_cz)*(PD(1)-sqrt_three*PD(2))+cx*cz2*&
    (d_PD(1,ix)-sqrt_three*d_PD(2,ix))-0.5_rp*d_cx*(cx2+cy2)*PD(1)-cx*(cx* & 
    d_cx+cy*d_cy)*PD(1)-0.5_rp*cx*(cx2+cy2)*d_PD(1,ix)
    h(3,5,ix) = (d_cx*cy2+cx*2._rp*cy*d_cy)*(sqrt_three*PD(1)-2._rp*PD(2))+cx*cy2*&
    (sqrt_three*d_PD(1,ix)-2._rp*d_PD(2,ix))+d_cx*PD(2)+cx*d_PD(2,ix)
    h(3,6,ix) = (d_cz*cy2+cz*2._rp*cy*d_cy)*(sqrt_three*PD(1)-2._rp*PD(2))+cz*cy2*&
    (sqrt_three*d_PD(1,ix)-2._rp*d_PD(2,ix))+d_cz*PD(2)+cz*d_PD(2,ix)
    h(3,7,ix) = h(2,6,ix)
    h(3,8,ix) = (d_cy*(cx2-cy2)+2._rp*cy*(cx*d_cx-cy*d_cy))*&
    (0.5_rp*sqrt_three*PD(1)-PD(2))+cy*(cx2-cy2)*(0.5_rp*sqrt_three*d_PD(1,ix)-&
    d_PD(2,ix))-d_cy*PD(2)-cy*d_PD(2,ix)
    h(3,9,ix) = (d_cy*cz2+2._rp*cy*cz*d_cz)*(PD(1)-sqrt_three*PD(2))+cy*cz2*&
    (d_PD(1,ix)-sqrt_three*d_PD(2,ix))-0.5_rp*d_cy*(cx2+cy2)*PD(1)-cy*&
    (cx*d_cx+cy*d_cy)*PD(1)-0.5_rp*cy*(cx2+cy2)*d_PD(1,ix)
    h(4,5,ix) = h(2,6,ix)
    h(4,6,ix) = (d_cy*cz2+cy*2._rp*cz*d_cz)*(sqrt_three*PD(1)-2._rp*PD(2))&
    +cy*cz2*(sqrt_three*d_PD(1,ix)-2._rp*d_PD(2,ix))+d_cy*PD(2)+cy*d_PD(2,ix)
    h(4,7,ix) = (d_cx*cz2+cx*2._rp*cz*d_cz)*(sqrt_three*PD(1)-2._rp*PD(2))&
    +cz2*cx*(sqrt_three*d_PD(1,ix)-2._rp*d_PD(2,ix))+d_cx*PD(2)+cx*d_PD(2,ix)
    h(4,8,ix) = (d_cz*(cx2-cy2)+cz*2._rp*(cx*d_cx-cy*d_cy))*&
    (0.5_rp*sqrt_three*PD(1)-PD(2))+cz*(cx2-cy2)*(0.5_rp*sqrt_three*d_PD(1,ix)-&
    d_PD(2,ix))
    h(4,9,ix) = 3._rp*cz2*d_cz*PD(1)+cz3*d_PD(1,ix)+(d_cz*(cx2+cy2)+cz*2._rp*&
    (cx*d_cx+cy*d_cy))*(sqrt_three*PD(2)-0.5_rp*PD(1))+cz*(cx2+cy2)*&
    (sqrt_three*d_PD(2,ix)-0.5_rp*d_PD(1,ix))
    
    do k=1,3
      do l=1,5
        h(l+4,k+1,ix)=-h(k+1,l+4,ix)
      end do
    end do

    ! dd coupling

    h(5,5,ix) = 2._rp*(cx*d_cx*cy2+cx2*cy*d_cy)*(3._rp*DD(1)-4._rp*DD(2)+&
    DD(3))+cx2*cy2*(3._rp*d_DD(1,ix)-4._rp*d_DD(2,ix)+d_DD(3,ix))+2._rp*&
    (cx*d_cx+cy*d_cy)*DD(2)+(cx2+cy2)*d_DD(2,ix)+2._rp*cz*d_cz*DD(3)+cz2*&
    d_DD(3,ix)
    h(5,6,ix) = (d_cx*cy2*cz+2._rp*cx*cy*d_cy*cz+cx*cy2*d_cz)*&
    (3._rp*DD(1)-4._rp*DD(2)+DD(3))+cx*cy2*cz*(3._rp*d_DD(1,ix)-4._rp*&
    d_DD(2,ix)+d_DD(3,ix))+(d_cx*cz+cx*d_cz)*(DD(2)-DD(3))+cx*cz*&
    (d_DD(2,ix)-d_DD(3,ix))
    h(5,7,ix) = (2._rp*cx*d_cx*cy*cz+cx2*d_cy*cz+cx2*cy*d_cz)*&
    (3._rp*DD(1)-4._rp*DD(2)+DD(3))+cx2*cy*cz*(3._rp*d_DD(1,ix)-4._rp*&
    d_DD(2,ix)+d_DD(3,ix))+(d_cy*cz+cy*d_cz)*(DD(2)-DD(3))+cy*cz*&
    (d_DD(2,ix)-d_DD(3,ix))
    h(5,8,ix) = (d_cx*cy*(cx2-cy2)+cx*d_cy*(cx2-cy2)+2._rp*cx*cy*&
    (cx*d_cx-cy*d_cy))*(1.5*DD(1)-2._rp*DD(2)+0.5_rp*DD(3))+cx*cy*&
    (cx2-cy2)*(1.5*d_DD(1,ix)-2._rp*d_DD(2,ix)+0.5_rp*d_DD(3,ix))
    h(5,9,ix) = sqrt_three*(d_cx*cy*cz2+cx*d_cy*cz2+2._rp*cx*cy*cz*d_cz)*&
    (DD(1)-2._rp*DD(2)+0.5_rp*DD(3))+sqrt_three*cx*cy*cz2*(d_DD(1,ix)-2._rp*&
    d_DD(2,ix)+0.5_rp*d_DD(3,ix))-0.5_rp*sqrt_three*(d_cx*cy+cx*d_cy)*&
    ((cx2+cy2)*DD(1)-DD(3))-0.5_rp*sqrt_three*cx*cy*(2._rp*(cx*d_cx+cy*d_cy)*&
    DD(1)+(cx2+cy2)*d_DD(1,ix)-d_DD(3,ix))
    h(6,6,ix) = 2._rp*cy*cz*(cz*d_cy+d_cz*cy)*(3._rp*DD(1)-4._rp*DD(2)+&
    DD(3))+cy2*cz2*(3._rp*d_DD(1,ix)-4._rp*d_DD(2,ix)+d_DD(3,ix))+2._rp*&
    (cy*d_cy+cz*d_cz)*DD(2)+(cy2+cz2)*d_DD(2,ix)+2._rp*cx*d_cx*DD(3)+cx2*&
    d_DD(3,ix)
    h(6,7,ix) = (d_cx*cy*cz2+cx*d_cy*cz2+2._rp*cx*cy*cz*d_cz)*&
    (3._rp*DD(1)-4._rp*DD(2)+DD(3))+cx*cy*cz2*(3._rp*d_DD(1,ix)-4._rp*&
    d_DD(2,ix)+d_DD(3,ix))+(d_cy*cx+cy*d_cx)*(DD(2)-DD(3))+cy*cx*&
    (d_DD(2,ix)-d_DD(3,ix))
    h(6,8,ix) = (d_cy*cz*(cx2-cy2)+cy*d_cz*(cx2-cy2)+2._rp*cy*cz*&
    (cx*d_cx-cy*d_cy))*(1.5*DD(1)-2._rp*DD(2)+0.5_rp*DD(3))+cy*cz*&
    (cx2-cy2)*(1.5*d_DD(1,ix)-2._rp*d_DD(2,ix)+0.5_rp*d_DD(3,ix))-&
    (d_cy*cz+cy*d_cz)*(DD(2)-DD(3))-cy*cz*(d_DD(2,ix)-d_DD(3,ix))
    h(6,9,ix) = sqrt_three*(d_cy*cz*(cx2+cy2)+cy*d_cz*(cx2+cy2)+2._rp*cy*cz*&
    (cx*d_cx+cy*d_cy))*(-0.5_rp*DD(1)+DD(2)-0.5_rp*DD(3))+sqrt_three*cy*cz*&
    (cx2+cy2)*(-0.5_rp*d_DD(1,ix)+d_DD(2,ix)-0.5_rp*d_DD(3,ix))+sqrt_three*&
    (d_cy*cz3+3._rp*cy*cz2*d_cz)*(DD(1)-DD(2))+sqrt_three*cy*cz*cz2*&
    (d_DD(1,ix)-d_DD(2,ix))
    h(7,7,ix) = 2._rp*cx*cz*(cz*d_cx+d_cz*cx)*(3._rp*DD(1)-4._rp*DD(2)+&
    DD(3))+cz2*cx2*(3._rp*d_DD(1,ix)-4._rp*d_DD(2,ix)+d_DD(3,ix))+2._rp*&
    (cz*d_cz+cx*d_cx)*DD(2)+(cz2+cx2)*d_DD(2,ix)+2._rp*cy*d_cy*DD(3)+&
    cy2*d_DD(3,ix)
    h(7,8,ix) = ((d_cz*cx+cz*d_cx)*(cx2-cy2)+cz*cx*2._rp*(cx*d_cx-cy*d_cy))*&
    (1.5*DD(1)-2._rp*DD(2)+0.5_rp*DD(3))+cz*cx*(cx2-cy2)*(1.5*d_DD(1,ix)-&
    2._rp*d_DD(2,ix)+0.5_rp*d_DD(3,ix))+(d_cz*cx+cz*d_cx)*(DD(2)-DD(3))+cz*cx*&
    (d_DD(2,ix)-d_DD(3,ix))
    h(7,9,ix) = sqrt_three*(d_cx*cz3+3._rp*cx*cz2*d_cz)*(DD(1)-DD(2))&
    +sqrt_three*cx*cz3*(d_DD(1,ix)-d_DD(2,ix))+sqrt_three*((d_cx*cz+cx*d_cz)*&
    (cx2+cy2)+2._rp*cx*cz*(cx*d_cx+cy*d_cy))*(-0.5_rp*DD(1)+DD(2)-0.5_rp*&
    DD(3))+sqrt_three*cx*cz*(cx2+cy2)*(-0.5_rp*d_DD(1,ix)+d_DD(2,ix)-0.5_rp*&
    d_DD(3,ix))
    h(8,8,ix) = 4._rp*(cx2-cy2)*(cx*d_cx-cy*d_cy)*(0.75_rp*DD(1)-DD(2)+0.25*&
    DD(3))+(cx2-cy2)**2._rp*(0.75_rp*d_DD(1,ix)-d_DD(2,ix)+0.25*d_DD(3,ix))+&
    2._rp*(cx*d_cx+cy*d_cy)*DD(2)+(cx2+cy2)*d_DD(2,ix)+2._rp*cz*d_cz*DD(3)&
    +cz2*d_DD(3,ix)
    h(8,9,ix) = 2._rp*sqrt_three*((cx*d_cx-cy*d_cy)*cz2+(cx2-cy2)*cz*d_cz)*&
    (0.5_rp*DD(1)-DD(2)+0.25_rp*DD(3))+sqrt_three*(cx2-cy2)*cz2*(0.5_rp*&
    d_DD(1,ix)-d_DD(2,ix)+0.25_rp*d_DD(3,ix))+0.5_rp*sqrt_three*&
    (cx*d_cx-cy*d_cy)*(-DD(1)*(cx2+cy2)+DD(3))+0.25_rp*sqrt_three*(cx2-cy2)*&
    (-d_DD(1,ix)*(cx2+cy2)-DD(1)*2._rp*(cx*d_cx+cy*d_cy)+d_DD(3,ix))
    h(9,9,ix) = 4._rp*(cz2-0.5_rp*(cx2+cy2))*(cz*d_cz-0.5_rp*&
    (cx*d_cx+cy*d_cy))*DD(1)+(cz2-0.5_rp*(cx2+cy2))**2._rp*d_DD(1,ix)+&
    6.0_rp*(cz*d_cz*(cx2+cy2)+cz2*(cx*d_cx+cy*d_cy))*DD(2)+3._rp*cz2*&
    (cx2+cy2)*d_DD(2,ix)+3._rp*(cx2+cy2)*(cx*d_cx+cy*d_cy)*DD(3)+0.75_rp*&
    (cx2+cy2)**2._rp*d_DD(3,ix)

    do k=2,5
      km=k-1
      do l=1,km
        h(k+4,l+4,ix)=h(l+4,k+4,ix)
      end do
    end do
  end subroutine d_slater_koster

  !> Read object in text format from file (default: 'in_atom_tb.txt')
  subroutine read_txt(obj,file)
    class(atom_tb),intent(inout) :: obj
    character(len=*),intent(in),optional :: file
    character(len=:),allocatable :: file_rt

    if (present(file)) then
      file_rt = trim(file)
    else
      file_rt = 'in_atom_tb.txt'
    endif

    ! Parent type procedure
    call obj%atom%read_txt(file_rt)
    ! Derived type procedure
    call obj%calculate_neighbours(obj%e_tb%r_c_max,obj%e_tb%tb_type)
    !deallocate(file_rt)
  end subroutine read_txt

end module atom_tb_mod
