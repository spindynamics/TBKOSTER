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
!  math.f90
!  TBKOSTER
module math_mod
#if defined(LAPACK95_FOUND)
  use lapack95, only: getrf, getri
#endif
  use precision_mod, only: ip, rp
  implicit none

  !> type for sorting routine
  type, public :: group
    integer :: order        ! original order of unsorted data
    real(rp) :: value       ! values to be sorted by
  end type group

  ! To do: precise the quantity (charge, energy, etc)
  real(rp), parameter :: epsilon = 1.0E-10_rp

  ! Irrational
  !> \f$ \sqrt{2} \f$
  real(rp), parameter :: sqrt_two = sqrt(2.0_rp)
  !> \f$ 1/\sqrt{2} \f$
  real(rp), parameter :: one_over_sqrt_two = 1.0_rp/sqrt_two
  !> \f$ \sqrt{3} \f$
  real(rp), parameter :: sqrt_three = sqrt(3.0_rp)
  !> \f$ 1/3 \f$
  real(rp), parameter :: one_third = 0.333333333333333333333_rp
  !> \f$ 2/3 \f$
  real(rp), parameter :: two_third = 0.666666666666666666666_rp

  ! Transcendental
  !> \f$ \pi \f$
  real(rp), parameter :: pi = 3.14159265358979323846_rp
  !> \f$ 2\times\pi \f$
  real(rp), parameter :: two_pi  = 2.0_rp*pi
  !> \f$ 4\times\pi \f$
  real(rp), parameter :: four_pi = 4.0_rp*pi
  !> \f$ \sqrt{\pi} \f$
  real(rp), parameter :: sqrt_pi = sqrt(pi)
  !> \f$ 1/\sqrt{\pi} \f$
  real(rp), parameter :: one_over_sqrt_pi = 1.0_rp/sqrt_pi
  !> \f$ \sqrt{2\pi} \f$
  real(rp), parameter :: sqrt_two_pi = sqrt_two*sqrt_pi
  !> \f$ 1/\sqrt{2\pi} \f$
  real(rp), parameter :: one_over_sqrt_two_pi = 1.0_rp/sqrt_two_pi
  !> \f$ \pi/180 \f$
  real(rp), parameter :: deg2rad = pi/180.0_rp
  !> \f$ 180/\pi \f$
  real(rp), parameter :: rad2deg = 180.0_rp/pi

  ! Complex
  !> Imaginary unit \f$ \mathrm{i} \f$
  complex(rp), parameter :: i_unit = cmplx(0.0_rp,1.0_rp,kind=rp)
  !> Identity matrix \f$ \mathbf{I} \f$ (rank 2)
  complex(rp), dimension(2,2), parameter :: sigma_0 &
    = reshape((/1.0_rp,0.0_rp,0.0_rp,1.0_rp/), (/2,2/))
  ! Pauli matrices
  !> Pauli matrix \f$ \sigma_x \f$
  complex(rp), dimension(2,2), parameter :: sigma_x &
    = reshape((/0.0_rp,1.0_rp,1.0_rp,0.0_rp/), (/2,2/))
  !> Pauli matrix \f$ \sigma_y \f$
  complex(rp), dimension(2,2), parameter :: sigma_y &
    = reshape((/0.0_rp,1.0_rp,-1.0_rp,0.0_rp/), (/2,2/))*i_unit
  !> Pauli matrix \f$ \sigma_z \f$
  complex(rp), dimension(2,2), parameter :: sigma_z &
    = reshape((/1.0_rp,0.0_rp,0.0_rp,-1.0_rp/), (/2,2/))
  complex(rp), dimension(9,9), parameter :: L_x &
    = reshape((/0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,-1.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,-1.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,-1.0_rp,-sqrt_three, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &    
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,sqrt_three,0.0_rp,0.0_rp,0.0_rp/),&                   
                (/9,9/))*i_unit
  complex(rp), dimension(9,9), parameter :: L_y &
    = reshape((/0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,-1.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,-1.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,-1.0_rp,sqrt_three, &    
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,-sqrt_three,0.0_rp,0.0_rp/),&                   
                (/9,9/))*i_unit
  complex(rp), dimension(9,9), parameter :: L_z &
    = reshape((/0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,-1.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,1.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,2.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,1.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,-1.0_rp,0.0_rp,0.0_rp,0.0_rp, &    
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,-2.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, &
                0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp, 0.0_rp,0.0_rp,0.0_rp/),&                   
                (/9,9/))*i_unit

contains
  !> Create vector of all zeros
  function zeros1(n) result(v)
    integer,intent(in) :: n
    real(rp),dimension(n) :: v

    v = 0.0_rp
  end function zeros1

  !> Create matrix of all zeros
  function zeros2(n1,n2) result(m)
    integer,intent(in) :: n1,n2
    real(rp),dimension(n1,n2) :: m

    m = 0.0_rp
  end function zeros2

  !> Create vector of all ones
  function ones1(n) result(v)
    integer,intent(in) :: n
    real(rp),dimension(n) :: v

    v = 1.0_rp
  end function ones1

  !> Create matrix of all ones
  function ones2(n1,n2) result(m)
    integer,intent(in) :: n1,n2
    real(rp),dimension(n1,n2) :: m

    m = 1.0_rp
  end function ones2

  !> Create diagonal matrix
  function diag1(v) result(m)
    real(rp),dimension(:),intent(in) :: v
    real(rp),dimension(size(v),size(v)) :: m
    integer :: i,n

    m = 0.0_rp
    n = size(v)

    do i=1,n
      m(i,i) = v(i)
    end do
  end function diag1

  !> Get diagonal element of matrix
  function diag2(m) result(v)
    real(rp),dimension(:,:),intent(in) :: m
    real(rp),dimension(min(size(m,1),size(m,2))) :: v
    integer :: i,n

    n = min(size(m,1),size(m,2))

    do i=1,n
      v(i) = m(i,i)
    end do
  end function diag2

  !> Identity matrix
  function eye(n) result(m)
    integer,intent(in) :: n
    real(rp),dimension(n,n) :: m
    integer :: i

    m = 0.0_rp

    do i=1,n
      m(i,i) = 1.0_rp
    end do
  end function eye

  !> Test if rank 4 matrix is Hermitian
  function is_hermitian_r4(m) result(l)
    complex(rp),dimension(:,:,:,:),intent(in) :: m
    logical :: l
    integer :: ia,io1,io2,is1,is2
    integer,dimension(2,2),parameter :: iss2is = reshape((/1,4,3,2/),(/2,2/))

    l = .true.
    do ia=1,size(m,1)
      do io1=1,size(m,2)
        do io2=1,size(m,3)
          do is1=1,2
            do is2=1,2
              l = l .and. abs(m(ia,io1,io2,iss2is(is1,is2)) &
               -conjg(m(ia,io2,io1,iss2is(is2,is1))))<epsilon
              if(.not. l) return
            end do !end is1
          end do !end is2
        end do !end io1
      end do !end io2
    end do !end ia
  end function is_hermitian_r4

  !> Vector cross product
  function cross_product(v1,v2) result(v3)
    real(rp),dimension(3),intent(in) :: v1,v2
    real(rp),dimension(3) :: v3

    v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
    v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
    v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
  end function cross_product

  !> Normalize a vector
  function normalize(v1) result(v2)
    real(rp), dimension(3), intent(in) :: v1
    real(rp), dimension(3) :: v2
    real(rp) :: norm
    norm=norm2(v1)
    v2(:) = v1(:)/norm
  end function normalize

  !> Matrix determinant
  function determinant(m) result(d)
    real(rp),dimension(3,3),intent(in) :: m
    real(rp) :: d
    real(rp), dimension(3) :: v1,v2,v3

    v1(:) = m(1,:)
    v2(:) = m(2,:)
    v3(:) = m(3,:)

    d = dot_product(v1,cross_product(v2,v3))

  end function determinant

  !> Matrix inverse
  function inverse(m1) result(m2)
    real(rp),dimension(:,:),intent(in) :: m1
    real(rp),dimension(size(m1,1),size(m1,2)) :: m2
    integer,dimension(min(size(m1,1),size(m1,2))) :: ipiv
#if !defined(LAPACK95_FOUND)
    integer :: m
    integer :: lda
    real(rp),dimension(max(1,size(m1,1))) :: work
    integer :: lwork
    integer :: info
#endif
#if defined(LAPACK95_FOUND)
    m2 = m1
    call getrf(m2,ipiv)
    call getri(m2,ipiv)
#else
    m = size(m1,1)
    lda = max(1,m)
    lwork = max(1,m)
    m2 = m1
    call dgetrf(m,m,m2,lda,ipiv,info)
    call dgetri(m,m2,lda,ipiv,work,lwork,info)
#endif
  end function inverse

  !> Matrix inverse of a 3x3 matrix
  function inverse_3x3(m1) result(m2)
    real(rp),dimension(3,3),intent(in) :: m1
    real(rp),dimension(3,3) :: m2
    real(rp) :: det

    det = determinant(m1)
    if (det >= epsilon) then

    m2(1,1) =  (m1(2,2) * m1(3,3) - m1(2,3) * m1(3,2)) / det;
    m2(1,2) = -(m1(1,2) * m1(3,3) - m1(1,3) * m1(3,2)) / det;
    m2(1,3) =  (m1(1,2) * m1(2,3) - m1(1,3) * m1(2,2)) / det;
    m2(2,1) = -(m1(2,1) * m1(3,3) - m1(2,3) * m1(3,1)) / det;
    m2(2,2) =  (m1(1,1) * m1(3,3) - m1(1,3) * m1(3,1)) / det;
    m2(2,3) = -(m1(1,1) * m1(2,3) - m1(1,3) * m1(2,1)) / det;
    m2(3,1) =  (m1(2,1) * m1(3,2) - m1(2,2) * m1(3,1)) / det;
    m2(3,2) = -(m1(1,1) * m1(3,2) - m1(1,2) * m1(3,1)) / det;
    m2(3,3) =  (m1(1,1) * m1(2,2) - m1(1,2) * m1(2,1)) / det;
    else
      write(*,*) 'math%inverse_3x3(): Warning - determinant is nearely zero'
      stop
    end if
     
  end function inverse_3x3

  !> Unique values c in an integer array a, index ia is such that c = a(ia) and
  !> index ic is such that a = c(ic)
  subroutine unique_int(a,c,ia,ic)
    integer,dimension(:),intent(in)  :: a
    integer,dimension(:),allocatable,intent(out) :: c,ia
    integer,dimension(size(a)),intent(out) :: ic
    integer :: sa,sc,ia1,ia2,ic1
    integer,dimension(size(a)) :: index_a
    logical,dimension(size(a)) :: mask

    sa = size(a)
    index_a = (/(ia1, ia1=1, sa)/)
    mask = .false.
    do ia1=1,sa-1
      if(.not. mask(ia1)) then
        do ia2=ia1+1,sa
          mask(ia2) = a(ia1) == a(ia2)
          if(mask(ia2)) index_a(ia2) = index_a(ia1)
        end do
      end if
    end do
    sc = sa - count(mask)
    allocate(c(sc),ia(sc))
    c = pack(a, .not. mask)
    ia = pack(index_a, .not. mask)
    do ia1=1,sa
      do ic1=1,sc
        if(index_a(ia1) == ia(ic1)) then
          ic(ia1) = ic1
          exit
        end if
      end do
    end do
  end subroutine unique_int

  ! Unique values c in a real array a down to precision epsilon, index ia is
  ! such that c = a(ia) andindex ic is such that a = c(ic)
  !subroutine unique_real(a,c,ia,ic)

  !> Returns a spherical vector \f$ v_{\mathrm{sph}} = (v_r,v_{\theta},v_{\phi})
  !> \f$ from a cartesian vector \f$ v_{\mathrm{cart}} = (v_x,v_y,v_z) \f$ where
  !> \f$ 0 \leq v_{\theta} \leq \pi \f$ and
  !> \f$ -\pi \leq v_{\phi} < \pi \f$
  function cart2sph(v_cart) result(v_sph)
    real(rp),dimension(3),intent(in) :: v_cart
    real(rp),dimension(3) :: v_sph

    v_sph(1)=norm2(v_cart)

    if(v_sph(1)<epsilon) then
      v_sph(2)=0.0_rp
    else
      v_sph(2)=acos(v_cart(3)/v_sph(1))
    end if

    if(sqrt(v_cart(1)*v_cart(1)+v_cart(2)*v_cart(2))<epsilon) then
      v_sph(3)=0.0_rp
    else
      if(abs(v_cart(2))>epsilon) then
        if(v_cart(2)>0.0_rp) then
          v_sph(3)= acos(v_cart(1)/sqrt(v_cart(1)*v_cart(1)+v_cart(2)*v_cart(2)))
        else
          v_sph(3)=-acos(v_cart(1)/sqrt(v_cart(1)*v_cart(1)+v_cart(2)*v_cart(2)))
        end if
      elseif(abs(v_cart(2))<epsilon) then
        if(v_cart(1)>0.0_rp) then
          v_sph(3)=0.0_rp
        else
          v_sph(3)=-pi
        end if
      end if
    end if
   ! write(6,*) "v_cart=",v_cart(:)
   ! write(6,*) "v_sph=",v_sph(:)
  end function cart2sph

  !> Returns a cartesian vector \f$ v_{\mathrm{cart}} = (v_x,v_y,v_z) \f$ from a
  !> spherical vector \f$ v_{\mathrm{sph}} = (v_r,v_{\theta},v_{\phi}) \f$
  function sph2cart(v_sph) result(v_cart)
    real(rp),dimension(3),intent(in) :: v_sph
    real(rp),dimension(3) :: v_cart

    v_cart(1)=v_sph(1)*sin(v_sph(2))*cos(v_sph(3))
    v_cart(2)=v_sph(1)*sin(v_sph(2))*sin(v_sph(3))
    v_cart(3)=v_sph(1)*cos(v_sph(2))
  end function sph2cart

  subroutine nm2rho(n,m_cart,rho)
    real(rp),intent(in) :: n
    real(rp),dimension(3),intent(in) :: m_cart
    complex(rp),dimension(2,2),intent(out) :: rho

    rho(1,1) = (n+m_cart(3))/2
    rho(2,2) = (n-m_cart(3))/2
    rho(1,2) = (m_cart(1)-i_unit*m_cart(2))/2
    rho(2,1) = (m_cart(1)+i_unit*m_cart(2))/2
  end subroutine nm2rho

  subroutine rho2nm(rho,n,m_cart)
    complex(rp),dimension(2,2),intent(in) :: rho
    real(rp),intent(out) :: n
    real(rp),dimension(3),intent(out) :: m_cart
    integer :: ispin,jspin

    do ispin=1,2
      do jspin=1,2
        if((abs(rho(ispin,jspin)-conjg(rho(jspin,ispin))))>epsilon) then
          write(*,*) 'math%rho2nm(): Warning - non Hermitian density'
        end if
      end do
    end do

    n = real(rho(1,1)+rho(2,2))
    m_cart(1) = real(rho(1,2)+rho(2,1))
    !m_cart(2) = i_unit*(rho(1,2)-rho(2,1))
    m_cart(2) = aimag(rho(2,1)-rho(1,2))
    m_cart(3) = real(rho(1,1)-rho(2,2))
  end subroutine rho2nm

  !> Compute Fermi(E,beta) occupation function 
  function fermi_function(E,beta) result(f)
    real(rp), intent(in) :: E,beta
    real(rp) :: f
    real(rp), parameter :: maxarg = 100._rp
    real(rp) :: x

    x=beta*E

    if(x<-maxarg) then
      f=1.0_rp
    elseif(x>maxarg) then
      f=0.0_rp
    else
      f=1.0_rp/(1.0_rp + exp(x))
    end if
  end function fermi_function

  !> Compute dFermi(E,beta)/dE
  function fermi_function_derivative(E,beta) result(f)
    real(rp), intent(in) :: E,beta
    real(rp) :: f
    real(rp), parameter :: maxarg = 100._rp
    real(rp) :: x

    x=beta*E

    if(x<-maxarg) then
      f=0.0_rp
    elseif(x>maxarg) then
      f=0.0_rp
    else
      f=-beta*exp(x)/(1.0_rp+exp(x))**2
    end if
  end function fermi_function_derivative

  !> fill an array arr and sort it to get the indx of sorted values
  subroutine indexx(n,arr,indx)
    ! INPUT
    integer,intent(in) :: n
    real(rp),dimension(0:n),intent(in) :: arr
    ! OUTPUT
    integer,dimension(0:n),intent(out) :: indx
    ! LOCAL
    integer,parameter :: m=7, nstack=50
    integer,dimension(:),allocatable :: istack
    integer :: i,indxt,ir,itemp,j,jstack,k,l
    real(rp) :: a
    if(.not. allocated(istack)) allocate(istack(nstack))
    indx = (/(j, j=0,n)/)
    jstack=0
    l=1
    ir=n
    do
      if(ir-l<m) then
        do j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do i=j-1,1,-1
            if(arr(indx(i))<=a) exit
            indx(i+1)=indx(i)
          end do
          if(i==1 .and. arr(indx(i))>a) i=0
          indx(i+1)=indxt
        end do
        if(jstack==0) return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1))>arr(indx(ir))) then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        end if
        if(arr(indx(l))>arr(indx(ir))) then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        end if
        if(arr(indx(l+1))>arr(indx(l))) then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        end if
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
        do
          i=i+1
          if(arr(indx(i))<a) cycle
          do
            j=j-1
            if(arr(indx(j))<=a) exit
          end do
          if(j<i) exit
          itemp=indx(i)
          indx(i)=indx(j)
          indx(j)=itemp
        end do
        indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack > nstack) then
          write(*,*) "nstack too small in indexx"
          stop
        end if
        if(ir-i+1>=j-l) then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        end if
      end if
    end do
    if(allocated(istack)) deallocate(istack)
  end subroutine indexx

  !> fill an array arr and sort it to get the indx of sorted values
  recursive subroutine QSort(na,a,indx)
 
    ! DUMMY ARGUMENTS
    integer, intent(in) :: nA
    type(group), dimension(nA), intent(in) :: A
    integer,dimension(nA),intent(out) :: indx
 
    ! LOCAL VARIABLES
    integer :: left, right
    real(rp) :: random
    real(rp) :: pivot
    type (group) :: temp
    type (group), dimension(nA) :: B
    integer :: i,marker

    forall(i=1:nA) B(i)=A(i)

    if (nA > 1) then
      call random_number(random)
      pivot = B(int(random*real(nA-1))+1)%value   ! random pivot (not best performance, but avoids worst-case)
      left = 0
      right = nA + 1
      do while (left < right)
        right = right - 1
        do while (B(right)%value > pivot)
          right = right - 1
        end do
        left = left + 1
        do while (B(left)%value < pivot)
          left = left + 1
        end do
        if (left < right) then
          temp = B(left)
          B(left) = B(right)
          B(right) = temp
        end if
      end do
 
      if (left == right) then
        marker = left + 1
      else
        marker = left
      end if
 
      call QSort(marker-1,B(:marker-1),indx)
      call QSort(nA-marker+1,B(marker:),indx)
    end if
    do i=1,nA
      indx(i)=B(i)%order
    end do

  end subroutine QSort

  real(rp) function theta_function(x,smearing)
    !-----------------------------------------------------------------------
    !
    !     This function computes the approximate theta function for the
    !     given order n, at the point x.
    !
    ! --> 'fd': Fermi-Dirac smearing.
    !       1.0/(1.0+exp(-x))
    !
    ! --> 'g':  Gaussian smearing (n=0).
    !
    ! --> 'mp': Methfessel-Paxton smearing (n=1). See PRB 40, 3616 (1989).
    !
    ! --> 'mv': Marzari-Vanderbilt (cold) smearing. See PRL 82, 3296 (1999).
    !       1/2*erf(x-1/sqrt(2)) + 1/sqrt(2*pi)*exp(-(x-1/sqrt(2))**2) + 1/2
    !
    real(rp), intent(in) :: x
    character(len=*), intent(in) :: smearing
    ! output: the value of the function
    ! input: the argument of the function
    ! input: the order of the function
    !
    !    the local variables
    !
    real(rp) :: a, hp, arg, hd, xp
    ! the coefficient a_n
    ! the hermitean function
    ! the argument of the exponential
    ! the hermitean function
    ! auxiliary variable (cold smearing)
    integer :: i, ni,ngauss
    ! counter on the n indices
    ! counter on 2n
    real(rp), parameter :: maxarg = 100._rp
    ! maximum value for the argument of the exponential

    ! Fermi-Dirac smearing
    if(smearing=='fd') then
      if(x<-maxarg) then
        theta_function = 0._rp
      elseif(x>maxarg) then
        theta_function = 1._rp
      else
        theta_function = 1.0_rp / (1.0_rp + exp ( - x) )
      end if
      return
    end if
    ! Cold smearing
    if(smearing=='mv') then
      xp = x - one_over_sqrt_two
      arg = min (maxarg, xp**2)
      theta_function = 0.5_rp * erf_qe (xp) + one_over_sqrt_two_pi * exp ( - &
      arg) + 0.5_rp
      return
    end if
    if(smearing=='mp') ngauss=1
    if(smearing=='g')  ngauss=0

    ! Methfessel-Paxton
    theta_function = 0.5_rp * erfc_qe (-x)

    if(ngauss==0) return
    hd = 0._rp
    arg = min (maxarg, x**2)
    hp = exp ( - arg)
    ni = 0
    a = one_over_sqrt_pi
    do i = 1, ngauss
      hd = 2.0_rp * x * hp - 2.0_rp * dble (ni) * hd
      ni = ni + 1
      a = - a / (dble (i) * 4.0_rp)
      theta_function = theta_function - a * hd
      hp = 2.0_rp * x * hd-2.0_rp * dble (ni) * hp
      ni = ni + 1
    end do
    return
  end function theta_function

  !***********************************************************************|
  real(rp) function delta_function(x,smearing)
    !
    !     The derivative of wgauss: an approximation to the delta function
    !
    ! --> 'fd': derivative of the Fermi-Dirac smearing.
    !       0.5/(1.0+cosh(x))
    !
    ! --> 'g':  derivative of the Gaussian smearing (n=0).
    !
    ! --> 'mp': derivative of the Methfessel-Paxton smearing.
    !
    ! --> 'mv': derivative of the Marzari-Vanderbilt (cold) smearing.
    !       1/sqrt(pi)*exp(-(x-1/sqrt(2))**2)*(2-sqrt(2)*x)
    !
    real(rp), intent(in) :: x
    character(len=*),intent(in) :: smearing
    ! output: the value of the function
    ! input: the point where to compute the function

    ! input: the order of the smearing function
    !
    !    here the local variables
    !
    real(rp) :: a, arg, hp, hd
    ! the coefficients a_n
    ! the argument of the exponential
    ! the hermite function
    ! the hermite function

    integer :: i, ni,ngauss
    ! counter on n values
    ! counter on 2n values


    ! Fermi-Dirac smearing
    if(smearing=='fd') then
      if(abs(x) <= 36.0_rp) then
        delta_function = 1.0_rp / (2.0_rp + exp ( - x) + exp ( + x) )
        ! in order to avoid problems for large values of x in the e
      else
        delta_function = 0.0_rp
      end if
      return
    end if
    ! cold smearing  (Marzari-Vanderbilt)
    if(smearing=='mv') then
      arg = min (200.0_rp, (x - one_over_sqrt_two ) **2)
      delta_function = one_over_sqrt_pi * exp ( - arg) * (2.0_rp - sqrt_two * x)
      return
    end if
    ! Methfessel-Paxton
    if(smearing=='mp') ngauss=1
    if(smearing=='g')  ngauss=0

    arg = min (200.0_rp, x**2)
    delta_function = exp ( - arg) * one_over_sqrt_pi
    if(ngauss==0) return
    hd = 0.0_rp
    hp = exp ( - arg)
    ni = 0
    a = one_over_sqrt_pi
    do i = 1, ngauss
      hd = 2.0_rp * x * hp - 2.0_rp * dble (ni) * hd
      ni = ni + 1
      a = - a / (dble (i) * 4.0_rp)
      hp = 2.0_rp * x * hd-2.0_rp * dble (ni) * hp
      ni = ni + 1
      delta_function = delta_function + a * hp
    end do
    return
  end function delta_function

  real(rp) function integrated_delta_function(x,smearing)
    !-----------------------------------------------------------------------
    !
    !    w1gauss(x,n) = \int_{-\infty}^x   y delta(y) dy
    !    where delta(x) is the current approximation for the delta function,
    !    as obtained from w0gauss(x,n)
    !
    ! --> 'fd.: Fermi-Dirac smearing (n=-99). In this case w1gauss corresponds
    !     to the negative of the electronic entropy.
    !
    ! --> 'mp': Methfessel-Paxton smearing (n>=0).
    !
    ! --> 'mv': Marzari-Vanderbilt (cold) smearing (n=-1).
    !       1/sqrt(2*pi)*(x-1/sqrt(2))*exp(-(x-1/sqrt(2))**2)
    !
    real(rp), intent(in) :: x
    character(len=*), intent(in) :: smearing
    ! output: the value of the function
    ! input: the point where to compute the function


    ! input: the order of the smearing function
    !
    !    here the local variables
    !

    real(rp) :: a, hp, arg, hpm1, hd, f, onemf, xp
    ! the coefficients a_n
    ! the hermite function
    ! the argument of the exponential
    ! the hermite function
    ! the hermite function
    ! Fermi-Dirac occupation number
    ! 1 - f
    ! auxiliary variable (cold smearing)

    integer :: i, ni,ngauss
    ! counter on n values
    ! counter on 2n values

    if(smearing=='fd') then
      if(abs (x) <= 36.0_rp) then
        f = 1.0_rp / (1.0_rp + exp ( - x) )
        onemf = 1.0_rp - f
        integrated_delta_function = f * log (f) + onemf * log (onemf)
        ! in order to avoid problems for large values of x
      else
        ! neglect w1gauss when abs(w1gauss) < 1.0d-14
        integrated_delta_function = 0.0_rp
      end if
      return

    end if
    ! Cold smearing
    if(smearing=='mv') then
      xp = x - one_over_sqrt_two
      arg = min (200._rp, xp**2)
      integrated_delta_function = one_over_sqrt_two_pi * xp * exp ( - arg)
      return

    end if
    ! Methfessel-Paxton
    if(smearing=='mp') ngauss=1
    if(smearing=='g')  ngauss=0

    arg = min (200.0_rp, x**2)
    integrated_delta_function = - 0.5_rp * one_over_sqrt_pi * exp ( - arg)
    if(ngauss==0) return
    hd = 0.0_rp
    hp = exp ( - arg)
    ni = 0
    a = one_over_sqrt_pi
    do i = 1, ngauss
      hd = 2.0_rp * x * hp - 2.0_rp * dble (ni) * hd
      ni = ni + 1
      hpm1 = hp
      hp = 2.0_rp * x * hd - 2.0_rp * dble (ni) * hp
      ni = ni + 1
      a = - a / (dble (i) * 4.0_rp)
      integrated_delta_function = integrated_delta_function - a * (0.5_rp * hp + dble (ni) * hpm1)
    end do
    return
  end function integrated_delta_function

  !***********************************************************************|
  !
  ! Copyright (C) 2002-2009 Quantum ESPRESSO group
  ! This file is distributed under the terms of the
  ! GNU General Public License. See the file `License'
  ! in the root directory of the present distribution,
  ! or http://www.gnu.org/copyleft/gpl.txt .
  !
  !---------------------------------------------------------------------
  real(rp) function erf_qe (x)
    !---------------------------------------------------------------------
    !
    !     Error function - computed from the rational approximations of
    !     W. J. Cody, Math. Comp. 22 (1969), pages 631-637.
    !
    !     for abs(x) le 0.47 erf is calculated directly
    !     for abs(x) gt 0.47 erf is calculated via erf(x)=1-erfc(x)
    !
    real(rp), intent(in) :: x
    real(rp) :: x2, p1 (4), q1 (4)
    data p1 / 2.426679552305318E2_rp, 2.197926161829415E1_rp, &
    6.996383488619136_rp,  -3.560984370181538E-2_rp /
    data q1 / 2.150588758698612E2_rp, 9.116490540451490E1_rp, &
    1.508279763040779E1_rp, 1.000000000000000_rp /
    !
    if(abs (x) > 6.0_rp) then
      !
      !  erf(6)=1-10^(-17) cannot be distinguished from 1
      !
      erf_qe = sign (1.0_rp, x)
    else
      if(abs(x) <= 0.47_rp) then
        x2 = x**2
        erf_qe=x *(p1 (1) + x2 * (p1 (2) + x2 * (p1 (3) + x2 * p1 (4) ) ) ) &
        / (q1 (1) + x2 * (q1 (2) + x2 * (q1 (3) + x2 * q1 (4) ) ) )
      else
        erf_qe = 1.0_rp - erfc_qe (x)
      end if
    end if
    !
    return
  end function erf_qe

  real(rp) function erfc_qe (x)
    !---------------------------------------------------------------------
    !
    !     Complementary error function
    !     erfc(x) = 1-erf(x)  - See comments in erf
    !
    real(rp),intent(in) :: x
    real(rp) :: ax, x2, xm2, p2 (8), q2 (8), p3 (5), q3 (5), pim1
    data p2 / 3.004592610201616E2_rp,  4.519189537118719E2_rp, &
    3.393208167343437E2_rp,  1.529892850469404E2_rp, &
    4.316222722205674E1_rp,  7.211758250883094_rp,   &
    5.641955174789740E-1_rp,-1.368648573827167E-7_rp /
    data q2 / 3.004592609569833E2_rp,  7.909509253278980E2_rp, &
    9.313540948506096E2_rp,  6.389802644656312E2_rp, &
    2.775854447439876E2_rp,  7.700015293522947E1_rp, &
    1.278272731962942E1_rp,  1.000000000000000_rp /
    data p3 /-2.996107077035422E-3_rp,-4.947309106232507E-2_rp, &
    -2.269565935396869E-1_rp,-2.786613086096478E-1_rp, &
    -2.231924597341847E-2_rp /
    data q3 / 1.062092305284679E-2_rp, 1.913089261078298E-1_rp, &
    1.051675107067932_rp,    1.987332018171353_rp,    &
    1.000000000000000_rp /

    data pim1 / 0.56418958354775629_rp /
    !        ( pim1= sqrt(1/pi) )
    ax = abs (x)
    if(ax > 26.0_rp) then
      !
      !  erfc(26.0)=10^(-296); erfc( 9.0)=10^(-37);
      !
      erfc_qe = 0.0_rp
    elseif(ax > 4.0_rp) then
      x2 = x**2
      xm2 = (1.0_rp / ax) **2
      erfc_qe = (1.0_rp / ax) * exp ( - x2) * (pim1 + xm2 * (p3 (1) &
      + xm2 * (p3 (2) + xm2 * (p3 (3) + xm2 * (p3 (4) + xm2 * p3 (5) &
      ) ) ) ) / (q3 (1) + xm2 * (q3 (2) + xm2 * (q3 (3) + xm2 * &
      (q3 (4) + xm2 * q3 (5) ) ) ) ) )
    elseif(ax > 0.47_rp) then
      x2 = x**2
      erfc_qe = exp ( - x2) * (p2 (1) + ax * (p2 (2) + ax * (p2 (3) &
      + ax * (p2 (4) + ax * (p2 (5) + ax * (p2 (6) + ax * (p2 (7) &
      + ax * p2 (8) ) ) ) ) ) ) ) / (q2 (1) + ax * (q2 (2) + ax * &
      (q2 (3) + ax * (q2 (4) + ax * (q2 (5) + ax * (q2 (6) + ax * &
      (q2 (7) + ax * q2 (8) ) ) ) ) ) ) )
    else
      erfc_qe = 1.0_rp - erf_qe (ax)
    end if
    !
    ! erf(-x)=-erf(x)  =>  erfc(-x) = 2-erfc(x)
    !
    if(x < 0.0_rp) erfc_qe = 2.0_rp - erfc_qe
    !
    return
  end function erfc_qe
end module math_mod
