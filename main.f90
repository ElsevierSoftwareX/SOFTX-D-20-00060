! EQUIVALENT POLYNOMIAL LIBRARY V 1.1
! BY GIULIO VENTURA
! giulio.ventura@polito.it
! www.equivalent-polynomials.net
!
! WHEN USING THIS LIBRARY PLEASE ALWAYS CITE
!
! Ventura G., On the elimination of quadrature subcells for discontinuous functions in the extended finite-element method, 
! International Journal for Numerical Methods in Engineering 66 (2006) 761–795.
!
! Ventura G., Benvenuti E., Equivalent polynomials for quadrature in heaviside function enriched elements, 
! International Journal for Numerical Methods in Engineering 102 (2015) 688–710.
!
! THIS IS AN EXAMPLE USAGE FILE WHERE EQUIVALENT POLYNOMIALS ARE USED TO COMPUTE THE VOLUMES AND MOMENTS OF
! INERTIA OF THE POSITIVE AND NEGATIVE ELEMENT SUBDOMAIN
!
! See the Library Application Note and the above papers for details

program eqpol_test

  implicit none
  integer etype,ngp
  double precision a,b,c,d,t,x,y,z
  double precision, allocatable :: w(:),eqcv(:),gp(:),gw(:)

  ! note: eqcv will contain the coefficients of the equivalent polynomial
  ! note: w is a vector of coefficients with only one component = 1 and the others = 0. 
  !       when HEqPol is called with w as coefficients instead of eqcv, it will give the value of
  !       the i-th monomial basis function at the evaluation point
  ! note: regardless of the particular case, the Gauss quadrature rule is set for the worst case (etype=31)

  integer i,j,k,fun

  interface

    double precision function HEqPol(x,y,z,eqcv,etype)
      implicit none
      integer etype
      double precision x,y,z,eqcv(:)
    end function HEqPol

    subroutine Heqpol_coefficients(a,b,c,d,eqcv,etype)
      implicit none
      integer etype
      double precision a,b,c,d,eqcv(:)
    end subroutine Heqpol_coefficients

  end interface

  write(*,*)
  write(*,'(1X,A)') 'EQUIVALENT POLYNOMIAL LIBRARY TEST PROGRAM'
  write(*,*)

  write(*,'(1X,A)',advance='no') 'etype (20,21,30,31) : '
  read (*,*) etype
  if (etype.lt.30) then
    ! 2D case
    write(*,'(1X,A)',advance='no') 'a,b,c : ' 
    read(*,*) a,b,c
    if (etype==20) ngp=1
    if (etype==21) ngp=3
  else
    ! 3D case
    write(*,'(1X,A)',advance='no') 'a,b,c,d : ' 
    read(*,*) a,b,c,d
    if (etype==30) ngp=1
    if (etype==31) ngp=3
  end if

  allocate(gp(ngp),gw(ngp))
  select case (etype)
  case (20)
    allocate(w(1),eqcv(1))
    gp=(/0.333333333333/); gw=(/0.5/)
  case (30)
    allocate(w(1),eqcv(1))
    gp=(/0.25/); gw=(/0.16666666666666666/)
  case (21)
    allocate(w(6),eqcv(6))
    gp=(/-0.77459667,0.,0.77459667/)
    gw=(/0.55555555,0.88888889,0.55555555/)
  case (31)
    allocate(w(23),eqcv(23))
    gp=(/-0.77459667,0.,0.77459667/)
    gw=(/0.55555555,0.88888889,0.55555555/)
  end select

  if (etype.lt.30) then

    ! 2D case
    z=0; d=0
    ! positive part
    call Heqpol_coefficients(a,b,c,d,eqcv,etype)
    write(*,*) 'Gauss integration result of the equivalent polynomial times basis functions - positive part'
    do fun=1,size(eqcv)
      t=0; w=0; w(fun)=1
      if (etype==20 .or. etype==30) then
        do i=1,size(gp)
          x=gp(i); y=gp(i)
          !   equivalent polynomial   Gauss weights fun-th monomial basis function
          t=t+HEqPol(x,y,z,eqcv,etype) *    gw(i) * HEqPol(x,y,z,w,etype)
        end do
      else if (etype==21 .or. etype==31) then
        do i=1,size(gp)
          x=gp(i)
          do j=1,size(gp)
            y=gp(j)
            !   equivalent polynomial   Gauss weights fun-th monomial basis function
            t=t+HEqPol(x,y,z,eqcv,etype)*gw(i)*gw(j)*HEqPol(x,y,z,w,etype)
          end do
        end do
      end if
      write(*,'(I2,5X,F10.7)') fun,t
    end do
    ! negative part
    call Heqpol_coefficients(-a,-b,-c,-d,eqcv,etype)
    write(*,*) 'Gauss integration result of the equivalent polynomial times basis functions - negative part'
    do fun=1,size(eqcv)
      t=0; w=0; w(fun)=1
      if (etype==20 .or. etype==30) then
        do i=1,size(gp)
          x=gp(i); y=gp(i)
          !   equivalent polynomial   Gauss weights fun-th monomial basis function
          t=t+HEqPol(x,y,z,eqcv,etype) *    gw(i) * HEqPol(x,y,z,w,etype)
        end do
      else if (etype==21 .or. etype==31) then
        do i=1,size(gp)
          x=gp(i)
          do j=1,size(gp)
            y=gp(j)
            !   equivalent polynomial   Gauss weights fun-th monomial basis function
            t=t+HEqPol(x,y,z,eqcv,etype)*gw(i)*gw(j)*HEqPol(x,y,z,w,etype)
          end do
        end do
      end if
      write(*,'(I2,5X,F10.7)') fun,t
    end do

  else
    ! 3D case
    ! positive part
    call Heqpol_coefficients(a,b,c,d,eqcv,etype)
    write(*,*) 'Gauss integration result of the equivalent polynomial times basis functions - positive part'
    do fun=1,size(eqcv)
      t=0; w=0; w(fun)=1
      if (etype==20 .or. etype==30) then
        do i=1,size(gp)
          x=gp(i); y=gp(i); z=gp(i)
          !   equivalent polynomial    Gauss weights     fun-th monomial basis function
          t=t+HEqPol(x,y,z,eqcv,etype) *    gw(i) *      HEqPol(x,y,z,w,etype)
        end do
      else if (etype==21 .or. etype==31) then
        do i=1,size(gp)
          x=gp(i)
          do j=1,size(gp)
            y=gp(j)
            do k=1,size(gp)
              z=gp(k)
              !   equivalent polynomial    Gauss weights     fun-th monomial basis function
              t=t+HEqPol(x,y,z,eqcv,etype)*gw(i)*gw(j)*gw(k)*HEqPol(x,y,z,w,etype)
            end do
          end do
        end do
      end if
      write(*,'(I2,5X,F10.7)') fun,t
    end do
    ! negative part
    call Heqpol_coefficients(-a,-b,-c,-d,eqcv,etype)
    write(*,*) 'Gauss integration result of the equivalent polynomial times basis functions - negative part'
    do fun=1,size(eqcv)
      t=0; w=0; w(fun)=1
      if (etype==20 .or. etype==30) then
        do i=1,size(gp)
          x=gp(i); y=gp(i); z=gp(i)
          !   equivalent polynomial    Gauss weights     fun-th monomial basis function
          t=t+HEqPol(x,y,z,eqcv,etype) *    gw(i) *      HEqPol(x,y,z,w,etype)
        end do
      else if (etype==21 .or. etype==31) then
        do i=1,size(gp)
          x=gp(i)
          do j=1,size(gp)
            y=gp(j)
            do k=1,size(gp)
              z=gp(k)
              !   equivalent polynomial    Gauss weights     fun-th monomial basis function
              t=t+HEqPol(x,y,z,eqcv,etype)*gw(i)*gw(j)*gw(k)*HEqPol(x,y,z,w,etype)
            end do
          end do
        end do
      end if
      write(*,'(I2,5X,F10.7)') fun,t
    end do

  end if
  
  write(*,*)
  
end program eqpol_test
