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
! USAGE
! Call first Heqpol_coefficients for a given element and discontinuity plane
! Then, at each Gauss point of coords x,y,z, evaluate HEqPol, passing the previously computed eqcv
! See the Library example file main.f90 for usage hints.

double precision function HEqPol(x,y,z,eqcv,etype)

  ! On input
  !
  ! x,y,z coordinate of the evaluation point in the PARENT element reference system
  ! x,y,z in [0,1] for triangular and tetrahedral elements
  ! x,y,z in [-1,1] for quad and hexas
  ! the z coordinate is unused for 2D elements
  !
  ! eqcv vector of polynomial coefficients as computed by eqpol_hex_coefficients
  !
  ! etype the element type. It is:
  ! 20  linear triangle, required length for eqcv is 1
  ! 21  linear quadrilateral, required length for eqcv is 6
  ! 30  linear tetrahedron
  ! 31  linear hexahedron, required length for eqcv is 23
  !
  ! On output
  ! the value of the equivalent polynomial at point x,y,z

  implicit none
  integer etype,l
  double precision x,y,z
  double precision eqcv(:)
  double precision v(size(eqcv))

  l = size(eqcv)
  select case (etype)

  case (20)
    if (l.ne.1) stop ' Invalid lenght for vector eqcv'
    v=(/ 1.D0 /)

  case (21)
    if (l.ne.6) stop ' Invalid lenght for vector eqcv'
    v=(/ 1.D0, x, x**2, y, x*y, y**2 /)

  case (30)
    if (l.ne.1) stop ' Invalid lenght for vector eqcv'
    v=(/ 1.D0 /)

  case (31)
    if (l.ne.23) stop ' Invalid lenght for vector eqcv'
    v=(/ 1.D0, x, x**2, y, x*y, x**2*y, y**2, x*y**2, x**2*y**2, z, x*z, x**2*z, y*z, x*y*z, &
    x**2*y*z, y**2*z, x*y**2*z, z**2, x*z**2, x**2*z**2, y*z**2, x*y*z**2, y**2*z**2 /)

  case default
    stop 'Invalid etype specified'

  end select

  HEqPol=dot_product(v,eqcv)

end function HEqPol


subroutine Heqpol_coefficients(a,b,c,d,eqcv,etype)

  ! On input: a,b,c,d, etype
  ! On output: eqcv 
  !
  ! This sub computed the vector of equivalent polynomial coefficients for given
  ! H discontinuity plane coefficients a,b,c,d in the parent element domain
  ! for 3D elements the discontinuity plane has equation a xi + b eta + c zeta + d = 0
  ! for 2D elements the discontinuity plane has equation a xi + b eta + c = 0
  !
  ! On output the vector eqcv of polinomial coefficients w.r.t. the base
  ! defined in the function Hex_HeqPol
  !
  ! etype the element type. It is:
  ! 20  linear triangle, required length for eqcv is 1
  ! 21  linear quadrilateral, required length for eqcv is 6
  ! 30  linear tetrahedron, required length for eqcv is 1
  ! 31  linear hexahedron, required length for eqcv is 23

  implicit none
  integer etype
  double precision a,b,c,d,eqcv(:)
  double precision BV(size(eqcv))
  real, parameter :: tol = 1E-4 ! tolerance for vanishing plane coefficients    

  double precision, parameter :: InvA20(1,1) = reshape( (/ 2.d0 /), shape(InvA20))

  double precision, parameter :: InvA21(6,6) = reshape( &
  (/ 7.d0/8.d0, 0.d0, -15.d0/16.d0, 0.d0, 0.d0, -15.d0/16.d0, &
  0.d0, 3.d0/4.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  -15.d0/16.d0, 0.d0, 45.d0/16.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 3.d0/4.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 9.d0/4.d0, 0.d0, &
  -15.d0/16.d0, 0.d0, 0.d0, 0.d0, 0.d0, 45.d0/16.d0 &
  /),shape(InvA21))

  double precision, parameter :: InvA30(1,1) = reshape( (/ 6.d0 /), shape(InvA30))

  double precision, parameter :: InvA31(23,23) = reshape( &
  (/ 151.d0/128.d0, 0.d0, -(105.d0/64.d0), 0.d0, 0.d0, 0.d0, -(105.d0/64.d0), 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, -(105.d0/64.d0), 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, 21.d0/16.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, &
  0.d0, -(105.d0/64.d0), 0.d0, 315.d0/64.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 21.d0/16.d0, &
  0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 81.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, -(135.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(105.d0/64.d0), 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, &
  315.d0/64.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, -(675.d0/128.d0), 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, &
  0.d0, -(675.d0/128.d0), 0.d0, 2025.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 21.d0/16.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, &
  0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  81.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(135.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 81.d0/32.d0, 0.d0, -(135.d0/32.d0), 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 27.d0/8.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, -(135.d0/32.d0), 0.d0, 405.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(135.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 405.d0/32.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, -(105.d0/64.d0), 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 315.d0/64.d0, 0.d0, -(675.d0/128.d0), 0.d0, &
  0.d0, -(675.d0/128.d0), 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 225.d0/128.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(675.d0/128.d0), 0.d0, 2025.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, -(45.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  135.d0/32.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(135.d0/32.d0), 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 405.d0/32.d0, 0.d0, 225.d0/128.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(675.d0/128.d0), &
  0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, 0.d0, -(675.d0/128.d0), 0.d0, 0.d0, 0.d0, 0.d0, 2025.d0/128.d0 &
  /),shape(InvA31))

  integer i
  double precision a2,a3,a4,a5,a6,b2,b3,b4,b5,b6,c2,c3,c4,c5,c6,d2,d3,d4,d5,d6
  double precision t,am,bm,cm
  double precision abs01,abs02,abs03,abs04,abs05,abs06,abs07,abs08,abs09,abs10,abs11,abs12
  double precision abs13,abs14,abs15,abs16

  ! plane coefficient normalization
  select case (etype)

  case (20,21)
    t=sqrt(a**2+b**2)
    a=a/t; b=b/t; c=c/t
    am=abs(a)
    bm=abs(b)

  case (30,31)
    t=sqrt(a**2+b**2+c**2)
    a=a/t; b=b/t; c=c/t; d=d/t
    am=abs(a)
    bm=abs(b)
    cm=abs(c)

  case default
    stop 'Invalid etype'

  end select


  select case (etype)

  case (20)
    a2=a**2; b2=b**2
    if ( am<bm*tol ) then
      ! a=0, b<>0
      ! b
      BV(1) = (b2-(2*b+c)*abs(c)+(b+c)*abs(b+c))/(4.d0*b2)
    else if ( bm<am*tol ) then
      ! a<>0, b=0
      ! a
      BV(1) = (a2-(2*a+c)*abs(c)+(a+c)*abs(a+c))/(4.d0*a2)
    else if (abs(a-b)<tol) then
      ! a<>0, b<>0, a=b
      ! a=b
      BV(1) = (a2+c*abs(c)+(a-c)*abs(a+c))/(4.d0*a2)
    else
      ! a<>0, b<>0
      ! a b
      BV(1) = ((a-b)*c*abs(c)+b*(a+c)*abs(a+c)+a*((a-b)*b-(b+c)*abs(b+c)))/(4.d0*a*(a-b)*b)
    end if
    eqcv=matmul(InvA20,BV)

  case (21)
    a2=a**2; b2=b**2; c2=c**2
    a3=a**3; b3=b**3; c3=c**3
    if ( am<bm*tol ) then
      ! a=0, b<>0
      ! b
      abs01=abs(b-c)
      abs02=abs(b+c)
      include 'quad_ph_6_b.f90'
    else if ( bm<am*tol ) then
      ! a<>0, b=0
      ! a
      abs01=abs(a-c)
      abs02=abs(a+c)
      include 'quad_ph_6_a.f90'
    else
      ! a<>0, b<>0
      ! a b
      abs01=abs(a-b-c)
      abs02=abs(a+b-c)
      abs03=abs(a-b+c)
      abs04=abs(a+b+c)
      include 'quad_ph_6_ab.f90'
    end if
    eqcv=matmul(InvA21,BV)

  case (30)
    a2=a**2; b2=b**2; c2=c**2; d2=d**2
    a3=a**3; b3=b**3; c3=c**3; d3=d**3
    if ( am<bm*tol .and. am<cm*tol .and. min(bm,cm)/max(bm,cm)>tol ) then
      ! a=0, b<>0, c<>0
      ! b c
      if (abs(b-c)>tol) then
        ! b<>c
        BV(1)=((b-c)*d*(c*d+b*(3*c+d))*abs(d)+c2*(b + d)**2*abs(b+d)+ &
              b2*((b-c)*c2-(c+d)**2*abs(c+d)))/(12.d0*b2*(b-c)*c2)
      else
        ! b=c
        BV(1) = (b3 + d*(3*b + 2*d)*abs(d) + (b2 - b*d - 2*d2)*Abs(b + d))/(12.*b3)
      end if
    else if ( bm<am*tol .and. bm<cm*tol .and. min(am,cm)/max(am,cm)>tol ) then
      ! a<>0, b=0, c<>0
      ! a c
      if (abs(a-c)>tol) then
        ! a<>c
        BV(1)=((a-c)*d*(c*d+a*(3*c+d))*abs(d)+c2*(a+d)**2*abs(a+d)+ &
              a2*((a-c)*c2-(c+d)**2*abs(c+d)))/(12.d0*a2*(a - c)*c2)
      else
        ! a=c
        BV(1) = (a3 + d*(3*a + 2*d)*abs(d) + (a2 - a*d - 2*d2)*abs(a + d))/(12.*a3)
      end if
    else if ( cm<am*tol .and. cm<bm*tol .and. min(am,bm)/max(am,bm)>tol ) then
      ! a<>0, b<>0, c=0
      ! a b
      if (abs(a-b)>tol) then
        ! a<>b
        BV(1)=((a-b)*d*(b*d+a*(3*b+d))*abs(d)+b2*(a+d)**2*abs(a+d)+ &
              a2*((a-b)*b2-(b+d)**2*abs(b+d)))/(12.d0*a2*(a-b)*b2)
      else
        ! a=b
        BV(1) = (a3 + d*(3*a + 2*d)*abs(d) + (a2 - a*d - 2*d2)*abs(a + d))/(12.*a3)
      end if
    else if ( am<cm*tol .and. bm<cm*tol ) then
      ! a=0, b=0, c<>0
      ! c
      BV(1)=(c3-(3*c2+3*c*d+d2)*abs(d)+(c+d)**2*abs(c+d))/(12.d0*c3)
    else if ( bm<am*tol .and. cm<am*tol ) then
      ! a<>0, b=0, c=0
      ! a
      BV(1)=(a3-(3*a2+3*a*d+d2)*abs(d)+(a+d)**2*abs(a+d))/(12.d0*a3)
    else if ( am<bm*tol .and. cm<bm*tol ) then
      ! a=0, b<>0, c=0
      ! b
      BV(1)=(b3-(3*b2+3*b*d+d2)*abs(d)+(b+d)**2*abs(b+d))/(12.d0*b3)
    else
      ! a<>0, b<>0, c<>0
      ! a b c
      if (abs(a-b)>tol .and. abs(a-c)>tol .and. abs(b-c)>tol) then
        ! a<>b, a<>c, b<>c
        BV(1)=((-a+b)*(a-c)*(b-c)*d2*abs(d)+b*(b-c)*c*(a+d)**2*abs(a+d)+a*(c*(-a+c)*(b+d)**2*abs(b+d)+ &
              (a-b)*b*((a-c)*(b-c)*c+(c+d)**2*abs(c+d))))/(12.d0*a*(a-b)*b*(a-c)*(b-c)*c)
      else if (abs(a-b)<tol .and. abs(a-c)>tol .and. abs(b-c)>tol) then
        ! a=b, a<>c, b<>c
        BV(1)=(-((a - c)**2*d2*abs(d)) + c*(a + d)*(a2 + c*d - 2*a*(c + d))*abs(a + d) + & 
              a2*((a - c)**2*c + (c + d)**2*abs(c + d)))/(12.*a2*(a - c)**2*c)
      else if (abs(a-b)>tol .and. abs(a-c)<tol .and. abs(b-c)>tol) then
        ! a<>b, a=c, b<>c
        BV(1)=(-((a - b)**2*d2*abs(d)) + b*(a + d)*(a2 + b*d - 2*a*(b + d))*abs(a + d) + & 
              a2*((a - b)**2*b + (b + d)**2*abs(b + d)))/(12.*a2*(a - b)**2*b)
      else if (abs(a-b)>tol .and. abs(a-c)>tol .and. abs(b-c)<tol) then
        ! a<>b, a<>c, b=c
        BV(1)=(-((a - b)**2*d2*abs(d)) + b2*(a + d)**2*abs(a + d) + &
              a*((a - b)**2*b2 + (b + d)*(b*(b - 2*d) + a*(-2*b + d))*abs(b + d)))/ &
              (12.*a*(a - b)**2*b2)
      end if
    end if
    eqcv=matmul(InvA30,BV)

  case (31)

    a2=a**2; b2=b**2; c2=c**2; d2=d**2
    a3=a**3; b3=b**3; c3=c**3; d3=d**3
    a4=a**4; b4=b**4; c4=c**4; d4=d**4
    a5=a**5; b5=b**5; c5=c**5; d5=d**5
    a6=a**6; b6=b**6; c6=c**6; d6=d**6
    if ( am<bm*tol .and. am<cm*tol .and. min(bm,cm)/max(bm,cm)>tol ) then
      ! a=0, b<>0, c<>0
      ! b c
      abs01=abs(b+c+d)
      abs02=abs(b+c-d)
      abs03=abs(b-c+d)
      abs04=abs(b-c-d)
      abs05=abs(-b+c+d)
      abs06=abs(-b+c-d)
      abs07=abs(-b-c+d)
      abs08=abs(-b-c-d)
      include 'hex_ph_23_bc.f90'
    else if ( bm<am*tol .and. bm<cm*tol .and. min(am,cm)/max(am,cm)>tol ) then
      ! a<>0, b=0, c<>0
      ! a c
      abs01=abs(a+c+d)
      abs02=abs(a+c-d)
      abs03=abs(a-c+d)
      abs04=abs(a-c-d)
      abs09=abs(-a+c+d)
      abs10=abs(-a+c-d)
      abs11=abs(-a-c+d)
      abs12=abs(-a-c-d)
      include 'hex_ph_23_ac.f90'
    else if ( cm<am*tol .and. cm<bm*tol .and. min(am,bm)/max(am,bm)>tol ) then
      ! a<>0, b<>0, c=0
      ! a b
      abs01=abs(a+b+d)
      abs02=abs(a+b-d)
      abs05=abs(a-b+d)
      abs06=abs(a-b-d)
      abs09=abs(-a+b+d)
      abs10=abs(-a+b-d)
      abs13=abs(-a-b+d)
      abs14=abs(-a-b-d)
      include 'hex_ph_23_ab.f90'
    else if ( am<cm*tol .and. bm<cm*tol ) then
      ! a=0, b=0, c<>0
      ! c
      abs01=abs(c+d)
      abs02=abs(c-d)
      abs03=abs(-c+d)
      abs04=abs(-c-d)
      include 'hex_ph_23_c.f90'
    else if ( bm<am*tol .and. cm<am*tol ) then
      ! a<>0, b=0, c=0
      ! a
      abs01=abs(a+d)
      abs02=abs(a-d)
      abs09=abs(-a+d)
      abs10=abs(-a-d)
      include 'hex_ph_23_a.f90'
    else if ( am<bm*tol .and. cm<bm*tol ) then
      ! a=0, b<>0, c=0
      ! b
      abs01=abs(b+d)
      abs02=abs(b-d)
      abs05=abs(-b+d)
      abs06=abs(-b-d)
      include 'hex_ph_23_b.f90'
    else
      ! a<>0, b<>0, c<>0
      ! a b c
      abs01=abs(a+b+c+d)
      abs02=abs(a+b+c-d)
      abs03=abs(a+b-c+d)
      abs04=abs(a+b-c-d)
      abs05=abs(a-b+c+d)
      abs06=abs(a-b+c-d)
      abs07=abs(a-b-c+d)
      abs08=abs(a-b-c-d)
      abs09=abs(-a+b+c+d)
      abs10=abs(-a+b+c-d)
      abs11=abs(-a+b-c+d)
      abs12=abs(-a+b-c-d)
      abs13=abs(-a-b+c+d)
      abs14=abs(-a-b+c-d)
      abs15=abs(-a-b-c+d)
      abs16=abs(-a-b-c-d)
      include 'hex_ph_23_abc.f90'
    end if
    eqcv=matmul(InvA31,BV)

  case default
    stop 'Invalid etype specified'

  end select

end subroutine Heqpol_coefficients
  