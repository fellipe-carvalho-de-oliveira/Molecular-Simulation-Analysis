!     
! File:   interpolation_fitting_methods.F90
! Author: fellipe
!
! Created on October 27, 2015, 6:31 PM
!

MODULE interpolation_fitting_methods
USE nrtype; USE nrutil 
IMPLICIT NONE

CONTAINS


!Given same-size arrays x and y containing a tabulated function yi = f(xi), this routine
!returns a same-size array of coefficients cj .
FUNCTION polcoe(x,y)
REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
REAL(DP), DIMENSION(size(x)) :: polcoe
INTEGER(I4B) :: i,k,n
REAL(DP), DIMENSION(size(x)) :: s
REAL(DP), DIMENSION(size(x),size(x)) :: a
n=assert_eq(size(x),size(y),"polcoe")

s=0.0 
s(n)=-x(1)
do i=2,n
s(n+1-i:n-1)=s(n+1-i:n-1)-x(i)*s(n+2-i:n)
s(n)=s(n)-x(i)
end do

a=outerdiff(x,x) 
polcoe=product(a,dim=2,mask=a/= 0.0) 
a(:,1)=-s(1)/x(:)
do k=2,n
a(:,k)=-(s(k)-a(:,k-1))/x(:)
end do
s=y/polcoe
polcoe=matmul(s,a) 
END FUNCTION polcoe


!Given arrays x and y of length N containing a tabulated function, i.e., yi = f (xi ), with x1 <
!x2 < . . . < xN , and given values yp1 and ypn for the first derivative of the interpolating
!function at points 1 and N , respectively, this routine returns an array y2 of length N
!that contains the second derivatives of the interpolating function at the tabulated points
!xi . If yp1 and/or ypn are equal to 1 × 10e30 or larger, the routine is signaled to set the
!corresponding boundary condition for a natural spline, with zero second derivative on that
!boundary.

SUBROUTINE spline(x,y,yp1,ypn,y2)

REAL(DP), DIMENSION(:), INTENT(IN) :: x,y
REAL(DP), INTENT(IN) :: yp1,ypn
REAL(DP), DIMENSION(:), INTENT(OUT) :: y2
INTEGER(I4B) :: n
REAL(DP), DIMENSION(size(x)) :: a,b,c,r
n=assert_eq(size(x),size(y),size(y2),'spline')
c(1:n-1)=x(2:n)-x(1:n-1)
r(1:n-1)=6.0_dp*((y(2:n)-y(1:n-1))/c(1:n-1))
r(2:n-1)=r(2:n-1)-r(1:n-2)
a(2:n-1)=c(1:n-2)
b(2:n-1)=2.0_dp*(c(2:n-1)+a(2:n-1))
b(1)=1.0
b(n)=1.0
if (yp1 > 0.99e30_dp) then
r(1)=0.0
c(1)=0.0
else
r(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
c(1)=0.5
end if
if (ypn > 0.99e30_dp) then
r(n)=0.0
a(n)=0.0
else
r(n)=(-3.0_dp/(x(n)-x(n-1)))*((y(n)-y(n-1))/(x(n)-x(n-1))-ypn)
a(n)=0.5
end if
call tridag(a(2:n),b(1:n),c(1:n-1),r(1:n),y2(1:n))
END SUBROUTINE spline

!Given the arrays xa and ya, which tabulate a function (with the xai ’s in increasing or
!decreasing order), and given the array y2a, which is the output from spline above, and
!given a value of x, this routine returns a cubic-spline interpolated value. The arrays xa, ya
!and y2a are all of the same size.

FUNCTION splint(xa,ya,y2a,x)

REAL(DP), DIMENSION(:), INTENT(IN) :: xa,ya,y2a
REAL(DP), INTENT(IN) :: x
REAL(DP) :: splint
INTEGER(I4B) :: khi,klo,n
REAL(DP) :: a,b,h
n=assert_eq(size(xa),size(ya),size(y2a),'splint')
klo=max(min(locate(xa,x),n-1),1)
khi=klo+1
h=xa(khi)-xa(klo)
if (h == 0.0) call nrerror('bad xa input in splint')
a=(xa(khi)-x)/h
b=(x-xa(klo))/h
splint=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp
END FUNCTION splint

!Given an array xx(1:N ), and given a value x, returns a value j such that x is between
!xx(j ) and xx(j + 1). xx must be monotonic, either increasing or decreasing. j = 0 or
!1j = N is returned to indicate that x is out of range.

FUNCTION locate(xx,x)
USE nrtype
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(IN) :: xx
REAL(DP), INTENT(IN) :: x
INTEGER(I4B) :: locate
INTEGER(I4B) :: n,jl,jm,ju
LOGICAL :: ascnd
n=size(xx)
ascnd = (xx(n) >= xx(1))
jl=0
ju=n+1
do
if (ju-jl <= 1) exit
jm=(ju+jl)/2
if (ascnd .eqv. (x >= xx(jm))) then
jl=jm
else
ju=jm
end if
end do
if (x == xx(1)) then
locate=1
else if (x == xx(n)) then
locate=n-1
else
locate=jl
end if
END FUNCTION locate


!Solves for a vector u of size N the tridiagonal linear set given by equation (2.4.1) using a
!serial algorithm. Input vectors b (diagonal elements) and r (right-hand sides) have size N ,
!while a and c (off-diagonal elements) are size N − 1.

SUBROUTINE tridag(a,b,c,r,u)
USE nrtype; USE nrutil, ONLY : assert_eq,nrerror
IMPLICIT NONE
REAL(DP), DIMENSION(:), INTENT(IN) :: a,b,c,r
REAL(DP), DIMENSION(:), INTENT(OUT) :: u
REAL(DP), DIMENSION(size(b)) :: gam
INTEGER(I4B) :: n,j
REAL(DP) :: bet
n=assert_eq((/size(a)+1,size(b),size(c)+1,size(r),size(u)/),'tridag')
bet=b(1)
if (bet == 0.0) call nrerror('tridag_ser: Error at code stage 1')
u(1)=r(1)/bet
do j=2,n
gam(j)=c(j-1)/bet
bet=b(j)-a(j-1)*gam(j)
if (bet == 0.0) &
call nrerror('tridag_ser: Error at code stage 2')
u(j)=(r(j)-a(j-1)*u(j-1))/bet
end do
do j=n-1,1,-1
u(j)=u(j)-gam(j+1)*u(j+1)
end do
END SUBROUTINE tridag


END MODULE interpolation_fitting_methods
