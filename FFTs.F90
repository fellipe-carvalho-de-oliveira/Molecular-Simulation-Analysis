!     
! File:   FFTs.F90
! Author: fellipe
!
! Created on September 15, 2015, 3:41 PM
!

MODULE FFTs
USE nrtype
USE nrutil

IMPLICIT NONE

CONTAINS


!In place Cooley-Tukey FFT of a complex vector, replaces the original vector
recursive subroutine fft(x)
    complex(kind=dp), dimension(:), intent(inout)  :: x
    complex(kind=dp)                               :: t
    integer                                        :: N
    integer                                        :: i
    complex(kind=dp), dimension(:), allocatable    :: even, odd
 
    N=size(x)
 
    if(N .le. 1) return
 
    allocate(odd((N+1)/2))
    allocate(even(N/2))
 
    ! divide
    odd =x(1:N:2)
    even=x(2:N:2)
 
    ! conquer
    call fft(odd)
    call fft(even)
 
    ! combine
    do i=1,N/2
       t=exp(cmplx(0.0_dp,-2.0_dp*pi*real(i-1,dp)/real(N,dp),kind=dp))*even(i)
       x(i)     = odd(i) + t
       x(i+N/2) = odd(i) - t
    end do
 
    deallocate(odd)
    deallocate(even)
 
end subroutine fft

!Normalized IFFT of a complex vector,replaces the original vector
subroutine ifft(x)
    complex(kind=dp), dimension(:), intent(inout)  :: x 
    integer                                        :: N
    N=size(x)
    
    
    x = conjg(x)
    call fft(x)
    x = conjg(x)
    x = x / N
    
end subroutine ifft

!Performs the time correlation function of two real vectors and stores the values in correlation_real
!x_real and y_real must be size that is power of 2
SUBROUTINE correl(x_real,y_real,correlation_real)
  complex(8), dimension(:), ALLOCATABLE :: x_complex,y_complex
  real(8), dimension(:), intent(inout)  :: x_real,y_real,correlation_real
  integer                               :: Nt, i , status_open
  COMPLEX(8) :: aux1 , aux2 
  Nt=size(x_real)
   
  ALLOCATE(x_complex(2*Nt),y_complex(2*Nt))
   
  DO i = 1,Nt
    x_complex(i) = x_real(i)
    y_complex(i) = y_real(i)
  END DO 
    
  DO i = Nt + 1,2*Nt
    x_complex(i) = 0.0
    y_complex(i) = 0.0
  END DO
   
  CALL fft(x_complex)
  CALL fft(y_complex)
    
  !x_complex(1)=cmplx(real(x_complex(1))*real(y_complex(1))/no2,aimag(x_complex(1))*aimag(y_complex(1))/no2, kind=spc)   

  x_complex = x_complex*conjg(y_complex)
    
  CALL ifft(x_complex)
    
  DO i = 1,Nt 
    x_complex(i) = x_complex(i) / ( (real(Nt - i + 1)))
    correlation_real(i) = real(x_complex(i))      
  END DO
  
  DEALLOCATE(x_complex,y_complex)
END SUBROUTINE correl

!Performs the time correlation function of two real vectors and stores the values in correlation_real
!x_real and y_real must be size that is power of 2
SUBROUTINE correl_complex(x_input,y_input,correlation_real)
  complex(8), dimension(:), ALLOCATABLE :: x_complex,y_complex
  real(8), dimension(:), intent(inout)  :: correlation_real
  complex(8), dimension(:), intent(inout)  :: x_input,y_input
  integer                               :: Nt, i , status_open
  COMPLEX(8) :: aux1 , aux2 
  Nt=size(x_input)
   
  ALLOCATE(x_complex(2*Nt),y_complex(2*Nt))
   
  DO i = 1,Nt
    x_complex(i) = x_input(i)
    y_complex(i) = y_input(i)
  END DO 
    
  DO i = Nt + 1,2*Nt
    x_complex(i) = 0.0
    y_complex(i) = 0.0
  END DO
   
  CALL fft(x_complex)
  CALL fft(y_complex)
    
  !x_complex(1)=cmplx(real(x_complex(1))*real(y_complex(1))/no2,aimag(x_complex(1))*aimag(y_complex(1))/no2, kind=spc)   

  x_complex = x_complex*conjg(y_complex) 
    
  CALL ifft(x_complex)
    
  DO i = 1,Nt 
    x_complex(i) = x_complex(i) / ( (real(Nt - i + 1)))
    correlation_real(i) = abs(x_complex(i))      
  END DO
   
  DEALLOCATE(x_complex,y_complex)

END SUBROUTINE correl_complex


!pega um vetor dataaa REAL e aplica transformada de fourier aos seus valores retornando o mesmo vetor com os valores
!transformados no lugar dos originais MOSTRA APENAS A PARTE REAL, se isign = 1.
! se o vator COMPLEXO zdata( N / 2) é fornecido o vetor dataa não é mudado
!e ele conterá os valores transformados, N ainda tem que ser potencia de 2 em todo caso
!se isign = -1 calcula transformada inversa
SUBROUTINE realft_dp(dataa,isign,zdata)

REAL(DP), DIMENSION(:), INTENT(INOUT) :: dataa
INTEGER(I4B), INTENT(IN) :: isign
COMPLEX(DPC), DIMENSION(:), OPTIONAL, TARGET :: zdata
INTEGER(I4B) :: n,ndum,nh,nq
COMPLEX(DPC), DIMENSION(size(dataa)/4) :: w
COMPLEX(DPC), DIMENSION(size(dataa)/4-1) :: h1,h2
COMPLEX(DPC), DIMENSION(:), POINTER :: cdata
COMPLEX(DPC) :: z
REAL(DP) :: c1=0.5_dp,c2
n=size(dataa)
call assert(iand(n,n-1)==0, 'n must be a power of 2 in realft_dp')
nh=n/2
nq=n/4
if (present(zdata)) then
ndum=assert_eq(n/2,size(zdata),'realft_dp')
cdata=>zdata
if (isign == 1) cdata=cmplx(dataa(1:n-1:2),dataa(2:n:2),kind=spc)
else
allocate(cdata(n/2))
cdata=cmplx(dataa(1:n-1:2),dataa(2:n:2),kind=spc)
end if
if (isign == 1) then
c2=-0.5_dp
call four1_dp(cdata,+1)
else
c2=0.5_dp
end if
w=zroots_unity(sign(n,isign),n/4)
w=cmplx(-aimag(w),real(w),kind=dpc)
h1=c1*(cdata(2:nq)+conjg(cdata(nh:nq+2:-1)))
h2=c2*(cdata(2:nq)-conjg(cdata(nh:nq+2:-1)))
cdata(2:nq)=h1+w(2:nq)*h2
cdata(nh:nq+2:-1)=conjg(h1-w(2:nq)*h2)
z=cdata(1)
if (isign == 1) then
cdata(1)=cmplx(real(z)+aimag(z),real(z)-aimag(z),kind=dpc)
else
cdata(1)=cmplx(c1*(real(z)+aimag(z)),c1*(real(z)-aimag(z)),kind=dpc)
call four1_dp(cdata,-1)
end if
if (present(zdata)) then
if (isign /= 1) then
dataa(1:n-1:2)=real(cdata)
dataa(2:n:2)=aimag(cdata)
end if
else
dataa(1:n-1:2)=real(cdata)
dataa(2:n:2)=aimag(cdata)
deallocate(cdata)
end if
END SUBROUTINE realft_dp

!pega um vetor dataaa COMPLEX e aplica transformada de fourier aos seus valores retornando o mesmo vetor com os valores
!transformados no lugar dos originais
SUBROUTINE four1_dp(dataa,isign)

COMPLEX(DPC), DIMENSION(:), INTENT(INOUT) :: dataa
INTEGER(I4B), INTENT(IN) :: isign
COMPLEX(DPC), DIMENSION(:,:), ALLOCATABLE :: dat,temp
COMPLEX(DPC), DIMENSION(:), ALLOCATABLE :: w,wp
REAL(DP), DIMENSION(:), ALLOCATABLE :: theta
INTEGER(I4B) :: n,m1,m2,j
n=size(dataa)
call assert(iand(n,n-1)==0, 'n must be a power of 2 in four1_dp')
m1=2**ceiling(0.5_sp*log(real(n,sp))/0.693147_sp)
m2=n/m1
allocate(dat(m1,m2),theta(m1),w(m1),wp(m1),temp(m2,m1))
dat=reshape(dataa,shape(dat))
call fourrow_dp(dat,isign)
theta=arth(0,isign,m1)*TWOPI_D/n
wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
w=cmplx(1.0_dp,0.0_dp,kind=dpc)
do j=2,m2
w=w*wp+w
dat(:,j)=dat(:,j)*w
end do
temp=transpose(dat)
call fourrow_dp(temp,isign)
dataa=reshape(temp,shape(dataa))
deallocate(dat,w,wp,theta,temp)
END SUBROUTINE four1_dp

!pega uma matriz dataaa COMPLEX e aplica transformada de fourier Às suas linhas, retornando a mesma matriz com as linhas
!agora com os valores transformados se isign = 1, dataa deve ter o número de dados que é potência de 2
SUBROUTINE fourrow_dp(dataaa,isign)
COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: dataaa
INTEGER(I4B), INTENT(IN) :: isign
INTEGER(I4B) :: n,i,istep,j,m,mmax,n2
REAL(DP) :: theta
COMPLEX(DPC), DIMENSION(size(dataaa,1)) :: temp
COMPLEX(DPC) :: w,wp
COMPLEX(DPC) :: ws
n=size(dataaa,2)
call assert(iand(n,n-1)==0, 'n must be a power of 2 in fourrow_dp')
n2=n/2
j=n2
do i=1,n-2
if (j > i) call swap(dataaa(:,j+1),dataaa(:,i+1))
m=n2
do
if (m < 2 .or. j < m) exit
j=j-m
m=m/2
end do
j=j+m
end do
mmax=1
do
if (n <= mmax) exit
istep=2*mmax
theta=PI_D/(isign*mmax)
wp=cmplx(-2.0_dp*sin(0.5_dp*theta)**2,sin(theta),kind=dpc)
w=cmplx(1.0_dp,0.0_dp,kind=dpc)
do m=1,mmax
ws=w
do i=m,n,istep
j=i+mmax
temp=ws*dataaa(:,j)
dataaa(:,j)=dataaa(:,i)-temp
dataaa(:,i)=dataaa(:,i)+temp
end do
w=w*wp+w
end do
mmax=istep
end do
END SUBROUTINE fourrow_dp
    


END MODULE FFTs
