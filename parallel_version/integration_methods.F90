!     
! File:   integration_methods.F90
! Author: fellipe
!
! Created on October 2, 2015, 10:50 AM
!

MODULE integration_methods
USE nrtype; USE nrutil; USE global_variables
IMPLICIT NONE

CONTAINS

   FUNCTION trapezios(y,x,tamanho)
   IMPLICIT NONE
   REAL(8), dimension(tamanho), intent(inout)  :: y,x
   INTEGER :: i,tamanho
   REAL(8) :: trapezios, aux, volume
   trapezios = 0.0 
   volume = L * L * L
   
!   OPEN(unit=1,file="viscosity.out",status="replace",iostat = status_open) 
!   IF (status_open > 0) STOP "error opening viscosity.out"
!   WRITE(1,*) "timestep", " Viscosity"
   
   DO i = 1,tamanho-1  
     trapezios = trapezios + (y(i) + y(i+1))*(x(i+1) - x(i))/2.0
!     IF (y(i) < 0.0) THEN 
!       RETURN 
!     END IF    
!     WRITE(1,*) i, trapezios * volume / kb / Temp
     
   END DO   
   
!   CLOSE(1)
   
END 

FUNCTION bacle(y,x,tamanho)
   IMPLICIT NONE
   REAL(8), DIMENSION(tamanho)  :: y,x
   INTEGER :: i,tamanho, interval
   REAL(8) :: bacle, aux , h
   bacle = 0.0
   
   IF (mod(tamanho,5) /= 0) THEN
       STOP " Bacle's rule needs a multiple of 5 points to integrate"
   END IF    
   
   interval = tamanho / 5
   
   h = x(2) - x(1)
   
   DO i = 1,interval
     bacle = bacle+( 7*y(1+5*(i-1)) + 32*y(2+5*(i-1)) + 12*y(3+5*(i-1)) + 32*y(4+5*(i-1)) + 7*y(5+5*(i-1)) )*h*2./45.
     IF (y(i) < 0.0) THEN 
       RETURN 
     END IF 
   END DO   
   
END 

FUNCTION trapezios_complex(y,x,tamanho)
   IMPLICIT NONE
   COMPLEX(8), DIMENSION(tamanho)  :: y
   REAL(8), DIMENSION(tamanho)  :: x
   INTEGER :: i,tamanho
   COMPLEX(8) :: trapezios_complex
   trapezios_complex = 0.0
   
   DO i = 1,tamanho - 1  
     trapezios_complex = trapezios_complex + (y(i) + y(i+1))*(x(i+1) - x(i))/2.
     IF (abs(y(i)) < 0.0) THEN 
       RETURN 
     END IF 
   END DO   
   
END

!rend la derivée du vecteur y dans le vecteur deriv_y de taille size(y) - 2
SUBROUTINE derivative(y,x,deriv_y)
IMPLICIT NONE    
REAL(8), DIMENSION(:), INTENT(INOUT) :: y,x
REAL(8), DIMENSION(size(y) - 2), INTENT(OUT) :: deriv_y
INTEGER :: i
REAL(8) :: h

h = x(2) - x(1)
deriv_y(1) = (y(2) - y(1)) / h
deriv_y(2) = (y(3) - y(1)) / (2 * h)  

DO i=3,size(y) - 2
deriv_y(i) = ( 8.0*y(i+1) - 8.0*y(i-1) + 1.0*y(i-2) - 1.0*y(i+2)  ) / ( 12 * h )
END DO    
    
END SUBROUTINE derivative




END MODULE
