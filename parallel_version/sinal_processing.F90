!     
! File:   sinal_processing.F90
! Author: fellipe
!
! Created on October 20, 2015, 4:16 PM
!

MODULE sinal_processing
USE nrtype; USE nrutil ; USE global_variables
IMPLICIT NONE

CONTAINS

SUBROUTINE smooth_rdf_average
INTEGER :: i,k

OPEN(unit=10,file="rdf_smoothed_average.out",status="replace",iostat = status_open ) 
IF (status_open > 0) STOP "error opening file rdf_smoothed_average.out"

DO i = 1,nframes
  
  DO k = 1, number_of_bins_gr - num_poin_appro_rdf + 1  
    rdf_smooth(k,i) = ( rdf(k,i) + rdf(k+1,i) + rdf(k+2,i) + rdf(k+3,i) + rdf(k+4,i) ) / 5.0      
  END DO 
  
  IF (i == 1 .or. i == nframes .or. (mod(i-1,gr_Tgap) == 0) ) THEN
    WRITE(10,*) " timestep ", i - 1
    WRITE(10,*) " distance                   g(r) smoothed"
      DO k = 1, number_of_bins_gr - num_poin_appro_rdf + 1 
        WRITE(10,*) distance_rdf(k) , rdf_smooth(k,i)
      END DO
    WRITE(10,*)
    WRITE(10,*)
  END IF
  
END DO 

CLOSE(10)    
END SUBROUTINE smooth_rdf_average 

!SUBROUTINE y_smooth_rdf_average(y)
!REAL(8), DIMENSION(:), INTENT(INOUT) :: y    
!INTEGER :: i,k
!
!OPEN(unit=10,file="y_smoothed_average.out",status="replace",iostat = status_open ) 
!IF (status_open > 0) STOP "error opening file y_smoothed_average.out"
!
!DO i = 1,nframes
!  
!  DO k = 1, number_of_bins_gr - num_poin_appro_rdf + 1  
!    rdf_smooth(k,i) = ( rdf(k,i) + rdf(k+1,i) + rdf(k+2,i) + rdf(k+3,i) + rdf(k+4,i) ) / 5.0      
!  END DO 
!  
!  IF (i == 1 .or. i == nframes .or. (mod(i-1,gr_Tgap) == 0) ) THEN
!    WRITE(10,*) " timestep ", i - 1
!    WRITE(10,*) " distance                   g(r) smoothed"
!      DO k = 1, number_of_bins_gr - num_poin_appro_rdf + 1 
!        WRITE(10,*) distance_rdf(k) , rdf_smooth(k,i)
!      END DO
!    WRITE(10,*)
!    WRITE(10,*)
!  END IF
!  
!END DO 
!
!CLOSE(10)    
!END SUBROUTINE y_smooth_rdf_average

SUBROUTINE smooth_rdf_allen_tilde
INTEGER :: i , k , n

n = number_of_bins_gr

OPEN(unit=10,file="rdf_smoothed_allen_tilde.out",status="replace",iostat = status_open ) 
IF (status_open > 0) STOP "error opening file rdf_smoothed_allen_tilde.out"

DO i = 1,nframes
  
  rdf_smooth(1,i) = (1.0 / 70.0) * (69.0 * rdf(1,i) + 4.0 * rdf(2,i) - 6.0 * rdf(3,i) + 4.0 * rdf(4,i) - rdf(5,i) )
  rdf_smooth(2,i) = (1.0 / 35.0) * (2.0 * rdf(1,i) + 27.0 * rdf(2,i) + 12.0 * rdf(3,i) - 8.0 * rdf(4,i) + 2.0 * rdf(5,i) )
      
  DO k = 3, n - 2  
    rdf_smooth(k,i) = (1.0 / 35.0) * &
    (-3.0 * rdf(k-2,i) + 12.0 * rdf(k-1,i) + 17.0 * rdf(k,i) + 12.0 * rdf(k+1,i) - 3.0 * rdf(k+2,i) )    
  END DO 
  
  rdf_smooth(n - 1,i) = (1.0 / 35.0) * & 
  (2.0 * rdf(n,i) + 27.0 * rdf(n - 1,i) + 12.0 * rdf(n - 2,i) - 8.0 * rdf(n - 3,i) + 2.0 * rdf(n - 4,i) )
  
  rdf_smooth(n,i) = (1.0 / 35.0) * &
  (69.0 * rdf(n,i) + 4.0 * rdf(n-1,i) - 6.0 * rdf(n-2,i) + 4.0 * rdf(n-3,i) - rdf(n-4,i) )
  
  IF (i == 1 .or. i == nframes .or. (mod(i-1,gr_Tgap) == 0) ) THEN
    WRITE(10,*) " timestep ", i - 1
    WRITE(10,*) " distance                   g(r) smoothed"
      DO k = 1, number_of_bins_gr - num_poin_appro_rdf + 1 
        WRITE(10,*) distance_rdf(k) , rdf_smooth(k,i)
      END DO
    WRITE(10,*)
    WRITE(10,*)
  END IF
  
END DO 

CLOSE(10)
END SUBROUTINE smooth_rdf_allen_tilde

END MODULE sinal_processing
