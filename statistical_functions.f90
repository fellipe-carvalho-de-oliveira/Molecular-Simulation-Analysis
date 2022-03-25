! File:   statistical_functions.f90
! Author: fellipe
!
! Created on 20 de Novembro de 2017, 18:13
!

MODULE statistical_functions
   IMPLICIT NONE
  
   CONTAINS

! --------------------------------------------------------------------
! SUBROUTINE  statistical(data,mean,std_dev):
!    This subroutine computes the mean, variance, standard deviation,
! sum and sum-of-squares of the vetor "data" :
!    (1) Sum       - sum of input values
!    (2) SumSQR    - sun-of-squares
!    (3) n         - number of input data items
!    (4) Mean      - computed mean value
!    (5) Variance  - computed variance
!    (6) StdDev    - computed standard deviation
! --------------------------------------------------------------------

   SUBROUTINE  statistical(dados,mean,std_dev)
      IMPLICIT  NONE
      REAL(8), DIMENSION(:), INTENT(INOUT)  :: dados
      INTEGER :: n , i
      REAL(8) :: Sum, SumSQR , variance
      REAL(8), INTENT(INOUT) :: mean, std_dev
      n = size(dados)
      
      Sum = 0.0
      SumSQR = 0.0
      DO i = 1,n
        Sum    = Sum + dados(i)
        SumSQR = SumSQR + dados(i)*dados(i)
      END DO
      
      mean = Sum / real(n,8)
      variance = (SumSQR - Sum*Sum/real(n,8))/(real(n,8)-1.0)
      std_dev   = SQRT(Variance)
      
   END SUBROUTINE  statistical

END MODULE statistical_functions
