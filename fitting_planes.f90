!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   fitting_planes.f90
! Author: fellipe
!
! Created on 21 de Novembro de 2017, 13:52
!

MODULE fitting_planes
    USE nrtype
    USE nrutil
    IMPLICIT NONE
    
    CONTAINS
    
    SUBROUTINE  fit_planes(x,y,z,a,b,c,SUM_SQ)
      REAL(8), DIMENSION(:), INTENT(INOUT)  :: x,y,z
      REAL(8), INTENT(INOUT)  :: a,b,c,SUM_SQ
      REAL(8) :: sumx, sumy, sumz, Sxx, Syy, Szz, Sxy, Sxz, Syz, c0, c1, c2, c3, r_zinho, s, t, p, q, R , a1,a2,a3,b1,b2,b3
      REAL(8) :: SUM_SQ1, SUM_SQ2, SUM_SQ3, rho, phi
      REAL(8) :: pi = 3.14159265359, lambda, delta
      INTEGER :: i,j,k,n,pos
      REAL(8), DIMENSION(:), ALLOCATABLE :: x_red, y_red, z_red
      REAL(8), DIMENSION(3) :: compara_vector3
      REAL(8), DIMENSION(2) :: compara_vector2
      
!      REAL(8), DIMENSION(4) :: aaa
!      REAL(8), DIMENSION(3) :: aa
!      REAL(8), DIMENSION(3) :: rti,rtr
      
      n = size(x)
      ALLOCATE(x_red(n),y_red(n),z_red(n))     
      
      !do i=1,n
      !    write(*,*) x(i),y(i),z(i)
      !end do    
      
      x_red = 0.0
      y_red = 0.0
      z_red = 0.0 
      sumx = 0.0
      sumy = 0.0
      sumz = 0.0
      Sxx = 0.0
      Syy = 0.0
      Szz = 0.0
      Sxy = 0.0
      Sxz = 0.0
      Syz = 0.0
      SUM_SQ = 0.0
      SUM_SQ1 = 0.0
      SUM_SQ2 = 0.0
      SUM_SQ3 = 0.0
      
      
      DO i = 1,n
        sumx = sumx + x(i)
        sumy = sumy + y(i)
        sumz = sumz + z(i)
      END DO    
      
      sumx = sumx/real(n,8)
      sumy = sumy/real(n,8)
      sumz = sumz/real(n,8)
      
      !CALL swarm(x_red, y_red, z_red,a,b,SUM_SQ)
      DO i = 1,n
        x_red(i) = x(i) - sumx
        IF (abs(x_red(i)) < 1.0e-12) x_red = 0.0
        y_red(i) = y(i) - sumy
        IF (abs(y_red(i)) < 1.0e-12) y_red = 0.0
        z_red(i) = z(i) - sumz
        IF (abs(z_red(i)) < 1.0e-12) z_red = 0.0
      END DO 
      
      DO i = 1,n
        Sxx = Sxx + x_red(i)*x_red(i)
        Syy = Syy + y_red(i)*y_red(i)
        Szz = Szz + z_red(i)*z_red(i)
        Sxy = Sxy + x_red(i)*y_red(i)
        Sxz = Sxz + x_red(i)*z_red(i)
        Syz = Syz + y_red(i)*z_red(i)        
      END DO 

      c0 = Syz*(Sxy*Sxy-Sxz*Sxz) + Sxy*Sxz*(Szz-Syy)
      
      c1 = Sxy ** 3.0 + Sxy * (Sxz ** 2.0 - 2.0 * Syz ** 2.0 - Szz ** 2.0)
      c1 = c1 + Sxy * (Sxx * Szz + Syy * Szz - Sxx * Syy) + Sxz * Syz * (Syy + Szz - 2.0 * Sxx)
      
      c2 = Syz ** 3.0 + Syz * (Sxz ** 2.0 - 2.0 * Sxy ** 2.0 - Sxx ** 2.0)
      c2 = c2 + Syz * (Sxx * Szz + Sxx * Syy - Syy * Szz) + Sxy * Sxz * (Sxx + Syy - 2.0 * Szz)
      
      c3 = Sxy*(Syz*Syz-Sxz*Sxz) + Sxz*Syz*(Sxx-Syy)
     
      IF (c3 /= 0.0) THEN
         
       r_zinho = c2/c3
       s = c1/c3
       t = c0/c3
       p = s - r_zinho * r_zinho / 3.0
       q = 2 * (r_zinho ** 3.0) / 27.0 - r_zinho * s / 3.0 + t
       R = q * q / 4.0 + p * p * p / 27.0
      
       IF (R > 0.0) THEN
        
         a = - r_zinho / 3.0 + (sqrt(R)- q / 2.0) ** (1.0/3.0) + (-sqrt(R)- q / 2.0) ** (1.0/3.0) 
        
         b = Sxy*Syz*a*a + (Syz*Syz - Sxy*Sxy)*a - Sxy*Syz
        
         b = b / ((Syz*(Sxx-Syy) - Sxy*Sxz)*a + Sxy*(Syy-Szz) + Sxz*Syz)
        
         DO i = 1,n
           SUM_SQ = SUM_SQ + (a*x_red(i) + b*y_red(i) + z_red(i))**2.0  
         END DO
        
         SUM_SQ = SUM_SQ / (a*a+b*b+1.0)          
       END IF 
      
       IF (R < 0.0) THEN
        
         rho = sqrt(-p*p*p/27.0)
         phi = acos(-q/2.0/rho)
        
         a1 = - r_zinho / 3.0 + 2.0 * ( (rho) ** (1.0/3.0) ) * cos(phi/3.0)
         a2 = - r_zinho / 3.0 + 2.0 * ( (rho) ** (1.0/3.0) ) * cos((phi+2.0*pi)/3.0)
         a3 = - r_zinho / 3.0 + 2.0 * ( (rho) ** (1.0/3.0) ) * cos((phi+4.0*pi)/3.0)
        
         b1 = Sxy*Syz*a1*a1 + (Syz*Syz - Sxy*Sxy)*a1 - Sxy*Syz        
         b1 = b1 / ((Syz*(Sxx-Syy) - Sxy*Sxz)*a1 + Sxy*(Syy-Szz) + Sxz*Syz)
        
         b2 = Sxy*Syz*a2*a2 + (Syz*Syz - Sxy*Sxy)*a2 - Sxy*Syz        
         b2 = b2 / ((Syz*(Sxx-Syy) - Sxy*Sxz)*a2 + Sxy*(Syy-Szz) + Sxz*Syz)  
        
         b3 = Sxy*Syz*a3*a3 + (Syz*Syz - Sxy*Sxy)*a3 - Sxy*Syz        
         b3 = b3 / ((Syz*(Sxx-Syy) - Sxy*Sxz)*a3 + Sxy*(Syy-Szz) + Sxz*Syz)        
        
         DO i = 1,n
           SUM_SQ1 = SUM_SQ1 + (a1*x_red(i) + b1*y_red(i) + z_red(i))**2.0 
           SUM_SQ2 = SUM_SQ2 + (a2*x_red(i) + b2*y_red(i) + z_red(i))**2.0  
           SUM_SQ3 = SUM_SQ3 + (a3*x_red(i) + b3*y_red(i) + z_red(i))**2.0  
         END DO
        
         SUM_SQ1 = SUM_SQ1 / (a1*a1+b1*b1+1.0)  
         SUM_SQ2 = SUM_SQ2 / (a2*a2+b2*b2+1.0)     
         SUM_SQ3 = SUM_SQ3 / (a3*a3+b3*b3+1.0)
        
         compara_vector3(1) = SUM_SQ1
         compara_vector3(2) = SUM_SQ2
         compara_vector3(3) = SUM_SQ3 
        
         pos = minloc(compara_vector3,DIM=1)
        
         IF (pos == 1) THEN
           a = a1
           b = b1
           SUM_SQ = SUM_SQ1
         END IF   
        
         IF (pos == 2) THEN
           a = a2
           b = b2
           SUM_SQ = SUM_SQ2
         END IF   
        
         IF (pos == 3) THEN
           a = a3
           b = b3
           SUM_SQ = SUM_SQ3
         END IF   
        
       END IF  
      
      !aaa(1) = c3;aaa(2) = c2;aaa(3) = c1;aaa(4) = c0
      !CALL zrhqr(aaa,rtr,rti)    
          
      lambda = (1.0/real(n,8)) * (a*sumx+b*sumy+sumz)
      c = -1.0/lambda
      b = -b/lambda
      a = -a/lambda
      RETURN
      
      END IF ! if do c3
      
      IF ((c3 == 0.0) .AND. (c2 /= 0.0)) THEN
        delta = c1*c1-4.0*c2*c0

        a1 = (-c1+sqrt(delta) )/ (2.0*c2)  
        a2 = (-c1-sqrt(delta) )/ (2.0*c2) 
        b1 = Sxy*Syz*a1*a1 + (Syz*Syz - Sxy*Sxy)*a1 - Sxy*Syz        
        b1 = b1 / ((Syz*(Sxx-Syy) - Sxy*Sxz)*a1 + Sxy*(Syy-Szz) + Sxz*Syz)
        b2 = Sxy*Syz*a2*a2 + (Syz*Syz - Sxy*Sxy)*a2 - Sxy*Syz        
        b2 = b2 / ((Syz*(Sxx-Syy) - Sxy*Sxz)*a2 + Sxy*(Syy-Szz) + Sxz*Syz)          
        
        DO i = 1,n
          SUM_SQ1 = SUM_SQ1 + (a1*x_red(i) + b1*y_red(i) + z_red(i))**2.0 
          SUM_SQ2 = SUM_SQ2 + (a2*x_red(i) + b2*y_red(i) + z_red(i))**2.0  
        END DO
        
        SUM_SQ1 = SUM_SQ1 / (a1*a1+b1*b1+1.0)  
        SUM_SQ2 = SUM_SQ2 / (a2*a2+b2*b2+1.0)     
        
        compara_vector2(1) = SUM_SQ1
        compara_vector2(2) = SUM_SQ2
        
        pos = minloc(compara_vector2,DIM=1)
        
        IF (pos == 1) THEN
          a = a1
          b = b1
          SUM_SQ = SUM_SQ1
        END IF   
        
        IF (pos == 2) THEN
          a = a2
          b = b2
          SUM_SQ = SUM_SQ2
        END IF    
      
      lambda = (1.0/real(n,8)) * (a*sumx+b*sumy+sumz)
      c = -1.0/lambda
      b = -b/lambda
      a = -a/lambda
      RETURN     
        
      END IF    
      
      IF ((c3 == 0.0) .AND. (c2 == 0.0) .AND. (c1 /= 0.0)) THEN      
        a = -c0/c1
        b = Sxy*Syz*a*a + (Syz*Syz - Sxy*Sxy)*a - Sxy*Syz       
        b = b / ((Syz*(Sxx-Syy) - Sxy*Sxz)*a + Sxy*(Syy-Szz) + Sxz*Syz)
        
        DO i = 1,n
          SUM_SQ = SUM_SQ + (a*x_red(i) + b*y_red(i) + z_red(i))**2.0  
        END DO
        
        SUM_SQ = SUM_SQ / (a*a+b*b+1.0)    
        lambda = (1.0/real(n,8)) * (a*sumx+b*sumy+sumz)
        c = -1.0/lambda
        b = -b/lambda
        a = -a/lambda
          
        RETURN
      END IF
      
      IF (c1 == 0.0) STOP "All coefficients of the polynom fitting are 0.0"
      RETURN
      
      
      
    END SUBROUTINE fit_planes
    
    

END MODULE fitting_planes
