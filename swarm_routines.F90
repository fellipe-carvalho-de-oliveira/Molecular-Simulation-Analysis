!     
! File:   swarm_routines.F90
! 
!
! Created on September 21, 2015, 8:58 AM
!

MODULE swarm_routines
USE global_variables
USE nrtype
USE nrutil

    
CONTAINS   

SUBROUTINE swarm
  IMPLICIT NONE  
  REAL(8) :: c1, c2, Sglobal , num , aux
  REAL(8), DIMENSION(3) :: omega
  REAL(8), DIMENSION(:), ALLOCATABLE  :: xglobal , Sbest
  REAL(8), DIMENSION(:,:), ALLOCATABLE  :: Const, xbest, V  
  INTEGER :: i, j, t 
  INTEGER :: k = 1500 !Number of particles
  INTEGER :: nv = 2 !Number of parameters
  INTEGER :: it = 6 !Number of iterations
  REAL(8), DIMENSION(2,2) ::  minmax !first line stores the minimum value of each variable
                                                 !second line stores the maximum value of each variable   
  minmax(1,1) = -1e2
  minmax(2,1) = 1e2
  minmax(1,2) = -pi/2.0
  minmax(2,2) = pi/2.0
  
  ALLOCATE(x(k,nv) , xbest(k,nv), xglobal(nv), const(2,nv), V(k,nv) , Sbest(k) )
  x = 0.0 
  xbest = 0.0 
  xglobal = 0.0
  const = 0.0
  V = 0.0
  Sbest = 0.0
  
  OPEN(unit=2,file="fitting_stress_NEMD.out",status="old",iostat = status_open ,position="append" ) 
  IF (status_open > 0) STOP "error opening fitting_stress_NEMD.out"
  
  WRITE(2,*) "Swarm"
  WRITE(2,"(10X,A3,18X,A5,17X,11A,12X,A7)") "tau","phase","soma erro^2","iteração"
  
  c1 = 1.5
  c2 = 1.5
  ! omega[2] must be greater than omega[3]
  omega(2) = 1.0
  omega(3) = 0.4
  
  
  !first iteration
  Sglobal = 1e30 !Sglobal initialization
  t=0
  
  ! i refers to the particle, j to the parameter
  DO i=1,k

   DO j =1, nv 
     CALL init_random_seed()         
     CALL RANDOM_NUMBER(num)     
     Const(1,j) = num 
     x(i,j) = minmax(1,j)+Const(1,j)*(minmax(2,j)-minmax(1,j))
     xbest(i,j) = x(i,j) 
   
     V(i,j)=0.0
   END DO

   CALL Objfun(i)
   Sbest(i)=S

   IF (Sbest(i) < Sglobal) THEN
     Sglobal = Sbest(i)
     
     DO j = 1,nv
       xglobal(j) = xbest(i,j)
     END DO
     
   END IF
  END DO

  !Printing global best variables and objetive function values
  
  !escrever todos os xglobal,sglobal e t na mesma linha
  
  DO i = 1,nv 
    WRITE(2,"(2X,EN20.3)",ADVANCE = "NO") xglobal(i)     
  END DO    
  WRITE(2,"(2X,EN20.3)",ADVANCE = "NO") Sglobal
  WRITE(2,"(2X,I10)") t
  
  !other iterations
    
  DO t = 1,it
    
    omega(1) = omega(2) + (omega(3)-omega(2))*((t-1)/(it-1))

    DO i = 1,k
  
      DO j = 1,nv
        CALL init_random_seed()  
        CALL RANDOM_NUMBER(num)
        Const(1,j) = num
        CALL init_random_seed()  
        CALL RANDOM_NUMBER(num)
        Const(2,j) = num

        V(i,j) = omega(1) * V(i,j)+c1*Const(1,j)*(xbest(i,j)-x(i,j))+c2*Const(2,j)*(xglobal(j)-x(i,j))
     
        x(i,j) = x(i,j) + V(i,j)
        IF (x(i,j) > minmax(2,j)) THEN
          x(i,j) = minmax(2,j)
          V(i,j) = 0
        END IF  
        IF (x(i,j) < minmax(1,j)) THEN
          x(i,j) = minmax(1,j)
          V(i,j) = 0
        END IF
      END DO

      CALL Objfun(i)
      IF (S<Sbest(i)) THEN
    
        Sbest(i)=S
        DO j=1,nv
          xbest(i,j) = x(i,j)
        END DO
        IF (Sbest(i)<Sglobal) THEN
      
          Sglobal=Sbest(i)
          DO j = 1,nv
            xglobal(j) = xbest(i,j)
          END DO
        END IF
      END IF
    END DO
    
    DO j = 1,nv 
      WRITE(2,"(2X,EN20.3)",ADVANCE = "NO") xglobal(j)     
    END DO    
    WRITE(2,"(2X,EN20.3)",ADVANCE = "NO") Sglobal
    WRITE(2,"(2X,I10)") t

  END DO

  stress_ampli = xglobal(1) ! nv is 1, that is there's just one parameter
  fase = xglobal(2)
  
  WRITE(2,*); WRITE(2,*);  
  CLOSE(2)  
  DEALLOCATE(x, xbest, xglobal, const, V, Sbest)
END SUBROUTINE swarm

SUBROUTINE Objfun(k) 
  INTEGER ::  i
  INTEGER, INTENT(IN) :: k
  REAL(8) :: aux
  
  !k refers to the particle, 1 refers to the only parameter
  
  S = 0.0  
  
  DO i = 1,correl_leng
    S = S + ( stress_acf_ave(i,1) - x(k,1) * sin(2.0 * pi * time_stress(i) / ( tp * delta_t_stress) + x(k,2) ) )  ** 2.0
  END DO    
  
END SUBROUTINE Objfun


 subroutine init_random_seed()
            use iso_fortran_env, only: int64
            implicit none
            integer, allocatable :: seed(:)
            integer :: i, n, un, istat, dt(8), pid
            integer(int64) :: t
          
            call random_seed(size = n)
            allocate(seed(n))
            ! First try if the OS provides a random number generator
            open(newunit=un, file="/dev/urandom", access="stream", &
                 form="unformatted", action="read", status="old", iostat=istat)
            if (istat == 0) then
               read(un) seed
               close(un)
            else
               ! Fallback to XOR:ing the current time and pid. The PID is
               ! useful in case one launches multiple instances of the same
               ! program in parallel.
               call system_clock(t)
               if (t == 0) then
                  call date_and_time(values=dt)
                  t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                       + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                       + dt(3) * 24_int64 * 60 * 60 * 1000 &
                       + dt(5) * 60 * 60 * 1000 &
                       + dt(6) * 60 * 1000 + dt(7) * 1000 &
                       + dt(8)
               end if
               pid = getpid()
               t = ieor(t, int(pid, kind(t)))
               do i = 1, n
                  seed(i) = lcg(t)
               end do
            end if
            call random_seed(put=seed)
          contains
            ! This simple PRNG might not be good enough for real work, but is
            ! sufficient for seeding a better PRNG.
            function lcg(s)
              integer :: lcg
              integer(int64) :: s
              if (s == 0) then
                 s = 104729
              else
                 s = mod(s, 4294967296_int64)
              end if
              s = mod(s * 279470273_int64, 4294967291_int64)
              lcg = int(mod(s, int(huge(0), int64)), kind(0))
            end function lcg
end subroutine init_random_seed

SUBROUTINE newton
  INTEGER :: i , count = 0
  REAL(8) :: omega , stress_ampli_0 , fase_0 , residuo 
  REAL(8), DIMENSION(2) :: grad_Fobj
  REAL(8), DIMENSION(2,2) :: hess_Fobj
  
  OPEN(unit=2,file="fitting_stress_NEMD.out",status="old",iostat = status_open , position="append") 
  IF (status_open > 0) STOP "error opening fitting_stress_NEMD.out"
  
  WRITE(2,*) "Método de Newton"
  WRITE(2,"(10X,A3,18X,A5,17X,11A,12X,A7)") "tau","phase","soma erro^2","iteração"
  
  omega = 2.0 * pi / ( tp * delta_t_stress)
  
  DO count = 1,1000
  
    residuo = 0.0  
      
    stress_ampli_0  = stress_ampli
    fase_0 = fase
  
    grad_Fobj = 0.0
    hess_Fobj = 0.0
  
    DO i = 1,correl_leng
      grad_Fobj(1) = grad_Fobj(1) + (stress_acf_ave(i,1) - stress_ampli * sin(omega * time_stress(i) + fase )) * &
                   (-sin(omega * time_stress(i) + fase))
                   
      grad_Fobj(2) = grad_Fobj(2) + (stress_acf_ave(i,1) - stress_ampli * sin(omega * time_stress(i) + fase )) * &
                   (-stress_ampli) * (cos(omega * time_stress(i) + fase))   
                   
      hess_Fobj(1,1) = hess_Fobj(1,1) + sin(omega * time_stress(i) + fase) ** 2.0   
    
      hess_Fobj(1,2) = hess_Fobj(1,2) + (stress_acf_ave(i,1) - stress_ampli * sin(omega * time_stress(i) + fase )) * (-1.0) * &
                     cos(omega * time_stress(i) + fase) + (-stress_ampli) * cos(omega * time_stress(i) + fase) * (-1.0) * &
                     sin(omega * time_stress(i) + fase)
                     
      hess_Fobj(2,1) = hess_Fobj(2,1) + cos(omega * time_stress(i) + fase) * (-sin(omega * time_stress(i) + fase) *(-stress_ampli) &
                     + (stress_acf_ave(i,1) - stress_ampli * sin(omega * time_stress(i) + fase )) * (-1.0) )   
                     
      hess_Fobj(2,2) = hess_Fobj(2,2) + (stress_acf_ave(i,1) - stress_ampli * sin(omega * time_stress(i) + fase )) * & 
                     (-sin(omega * time_stress(i) + fase)) - stress_ampli * cos(omega * time_stress(i) + fase) ** 2.0
      
    END DO    
  
    grad_Fobj(1) = grad_Fobj(1) * 2.0
    grad_Fobj(2) = grad_Fobj(2) * 2.0
    hess_Fobj(1,1) = hess_Fobj(1,1) * 2.0
    hess_Fobj(1,2) = hess_Fobj(1,2) * 2.0
    hess_Fobj(2,1) = hess_Fobj(2,1) * 2.0
    hess_Fobj(2,2) = hess_Fobj(2,2) * (-2.0) * stress_ampli
  
    stress_ampli = stress_ampli - (hess_Fobj(1,1)*grad_Fobj(1) + hess_Fobj(1,2)*grad_Fobj(2))
    fase = fase - (hess_Fobj(2,1)*grad_Fobj(1) + hess_Fobj(2,2)*grad_Fobj(2))
    
    residuo = residuo + abs(stress_ampli - stress_ampli_0)
    residuo = residuo + abs(fase - fase_0) 
    
    S = 0.0
    DO i = 1,correl_leng
      S = S + ( stress_acf_ave(i,1) - stress_ampli * sin(omega * time_stress(i) + fase ) )  ** 2.0
    END DO
    
    WRITE(2,"(2X,EN20.3)",ADVANCE = "NO") stress_ampli    
    WRITE(2,"(2X,EN20.3)",ADVANCE = "NO") fase   
    WRITE(2,"(2X,EN20.3)",ADVANCE = "NO") S
    WRITE(2,"(2X,I10)") count
    

    IF (residuo < 1e-2) THEN
      CLOSE(2)  
      WRITE(*,*) "Newton convergiu com  ", count, " iterações ", ", resíduo = ", residuo
      RETURN
    END IF    
  
  END DO
  
  
  CLOSE(2)
  IF (count > 999) STOP "Método de Newton não convergiu"
  
END SUBROUTINE newton


END MODULE swarm_routines
