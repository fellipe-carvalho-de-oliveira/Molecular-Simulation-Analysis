!     
! File:   rheology_routines.F90
! Author: fellipe
!
! Created on September 18, 2015, 5:33 PM
!

MODULE rheology_routines
USE global_variables
USE FFTs
USE swarm_routines
USE functions


CONTAINS

SUBROUTINE rheological_properties_NEMD
  IMPLICIT NONE
  INTEGER :: s , i , pos, n_swarm 
  REAL(8) :: aux1, aux2, Gd_correct , phase_correct , tau_correct , G_prime, G_prime_2 , omega
  REAL(8), DIMENSION(correl_leng) :: strain_rate
  
  OPEN(unit=10,file="fitting_stress_NEMD.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening fitting_stress_NEMD.out"
  CLOSE(10)
  
  aux1 = 0.0; aux2 = 0.0;
  n_swarm = 3                           ! Número de vezes pra rodar o swarm e fazer média
  DO i = 1, n_swarm
   CALL swarm
   aux1 = aux1 + stress_ampli
   aux2 = aux2 + fase
  END DO
  !CALL newton  
  
  stress_ampli = aux1 / REAL(n_swarm,8)
  fase = aux2 / REAL(n_swarm,8)
  
  OPEN(unit=10,file="fitting_stress_NEMD.out",status="old",iostat = status_open, position="append") 
  IF (status_open > 0) STOP "error opening fitting_stress_NEMD.out"
  WRITE(10,*) "Amplitude média da tensão = ", stress_ampli
  WRITE(10,*) "Fase média = ", fase
  CLOSE(10)
  
  Gd_correct = abs(stress_ampli) / strain_ampli  
  phase_correct = fase
  tau_correct = stress_ampli
  
  G_prime = Gd_correct * cos(abs(phase_correct))
  G_prime_2 = Gd_correct * sin(abs(phase_correct))
  
  omega = 2.0*pi/REAL(tp,8)/delta_t_stress
  
  DO s = 1,  correl_leng
    stress_fit(s) = tau_correct *  sin(omega * time_stress(s)  + phase_correct) 
  END DO
  
  OPEN(unit=1,file="stress_fitted.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "stress_fitted.out"
  
  DO i = 1, correl_leng
    strain_rate(i) = (strain_ampli*omega)*cos(omega*time_stress(i))
  END DO
  
  WRITE(1,"(3x,a4,22x,a6,19x,a6,20x,a13,14x,a11,15x,a9)") "time", "Strain", "Stress", "Stress fitted" , "Strain rate" , &
  "Viscosity"
  DO i = 1, correl_leng
    WRITE(1,*) time_stress(i) , strain(i) , stress_acf_ave(i,1) , stress_fit(i) &
    , strain_rate(i) , abs(stress_fit(i) / strain_rate(i) ) 
  END DO
  
  OPEN(unit=2,file="rheological_properties_NEMD.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "rheological_properties_NEMD.out"
  
    WRITE(2,*) "G' : " , G_prime 
    WRITE(2,*) "G'' : " , G_prime_2
    WRITE(2,*) "Dynamic viscosity : " , G_prime_2 / omega
  
  CLOSE(1)
  CLOSE(2)
  
END SUBROUTINE  rheological_properties_NEMD

SUBROUTINE stress_ave_NEMD
INTEGER :: i
REAL(8), DIMENSION(correl_leng) :: a

a = 0.0

DO i = 1,numb_perio
  DO j = 1,correl_leng
    a(j) = a(j) + stress(j + correl_leng * (i - 1),1)
  END DO    
END DO    

DO j = 1,correl_leng
  stress_acf_ave(j,1) = a(j) / numb_perio
END DO  

END SUBROUTINE stress_ave_NEMD    

SUBROUTINE rheological_properties_GK
  IMPLICIT NONE
  INTEGER :: j , jj
  REAL(8) :: distance , a , i , ang_velo , delta_phase , strain_0 , stress_0 , G_prime, G_prime_2, volume , orale 
  REAL(8), DIMENSION(correl_leng) :: aux1, aux2,aux3, aux4,aux5, aux6
  REAL(8), DIMENSION(6) ::viscosity_no_shear
  COMPLEX(8), DIMENSION(6) :: complex_modulus 
  COMPLEX(8), DIMENSION(correl_leng,6) :: stress_acf_complex 
  COMPLEX(8), DIMENSION(correl_leng) :: aux_complex_1, aux_complex_2,aux_complex_3,aux_complex_4,aux_complex_5,aux_complex_6

  volume = L * L * L
  viscosity_no_shear = volume / kb / Temp
  DO j = 1,correl_leng
    aux1(j) = stress_acf_ave(j,1)
    aux2(j) = stress_acf_ave(j,2)  
    aux3(j) = stress_acf_ave(j,3)  
    aux4(j) = stress_acf_ave(j,4)  
    aux5(j) = stress_acf_ave(j,5)  
    aux6(j) = stress_acf_ave(j,6)   
  END DO    
  viscosity_no_shear(1) = viscosity_no_shear(1) * trapezios(aux1,time_stress,correl_leng)
  viscosity_no_shear(2) = viscosity_no_shear(2) * trapezios(aux2,time_stress,correl_leng)
  viscosity_no_shear(3) = viscosity_no_shear(3) * trapezios(aux3,time_stress,correl_leng)
  viscosity_no_shear(4) = viscosity_no_shear(4) * trapezios(aux4,time_stress,correl_leng)
  viscosity_no_shear(5) = viscosity_no_shear(5) * trapezios(aux5,time_stress,correl_leng)
  viscosity_no_shear(6) = viscosity_no_shear(6) * trapezios(aux6,time_stress,correl_leng)
  
  
  ang_velo = 2 * pi / (tp * delta_t_stress) ! 2*pi/ (time of 1 period)  
  
  ! G' G'' by complex modulus
  
  DO j = 1,correl_leng
    stress_acf_complex(j,1) = exp(-1.0 * (0.0,1.0) * ang_velo * time_stress(j)) * stress_acf_ave(j,1)
    aux_complex_1(j) = stress_acf_complex(j,1)
    stress_acf_complex(j,2) = exp(-1.0 * (0.0,1.0) * ang_velo * time_stress(j)) * stress_acf_ave(j,2) 
    aux_complex_2(j) = stress_acf_complex(j,2)
    stress_acf_complex(j,3) = exp(-1.0 * (0.0,1.0) * ang_velo * time_stress(j)) * stress_acf_ave(j,3) 
    aux_complex_3(j) = stress_acf_complex(j,3)
    stress_acf_complex(j,4) = exp(-1.0 * (0.0,1.0) * ang_velo * time_stress(j)) * stress_acf_ave(j,4) 
    aux_complex_4(j) = stress_acf_complex(j,4)
    stress_acf_complex(j,5) = exp(-1.0 * (0.0,1.0) * ang_velo * time_stress(j)) * stress_acf_ave(j,5) 
    aux_complex_5(j) = stress_acf_complex(j,5)
    stress_acf_complex(j,6) = exp(-1.0 * (0.0,1.0) * ang_velo * time_stress(j)) * stress_acf_ave(j,6) 
    aux_complex_6(j) = stress_acf_complex(j,6)

  END DO
  
  complex_modulus(1) = (volume / kb / Temp) * (0.0, 1.0) * ang_velo * &
  trapezios_complex(aux_complex_1,time_stress,correl_leng)
  complex_modulus(2) = (volume / kb / Temp) * (0.0, 1.0) * ang_velo * &
  trapezios_complex(aux_complex_2,time_stress,correl_leng)
  complex_modulus(3) = (volume / kb / Temp) * (0.0, 1.0) * ang_velo * &
  trapezios_complex(aux_complex_3,time_stress,correl_leng)
  complex_modulus(4) = (volume / kb / Temp) * (0.0, 1.0) * ang_velo * &
  trapezios_complex(aux_complex_4,time_stress,correl_leng)
  complex_modulus(5) = (volume / kb / Temp) * (0.0, 1.0) * ang_velo * &
  trapezios_complex(aux_complex_5,time_stress,correl_leng)
  complex_modulus(6) = (volume / kb / Temp) * (0.0, 1.0) * ang_velo * &
  trapezios_complex(aux_complex_6,time_stress,correl_leng)
  
  OPEN(unit=1,file="rheological_properties_gk.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "rheological_properties_gk.out"
  
  WRITE(1,*) "G* em xx : " , complex_modulus(1)
  WRITE(1,*) "G* em yy : " , complex_modulus(2)
  WRITE(1,*) "G* em zz : " , complex_modulus(3)
  WRITE(1,*) "G* em xy : " , complex_modulus(4)
  WRITE(1,*) "G* em xz : " , complex_modulus(5)
  WRITE(1,*) "G* em yz : " , complex_modulus(6)
  WRITE(1,*) "viscosidade no shear em xx : " ,viscosity_no_shear(1)
  WRITE(1,*) "viscosidade no shear em yy : " ,viscosity_no_shear(2)
  WRITE(1,*) "viscosidade no shear em zz : " ,viscosity_no_shear(3)
  WRITE(1,*) "viscosidade no shear em xy : " ,viscosity_no_shear(4)
  WRITE(1,*) "viscosidade no shear em xz : " ,viscosity_no_shear(5)
  WRITE(1,*) "viscosidade no shear em yz : " ,viscosity_no_shear(6)
  WRITE(1,*) "Shear Viscosity ( (xy + xz + yz) / 3 ) : ", &
             ( viscosity_no_shear(4) + viscosity_no_shear(5) + viscosity_no_shear(6) ) / 3.0
  
  
  CLOSE(1)

END SUBROUTINE  rheological_properties_GK

SUBROUTINE stress_acf_definition
  IMPLICIT NONE
  INTEGER :: i , j , k , kk , ii , initial_frame , final_frame
  REAL(8) :: distance  
  REAL(8), DIMENSION(6) :: a 
  
  OPEN(unit=10,file="stress_acf_definition.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "stress_acf_definition.out"
  
  WRITE(10,'(3x,a5,13x,a24)') "frame","stress_acf by definition"
  
  DO ii = 1, numb_perio 
      
    initial_frame  = (ii - 1) * correl_leng + 1
    final_frame = ii * correl_leng
    
    DO k = initial_frame,final_frame,stress_acf_Tgap   !delta t
      a = 0.0  
    
      !calculating autocorrelation function
      DO kk = initial_frame,final_frame - k + 1 + correl_leng*(ii - 1)      !number of origins
        a(1) = a(1) + stress(k + kk - 1 - correl_leng*(ii - 1),1) * stress(kk,1) 
        a(2) = a(2) + stress(k + kk - 1 - correl_leng*(ii - 1),2) * stress(kk,2) 
        a(3) = a(3) + stress(k + kk - 1 - correl_leng*(ii - 1),3) * stress(kk,3) 
        a(4) = a(4) + stress(k + kk - 1 - correl_leng*(ii - 1),4) * stress(kk,4) 
        a(5) = a(5) + stress(k + kk - 1 - correl_leng*(ii - 1),5) * stress(kk,5) 
        a(6) = a(6) + stress(k + kk - 1 - correl_leng*(ii - 1),6) * stress(kk,6) 
      END DO
      
      !End calculating autocorrelation function
      
      stress_acf_ave(k - correl_leng*(ii - 1),1) = stress_acf_ave(k - correl_leng*(ii - 1),1) + a(1) / (final_frame - k + 1) 
      stress_acf_ave(k - correl_leng*(ii - 1),2) = stress_acf_ave(k - correl_leng*(ii - 1),2) +  a(2) / (final_frame - k + 1) 
      stress_acf_ave(k - correl_leng*(ii - 1),3) = stress_acf_ave(k - correl_leng*(ii - 1),3) +  a(3) / (final_frame - k + 1) 
      stress_acf_ave(k - correl_leng*(ii - 1),4) = stress_acf_ave(k - correl_leng*(ii - 1),4) +  a(4) / (final_frame - k + 1) 
      stress_acf_ave(k - correl_leng*(ii - 1),5) = stress_acf_ave(k - correl_leng*(ii - 1),5) + a(5) / (final_frame - k + 1) 
      stress_acf_ave(k - correl_leng*(ii - 1),6) = stress_acf_ave(k - correl_leng*(ii - 1),6) + a(6) / (final_frame - k + 1)  
      
    END DO
    
  END DO
  
  DO k = 1,correl_leng
    DO j = 1,6  
      stress_acf_ave(k,j) = stress_acf_ave(k,j) / numb_perio    
    END DO
    if (stress_acf_ave(k,1) /= 0.0) then
    WRITE(10,*) k - 1 , stress_acf_ave(k,1), stress_acf_ave(k,2), stress_acf_ave(k,3), stress_acf_ave(k,4), stress_acf_ave(k,5), &
                stress_acf_ave(k,6)
    end if            
  END DO
  
  
  
  CLOSE(10)
  

END SUBROUTINE  stress_acf_definition

SUBROUTINE stress_acf_FFT
  IMPLICIT NONE
  INTEGER :: i , ii, j , k , initial_frame, final_frame
  REAL(8), DIMENSION(correl_leng) :: aux1, aux2, aux3, aux4, aux5, aux6  
  REAL(8), DIMENSION(correl_leng) :: aux7, aux8, aux9, aux10, aux11, aux12
  
  aux1 = 0.0; aux2 = 0.0; aux3 = 0.0; aux4 = 0.0; aux5 = 0.0; aux6 = 0.0; aux7 = 0.0; aux8 = 0.0; aux9 = 0.0; aux10 = 0.0;
  aux11 = 0.0; aux12 = 0.0;
  
  OPEN(unit=10,file="stress_acf_FFT.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "stress_acf_FFT.out"
  
  DO ii = 1, numb_perio
      
    initial_frame  = (ii - 1) * correl_leng + 1
    final_frame = ii * correl_leng  
      
    DO j = initial_frame,final_frame
      aux1(j - correl_leng*(ii - 1)) = stress(j,1)
      aux2(j - correl_leng*(ii - 1)) = stress(j,2)  
      aux3(j - correl_leng*(ii - 1)) = stress(j,3)  
      aux4(j - correl_leng*(ii - 1)) = stress(j,4)  
      aux5(j - correl_leng*(ii - 1)) = stress(j,5)  
      aux6(j - correl_leng*(ii - 1)) = stress(j,6)  
    END DO
  
    CALL correl(aux1,aux1,aux7) 
    CALL correl(aux2,aux2,aux8) 
    CALL correl(aux3,aux3,aux9) 
    CALL correl(aux4,aux4,aux10) 
    CALL correl(aux5,aux5,aux11) 
    CALL correl(aux6,aux6,aux12) 
  
    DO j = 1,correl_leng
      stress_acf_ave(j,1) = stress_acf_ave(j,1) + aux7(j)
      stress_acf_ave(j,2) = stress_acf_ave(j,2) + aux8(j) 
      stress_acf_ave(j,3) = stress_acf_ave(j,3) + aux9(j)
      stress_acf_ave(j,4) = stress_acf_ave(j,4) + aux10(j)
      stress_acf_ave(j,5) = stress_acf_ave(j,5) + aux11(j)
      stress_acf_ave(j,6) = stress_acf_ave(j,6) + aux12(j)
    END DO    
  END DO
  
  WRITE(10,'(3x,a5,13x,a17)') "frame","stress_acf by FFT"
  
  DO k = 1,correl_leng
    DO j = 1,6  
      stress_acf_ave(k,j) = stress_acf_ave(k,j) / numb_perio    
    END DO
    WRITE(10,*) k - 1 , stress_acf_ave(k,1), stress_acf_ave(k,2), stress_acf_ave(k,3), stress_acf_ave(k,4), stress_acf_ave(k,5), &
                stress_acf_ave(k,6)
  END DO

  CLOSE(10)
END SUBROUTINE  stress_acf_FFT

SUBROUTINE alloc_rheology_vectors
  ALLOCATE(stress(nframes_stress,6), strain(correl_leng), time_stress(correl_leng))    
  ALLOCATE(correlation_stress(correl_leng) ) 
  ALLOCATE(stress_fit(correl_leng))
  ALLOCATE(stress_acf_ave(correl_leng,6))
  stress_acf_ave = 0.0
END SUBROUTINE alloc_rheology_vectors 

SUBROUTINE dealloc_rheology_vectors
  DEALLOCATE(stress, strain,time_stress,correlation_stress,stress_fit,stress_acf_ave)     

END SUBROUTINE dealloc_rheology_vectors 

    
SUBROUTINE getting_values_stress
  IMPLICIT NONE  
  
  INTEGER :: control = 0 , i , j  
  REAL(8) :: aux1,aux2
  REAL(8), DIMENSION(:) , ALLOCATABLE :: aux3,aux4
  
  nframes_stress = 0
  
  OPEN(unit=1,file=stress_file ,status="old",iostat = status_open ) 
  
  IF (status_open > 0) STOP "error opening rheological input file"

  READ (1,*) ; READ (1,*)
  
  DO WHILE (control == 0)
     nframes_stress =  nframes_stress + 1 !getting number of lines
    READ (1 , * , iostat= control) 
  END DO
  
  REWIND 1 
  
  READ (1,*) ; READ (1,*) 
  
  READ (1,*) aux1
  READ (1,*) aux2
  dump_interval_stress = aux2 - aux1
  
  REWIND 1 
  
  READ (1,*) ; READ (1,*)   
  
  nframes_stress = nframes_stress - 1
  
  ALLOCATE(aux3(nframes_stress),aux4(nframes_stress))
  
  CALL alloc_rheology_vectors
  
  IF (G_primes_option == 1) THEN
    DO i = 1,nframes_stress
      READ (1,*) aux3(i),stress(i,1),stress(i,2),stress(i,3),stress(i,4),stress(i,5),stress(i,6)
    END DO

    DO i = 1,correl_leng
      time_stress(i) = aux3(i) * delta_t_stress
    END DO
    
  END IF 
    
  IF (G_primes_option == 2) THEN
    DO i = 1,nframes_stress
      READ (1,*) aux3(i),aux4(i),stress(i,1)
    END DO  
    
    DO i = 1,correl_leng
      time_stress(i) = aux3(i) * delta_t_stress
      strain(i) = aux4(i) * L                     ! Estou multiplicando por L, porque no lammps pedi para printar deformação/L
    END DO
    
  END IF
  


  CLOSE(1)
  
  WRITE(*,*) ; WRITE(*,*) "Rheological Properties"  
  numb_perio = nframes_stress / correl_leng
  WRITE(*,*) "Number of periods is", numb_perio
  
  DEALLOCATE(aux3)
END SUBROUTINE getting_values_stress    


END MODULE rheology_routines
