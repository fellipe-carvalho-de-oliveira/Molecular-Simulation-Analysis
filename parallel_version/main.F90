! static_factor_1 calculates S(q) directly from atoms positions
! static_factor_2 calculates S(q) by 1 + 4.0 * pi * rho * trapezios(y,distance_rdf,number_of_bins_gr)
! static_factor_3 calculates S(q) by fourier transforming g(r), it doesn't work yet

PROGRAM main  
USE nrtype
USE nrutil
USE sinal_processing    
USE functions
USE global_variables
USE rheology_routines
USE interpolation_fitting_methods
USE cluster_analysis

!REAL(8):: start,finish
!CALL cpu_time(start)

CALL read_data
CALL get_particle_positions

IF (Pc_flag == 1) THEN
  CALL Pc_function
END IF

IF (radial_distribution_flag == 1) THEN
  CALL radial_distri
END IF

  

IF (static_factor_flag == 1) THEN
    
  IF (static_factor_option == 1) THEN  
    CALL static_factor_1
  END IF
  
  IF (static_factor_option == 2) THEN  
    CALL static_factor_2
  END IF  
  
  IF (static_factor_option == 3) THEN  
    CALL static_factor_3
  END IF
  
END IF

IF (Fs_flag == 1) THEN
  CALL self_intermediate_scattering_function
END IF

IF (Fc_flag == 1) THEN
  CALL collective_intermediate_scattering_function
END IF

IF (MSD_flag == 1) THEN
  CALL MSD_total
END IF

IF (MSD_by_contact_number_flag == 1) THEN
  CALL MSD_by_contact_number
END IF

IF (Nc_by_contact_number_flag  == 1) THEN
  CALL Nc_by_contact_number 
END IF

IF (cluster_flag == 1) THEN
  CALL clustering
END IF


IF (flag_rheology == 1) THEN
  
  CALL getting_values_stress 
  
  IF (stress_acf_flag == 1) THEN
    
    !option 1 by definition, option 1 by fourier transform  
    IF (stress_acf_option == 1) THEN
      CALL stress_acf_definition
    END IF    
  
    IF (stress_acf_option == 2) THEN
      CALL stress_acf_FFT  
    END IF
  
  END IF
  
  IF (G_primes_flag == 1) THEN
    
    !option 1 GK, option 2 NEMD
    IF (G_primes_option == 1) THEN
        
      IF (stress_acf_flag == 0) THEN
        STOP "There is no Stress ACF for calculating G' , G''"    
      END IF   
        
      CALL rheological_properties_GK        
    END IF 
    
    IF (G_primes_option == 2) THEN
      CALL  stress_ave_NEMD
      CALL rheological_properties_NEMD       
    END IF
    
    
  END IF  
  CALL dealloc_rheology_vectors  
END IF    


!CALL static_factor_3
!CALL contacs_number
!CALL MSD_by_contact_number
!CALL self_intermediate_scattering_function_FFT
CALL dealloc

!WRITE(*,*)
!CALL cpu_time(finish)
!PRINT '("Time = ",f12.4," minutes.")',(finish-start)/60.0
END PROGRAM main    
