!     
! File:   functions.F90
! Author: fellipe
!
! r: is the cuttof for  calculating the correlation (half of box size)  
!
! Created on May 20, 2015, 9:16 AM
!

MODULE functions
USE nrtype
USE nrutil
USE sinal_processing     
USE FFTs    
USE global_variables
USE swarm_routines
USE integration_methods
USe interpolation_fitting_methods
 
CONTAINS

SUBROUTINE contacs_number 
  IMPLICIT NONE
  INTEGER :: i , ii , j , jj , k , kk  , posi1 , posi2
  REAL(8) :: dist_neigh
  REAL(8), DIMENSION(:), ALLOCATABLE  :: rdf_bond , aux_bond 
  
  ALLOCATE(rdf_bond(number_of_bins_gr - 1))
  
  contacs_number_flag = 1
  
  IF (radial_distribution_flag == 0) THEN
    CALL radial_distri
  END IF
  
  DO i = 1,nframes     
      
    DO j = 1,number_of_bins_gr-1   
      rdf_bond(j) = rdf(j,i) 
    END DO
    posi1 = 0
    posi2 = 0
    
    CALL ordenar_cres(rdf_bond,posi1) 
    
    ALLOCATE(aux_bond(number_of_bins_gr-1 - posi1))
      
    DO j = posi1 + 1,number_of_bins_gr-1    
      aux_bond(j - posi1) = rdf(j,i)
    END DO    
    
    CALL ordenar_decres(aux_bond,posi2)
    
    posi2 = 2 !!!!
    
    bond_length(i) = distance_rdf(posi1 + posi2)
    
    DO k = 1,natoms
      DO kk = 1,natoms 
        dist_neigh = sqrt( (rxyz(kk,1,i) - rxyz(k,1,i)) ** 2 + (rxyz(kk,2,i) - rxyz(k,2,i)) ** 2 + &
        (rxyz(kk,3,i) - rxyz(k,3,i)) ** 2 )   
        IF ( (dist_neigh <= bond_length(i)) .and. (dist_neigh /= 0.0) .and. (bond_length(i) < cutoff) ) THEN
            Nc(k,i) = Nc(k,i) + 1
        END IF    
      END DO
    END DO
    
    
    DEALLOCATE(aux_bond)
  END DO 
  
  DEALLOCATE(rdf_bond)
END SUBROUTINE contacs_number


SUBROUTINE Pc_function
  IMPLICIT NONE
  INTEGER :: i , j , k , kk , status_open
  
  IF (contacs_number_flag == 0) THEN
    CALL  contacs_number
  END IF 
  
  Pc_flag = 1
  
  OPEN(unit=4,file="Probability_contacts.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file Probability_contacts.out"
    
  Pc = 0.0
  !Calculating Probabilities
  DO i = 1,nframes       
    DO j = 1,natoms
      SELECT CASE (Nc(j,i))
        CASE (0)
          Pc(1,i) =  Pc(1,i) + 1.0  
        CASE (1)
          Pc(2,i) =  Pc(2,i) + 1.0  
        CASE (2)
          Pc(3,i) =  Pc(3,i) + 1.0  
        CASE (3)
          Pc(4,i) =  Pc(4,i) + 1.0
        CASE (4)
          Pc(5,i) =  Pc(5,i) + 1.0  
        CASE (5)
          Pc(6,i) =  Pc(6,i) + 1.0
        CASE (6)
          Pc(7,i) =  Pc(7,i) + 1.0          
        CASE (7)
          Pc(8,i) =  Pc(8,i) + 1.0              
        CASE (8)
          Pc(9,i) =  Pc(9,i) + 1.0  
        CASE (9)
          Pc(10,i) =  Pc(10,i) + 1.0                      
        CASE (10)
          Pc(11,i) =  Pc(11,i) + 1.0
        CASE (11)
          Pc(12,i) =  Pc(12,i) + 1.0                             
        CASE (12)                                 
          Pc(13,i) =  Pc(13,i) + 1.0             
      END SELECT
    END DO
    
    DO j = 1,13
      Pc(j,i) = Pc(j,i) / natoms
    END DO   
    
    WRITE(4,*) "frame    ",i
    WRITE(4,*) "Contacts number    " , "Pc"
    DO j = 1,13
      WRITE(4,*) j - 1 , "      " , Pc(j,i) 
    END DO
    WRITE(4,*) 
    WRITE(4,*)
    
  END DO 
  !End Calculating Probabilities
  
  CLOSE(4)     
  
END SUBROUTINE Pc_function

!SUBROUTINE relative_proba_by_contact_number
!  OPEN(unit=10,file="Relative_proba(group_size).out",status="replace",iostat = status_open )
!  IF (status_open > 0) STOP "error opening Relative_proba(group_size).out"
!  
!  WRITE(10,*) "frame    ",i
!  WRITE(10,*) "Contacts number    " , "P(Nc+1) / P(Nc)"
!  DO j = 1,12
!    WRITE(10,*) j - 1 , "      " , Pc(j+1,i) / Pc(j,i)
!  END DO
!  WRITE(10,*) 
!  WRITE(10,*)
!    
!  CLOSE(10)   
!  
!END SUBROUTINE relative_proba_by_contact_number



SUBROUTINE create_contact_group !creates the groups and calculates the number of particle in each of them
  IMPLICIT NONE
  INTEGER :: i , j , k , kk , posi 
  
  create_contact_group_flag = 1
  
  !Calculating the frequency of occurence for each atom
  DO j = 1,natoms
   count_part = 0.0
    DO i = 1,nframes  
      SELECT CASE (Nc(j,i))
      CASE (0)
        count_part(1) = count_part(1) + 1   
      CASE (1)
        count_part(2) = count_part(2) + 1  
      CASE (2)
        count_part(3) = count_part(3) + 1  
      CASE (3)
        count_part(4) = count_part(4) + 1  
      CASE (4)
        count_part(5) = count_part(5) + 1  
      CASE (5)
        count_part(6) = count_part(6) + 1  
      CASE (6)
        count_part(7) = count_part(7) + 1  
      CASE (7)
        count_part(8) = count_part(8) + 1  
      CASE (8)
        count_part(9) = count_part(9) + 1  
      CASE (9)
        count_part(10) = count_part(10) + 1  
      CASE (10)
        count_part(11) = count_part(11) + 1  
      CASE (11)
        count_part(12) = count_part(12) + 1  
      CASE (12)
        count_part(13) = count_part(13) + 1 
      END SELECT  
    END DO
    
    CALL ordenar_cres(count_part,posi)
    index_part(j) = posi - 1
    
  END DO
  !End calculating the frequency of occurence for each atom
   
  !creating count_part vector, number of particles in a bin
  count_part = 0.0
  DO j = 1,natoms
    SELECT CASE (index_part(j))
    CASE (0)
      count_part(1) = count_part(1) + 1   
    CASE (1)
      count_part(2) = count_part(2) + 1  
    CASE (2)
      count_part(3) = count_part(3) + 1  
    CASE (3)
      count_part(4) = count_part(4) + 1  
    CASE (4)
      count_part(5) = count_part(5) + 1  
    CASE (5)
      count_part(6) = count_part(6) + 1  
    CASE (6)
      count_part(7) = count_part(7) + 1  
    CASE (7)
      count_part(8) = count_part(8) + 1  
    CASE (8)
      count_part(9) = count_part(9) + 1  
    CASE (9)
      count_part(10) = count_part(10) + 1  
    CASE (10)
      count_part(11) = count_part(11) + 1  
    CASE (11)
      count_part(12) = count_part(12) + 1  
    CASE (12)
      count_part(13) = count_part(13) + 1 
    END SELECT  
  END DO

  OPEN(unit=1,file="Number_particles_by_group.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file Number_particles_by_group.out"
  
  WRITE(1,*) "numero de particulas por grupo"
  DO i=1,13
    WRITE(1,*) "group " , i - 1, count_part(i) 
  END DO    
  
  CLOSE(1)
END SUBROUTINE create_contact_group

SUBROUTINE MSD_by_contact_number
  IMPLICIT NONE
  INTEGER :: i , j , k , kk , posi , status_open
  REAL(8) :: distance , norm  
  REAL(8), DIMENSION(nframes,13) :: MSD_Nc 
  REAL(8), DIMENSION(13) :: a 
  REAL(8), DIMENSION(nframes - 1,13) :: MSD_derivative
  MSD_Nc = 0.0
  
  IF (contacs_number_flag == 0) THEN
    !This function calculates Nc(j,i) that gives the number of contacts of the particle (j) at each frame (i) 
    CALL contacs_number
  END IF
  
  IF (create_contact_group_flag == 0) THEN
    CALL create_contact_group  !creates count_part that countains the number of particles by contact number
  END IF   
  
  MSD_by_contact_number_flag = 1
    
  OPEN(unit=11,file="MSD_by_contact_number.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file MSD_by_contact_number.out" 
  
  WRITE(11,*) "time simulation   ", &
  "group 0   ","group 1   ","group 2   ","group 3   ","group 4   ","group 5   ","group 6   ","group 7   ","group 8   ",&
  "group 9   ","group 10   ","group 11   ",  "group 12   "  
  
  !calculating correlation
  DO k = 1,nframes,MSD_Tgap  !delta t
    a = 0.0 
    DO kk = 1,nframes - k + 1       !number of origins
      DO j = 1,natoms
        SELECT CASE (index_part(j)) 
        CASE (0)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(1) = a(1) + distance 
        CASE (1)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(2) = a(2) + distance
        CASE (2)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(3) = a(3) + distance  
        CASE (3)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(4) = a(4) + distance
        CASE (4)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(5) = a(5) + distance 
        CASE (5)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(6) = a(6) + distance
        CASE (6)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(7) = a(7) + distance    
        CASE (7)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(8) = a(8) + distance 
        CASE (8)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(9) = a(9) + distance  
        CASE (9)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(10) = a(10) + distance 
        CASE (10)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(11) = a(11) + distance
        CASE (11)
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(12) = a(12) + distance 
        CASE (12)                                 
          distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
          distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
          distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
          a(13) = a(13) + distance
        END SELECT  
      END DO
    END DO
    
    DO j=1,13
      IF (count_part(j) /= 0) THEN    
        MSD_Nc(k,j) =  a(j) / (nframes - k + 1)  / count_part(j)
      END IF
    END DO
    
    WRITE(11,*) (k - 1) * dump_interval * delta_t ,MSD_Nc(k,1),MSD_Nc(k,2),MSD_Nc(k,3),MSD_Nc(k,4),MSD_Nc(k,5),MSD_Nc(k,6),&
                MSD_Nc(k,7),MSD_Nc(k,8), MSD_Nc(k,9),MSD_Nc(k,10),MSD_Nc(k,11),MSD_Nc(k,12),MSD_Nc(k,13)
    
  END DO
  !end calculating MSD (correlating)   
  
  OPEN(unit=13,file="MSD_derivative.out",status="replace",iostat = status_open )
  IF (status_open > 0) STOP "error opening MSD_derivative.out"
  
  WRITE(13,*) "time simulation   ", &
  "group 0   ","group 1   ","group 2   ","group 3   ","group 4   ","group 5   ","group 6   ","group 7   ","group 8   ",&
  "group 9   ","group 10   ","group 11   ",  "group 12   "  
  
  !Calculating rates of particle migration and rates of gaining a contact 
  DO i = 1,nframes - 1
    DO k = 1,13  
      MSD_derivative(i,k) = ( MSD_Nc(i+1,k) - MSD_Nc(i,k) ) / delta_t     
    END DO  
    
    WRITE(13,*) i,MSD_derivative(i,1),MSD_derivative(i,2),MSD_derivative(i,3),MSD_derivative(i,4),MSD_derivative(i,5),&
                          MSD_derivative(i,6),MSD_derivative(i,7),MSD_derivative(i,8),&
                          MSD_derivative(i,9),MSD_derivative(i,10),MSD_derivative(i,11),MSD_derivative(i,12),MSD_derivative(i,13)                      
   
  END DO   
  !END Calculating rates of particle migration and rates of gaining a contact 
  
  OPEN(unit=15,file="MSD_derivative_weighted.out",status="replace",iostat = status_open )
  IF (status_open > 0) STOP "error opening file MSD_derivative_weighted.out"
  
  WRITE(15,*) "time simulation   ", &
  "group 0   ","group 1   ","group 2   ","group 3   ","group 4   ","group 5   ","group 6   ","group 7   ","group 8   ",&
  "group 9   ","group 10   ","group 11   ",  "group 12   "  
  
  !Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size   
  DO i = 1,nframes - 1    
    WRITE(15,*) i,MSD_derivative(i,1)*count_part(1),MSD_derivative(i,2)*count_part(2),&
                          MSD_derivative(i,3)*count_part(3),MSD_derivative(i,4)*count_part(4),&
                          MSD_derivative(i,5)*count_part(5),MSD_derivative(i,6)*count_part(6),&
                          MSD_derivative(i,7)*count_part(7),MSD_derivative(i,8)*count_part(8),&
                          MSD_derivative(i,9)*count_part(9),MSD_derivative(i,10)*count_part(10),&
                          MSD_derivative(i,11)*count_part(11),MSD_derivative(i,12)*count_part(12),&
                          MSD_derivative(i,13)*count_part(13)
 
  END DO   
  !END Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size   
 
  
  !Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size , both normalized  
  OPEN(unit=17,file="MSD_derivative_weighted_normalized.out",status="replace",iostat = status_open )
  IF (status_open > 0) STOP "error opening file MSD_derivative_weighted_normalized.out"
  
  WRITE(17,*) "time simulation   ", &
  "group 0   ","group 1   ","group 2   ","group 3   ","group 4   ","group 5   ","group 6   ","group 7   ","group 8   ",&
  "group 9   ","group 10   ","group 11   ",  "group 12   "  
  
  !Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size   
  DO i = 1,nframes - 1    
    
    norm = 0.0
    DO j = 1,13
      norm = norm + MSD_derivative(i,j)*count_part(j)  
    END DO    
      
    WRITE(17,*)"t[",i,"]",MSD_derivative(i,1)*count_part(1)/norm,MSD_derivative(i,2)*count_part(2)/norm,&
                          MSD_derivative(i,3)*count_part(3)/norm,MSD_derivative(i,4)*count_part(4)/norm,&
                          MSD_derivative(i,5)*count_part(5)/norm,MSD_derivative(i,6)*count_part(6)/norm,&
                          MSD_derivative(i,7)*count_part(7)/norm,MSD_derivative(i,8)*count_part(8)/norm,&
                          MSD_derivative(i,9)*count_part(9)/norm,MSD_derivative(i,10)*count_part(10)/norm,&
                          MSD_derivative(i,11)*count_part(11)/norm,MSD_derivative(i,12)*count_part(12)/norm,&
                          MSD_derivative(i,13)*count_part(13)/norm
  
  END DO   
  !END Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size, both normalized  

  CLOSE(11)
  CLOSE(13)
  CLOSE(15)
  CLOSE(17)
END SUBROUTINE MSD_by_contact_number

SUBROUTINE Nc_by_contact_number
  IMPLICIT NONE
  INTEGER :: i , j , k , kk 
  REAL(8) :: norm   
  REAL(8), DIMENSION(13) :: a 
  REAL(8), DIMENSION(nframes,13) :: NC_ave
  REAL(8), DIMENSION(nframes - 1,13) :: Nc_ave_derivative
  NC_ave = 0.0
  
  Nc_by_contact_number_flag = 1
  
  IF (contacs_number_flag == 0) THEN
    !This function calculates Nc(j,i) that gives the number of contacts of the particle (j) at each frame (i) 
    CALL contacs_number
  END IF
  
  IF (create_contact_group_flag == 0) THEN
    CALL create_contact_group  !creates count_part that countains the number of particles by contact number
  END IF   
  
  OPEN(unit=12,file="Nc_by_contact_number.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file Nc_average.out"
  
  WRITE(12,*) "time simulation   ", &
  "group 0   ","group 1   ","group 2   ","group 3   ","group 4   ","group 5   ","group 6   ","group 7   ","group 8   ",&
  "group 9   ","group 10   ","group 11   ",  "group 12   " 
  
  DO j = 1,13
    a(j) = (j - 1)
  END DO    
  
  !calculating number of contacts average
  DO k = 1,nframes  !delta t 
    DO j = 1,natoms    
      SELECT CASE (index_part(j))
      CASE (0)
        a(1) = a(1) + Nc(j,k) 
      CASE (1)
        a(2) = a(2) + Nc(j,k)
      CASE (2)
        a(3) = a(3) + Nc(j,k)  
      CASE (3)
        a(4) = a(4) + Nc(j,k)
      CASE (4)
        a(5) = a(5) + Nc(j,k) 
      CASE (5)
        a(6) = a(6) + Nc(j,k) 
      CASE (6)
        a(7) = a(7) + Nc(j,k) 
      CASE (7)
        a(8) = a(8) + Nc(j,k) 
      CASE (8)
        a(9) = a(9) + Nc(j,k) 
      CASE (9)
        a(10) = a(10) + Nc(j,k) 
      CASE (10)
        a(11) = a(11) + Nc(j,k) 
      CASE (11)
        a(12) = a(12) + Nc(j,k) 
      CASE (12)                                 
        a(13) = a(13) + Nc(j,k) 
      END SELECT
    END DO
    
    DO j=1,13
      IF (count_part(j) /= 0.0) THEN    
        Nc_ave(k,j) =  a(j) / k  / count_part(j) !average over number of frames taken into account and number of members of the contact groups
      END IF 
    END DO
    
    WRITE(12,*) (k - 1) * dump_interval * delta_t ,Nc_ave(k,1),Nc_ave(k,2),Nc_ave(k,3),Nc_ave(k,4),Nc_ave(k,5),Nc_ave(k,6),&
                Nc_ave(k,7),Nc_ave(k,8), Nc_ave(k,9),Nc_ave(k,10),Nc_ave(k,11),Nc_ave(k,12),Nc_ave(k,13)
  END DO
  !End calculating number of contacts average
  
  OPEN(unit=14,file="Nc_ave_derivative.out",status="replace",iostat = status_open )
  IF (status_open > 0) STOP "error opening file Nc_ave_derivative.out"
  
  WRITE(14,*) "time simulation   ", &
  "group 0   ","group 1   ","group 2   ","group 3   ","group 4   ","group 5   ","group 6   ","group 7   ","group 8   ",&
  "group 9   ","group 10   ","group 11   ",  "group 12   " 
  
  !Calculating rates of particle migration and rates of gaining a contact 
  DO i = 1,nframes - 1
    DO k = 1,13  
      Nc_ave_derivative(i,k) = ( Nc_ave(i+1,k) - Nc_ave(i,k) ) / delta_t      
    END DO  
                          
    WRITE(14,*) i,Nc_ave_derivative(i,1),Nc_ave_derivative(i,2),Nc_ave_derivative(i,3),Nc_ave_derivative(i,4),&
                          Nc_ave_derivative(i,5),Nc_ave_derivative(i,6),Nc_ave_derivative(i,7),Nc_ave_derivative(i,8),&
                          Nc_ave_derivative(i,9),Nc_ave_derivative(i,10),Nc_ave_derivative(i,11),Nc_ave_derivative(i,12),&
                          Nc_ave_derivative(i,13)                      
   
  END DO   
  !END Calculating rates of particle migration and rates of gaining a contact
  
  OPEN(unit=16,file="Nc_ave_derivative_weighted.out",status="replace",iostat = status_open )
  IF (status_open > 0) STOP "error opening file Nc_ave_derivative_weighted.out"
  
  WRITE(16,*) "time simulation   ", &
  "group 0   ","group 1   ","group 2   ","group 3   ","group 4   ","group 5   ","group 6   ","group 7   ","group 8   ",&
  "group 9   ","group 10   ","group 11   ",  "group 12   " 
  
  !Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size   
  DO i = 1,nframes - 1    
                       
    WRITE(16,*) i,Nc_ave_derivative(i,1)*count_part(1),Nc_ave_derivative(i,2)*count_part(2),&
                          Nc_ave_derivative(i,3)*count_part(3),Nc_ave_derivative(i,4)*count_part(4),&
                          Nc_ave_derivative(i,5)*count_part(5),Nc_ave_derivative(i,6)*count_part(6),&
                          Nc_ave_derivative(i,7)*count_part(7),Nc_ave_derivative(i,8)*count_part(8),&
                          Nc_ave_derivative(i,9)*count_part(9),Nc_ave_derivative(i,10)*count_part(10),&
                          Nc_ave_derivative(i,11)*count_part(11),Nc_ave_derivative(i,12)*count_part(12),&
                          Nc_ave_derivative(i,13)*count_part(13)                      
   
  END DO   
  !END Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size   
   
  OPEN(unit=18,file="Nc_ave_derivative_weighted_normalized.out",status="replace",iostat = status_open )
  IF (status_open > 0) STOP "error opening file 18"
  
  WRITE(18,*) "time simulation   ", &
  "group 0   ","group 1   ","group 2   ","group 3   ","group 4   ","group 5   ","group 6   ","group 7   ","group 8   ",&
  "group 9   ","group 10   ","group 11   ",  "group 12   " 
  
  !Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size   
  DO i = 1,nframes - 1    
                        
    norm = 0.0
    DO j = 1,13
      norm = norm + Nc_ave_derivative(i,j)*count_part(j)  
    END DO                      
                          
    WRITE(18,*) i,Nc_ave_derivative(i,1)*count_part(1)/norm,Nc_ave_derivative(i,2)*count_part(2)/norm,&
                          Nc_ave_derivative(i,3)*count_part(3)/norm,Nc_ave_derivative(i,4)*count_part(4)/norm,&
                          Nc_ave_derivative(i,5)*count_part(5)/norm,Nc_ave_derivative(i,6)*count_part(6)/norm,&
                          Nc_ave_derivative(i,7)*count_part(7)/norm,Nc_ave_derivative(i,8)*count_part(8)/norm,&
                          Nc_ave_derivative(i,9)*count_part(9)/norm,Nc_ave_derivative(i,10)*count_part(10)/norm,&
                          Nc_ave_derivative(i,11)*count_part(11)/norm,Nc_ave_derivative(i,12)*count_part(12)/norm,&
                          Nc_ave_derivative(i,13)*count_part(13)/norm                      
   
  END DO   
  !END Calculating rates of particle migration weighted by group size and rates of gaining a contact weighted by group size, both normalized  

  CLOSE(12)
  CLOSE(14) 
  CLOSE(16)
  CLOSE(18)
END SUBROUTINE Nc_by_contact_number    


SUBROUTINE MSD_total
  IMPLICIT NONE
  INTEGER :: i , j , k , kk 
  REAL(8) :: distance , sum_msd, drx, dry, drz
  REAL(8), DIMENSION(size(MSD)) :: time_msd 
  REAL(8), DIMENSION(size(MSD) - 2) :: derivative_y 
  
  
  OPEN(unit=10,file="MSD_total.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file MSD_total.out"
  
  WRITE(10,'(3x,a5,13x,a3)') "time","MSD"
  
  DO k = 1,nframes,MSD_Tgap  !delta t
    sum_msd = 0.0  
    
!    !calculating autocorrelation function
!    DO kk = 1,nframes - k + 1       !number of origins
!      DO j = 1,natoms    
!        distance = ( rxyz(j,1,k + kk - 1) -  rxyz(j,1,kk) ) ** 2.0
!        distance = distance + ( rxyz(j,2,k + kk - 1) -  rxyz(j,2,kk) ) ** 2.0
!        distance = distance + ( rxyz(j,3,k + kk - 1) -  rxyz(j,3,kk) ) ** 2.0    
!        sum_msd = sum_msd + distance 
!      END DO
!    END DO
!    !End calculating autocorrelation function
    
    !NEW MSD CALCULATION
    DO j = 1,natoms    
      drx = rxyz(j,1,k) - rxyz(j,1,1)
      dry = rxyz(j,2,k) - rxyz(j,2,1)
      drz = rxyz(j,3,k) - rxyz(j,3,1)
      sum_msd = sum_msd + (drx**2.0 + dry**2.0 + drz**2.0)    
    END DO    
      
!    MSD(k) =  sum_msd / (nframes - k + 1)  / natoms
    MSD(k) =  sum_msd / natoms
    time_msd(k) = (k - 1) * dump_interval * delta_t 
    WRITE(10,*) time_msd(k) , MSD(k)
    WRITE(*,*) "MSD calculation, timestep  ",k
  END DO
  
  CALL  derivative(MSD,time_msd,derivative_y)
  
  WRITE(10,*) 
  WRITE(10,*) 
  WRITE(10,*) "time" , "3D Diffusion coefficient"
  
  DO k = 1,size(MSD) - 2
    WRITE(10,*) time_msd(k) , derivative_y(k) / 6.0      
  END DO    
  
  CLOSE(10)
  
END SUBROUTINE MSD_total


SUBROUTINE collective_intermediate_scattering_function 
  IMPLICIT NONE
  INTEGER :: i , j , jj , k , kk   
  COMPLEX(8) :: a , b , c
  
  OPEN(unit=6,file="Fc.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file Fc.out"
  
  WRITE(6,'(3x,a5,13x,a3,23x,a3,23x,a3)') "time","Fcx","Fcy","Fcz"
  
  DO k = 1,nframes,Fc_Tgap  !delta t
    a = (0.0,0.0)  
    b = (0.0,0.0)  
    c = (0.0,0.0)
    
    !calculating autocorrelation function
    DO kk = 1,nframes - k + 1
      DO j = 1,natoms      
        DO jj = 1,natoms  
           a =  a + exp(-1.0 * (0.0, 1.0) * Fc_q * ( rxyz(j,1,k + kk - 1)  - rxyz(jj,1,kk) ) )
           b =  b + exp(-1.0 * (0.0, 1.0) * Fc_q * ( rxyz(j,2,k + kk - 1)  - rxyz(jj,2,kk) ) )
           c =  c + exp(-1.0 * (0.0, 1.0) * Fc_q * ( rxyz(j,3,k + kk - 1)  - rxyz(jj,3,kk) ) )
        END DO 
      END DO
    END DO
    !End calculating autocorrelation function
      
    Fsx(k) =  abs(a) / (nframes - k + 1)  / natoms
    Fsy(k) =  abs(b) / (nframes - k + 1)  / natoms
    Fsz(k) =  abs(c) / (nframes - k + 1)  / natoms
    WRITE(6,*) (k - 1) * dump_interval * delta_t , Fsx(k) , Fsy(k) , Fsz(k)
  END DO

  CLOSE(6)
END SUBROUTINE collective_intermediate_scattering_function

!SUBROUTINE self_intermediate_scattering_function_FFT 
!  IMPLICIT NONE
!  INTEGER :: i , j , k , kk , dim = 3 , m , n , count
!  REAL(8) :: L , aux1 , aux2 , box_lower_limit , box_upper_limit , q_max , dq , q_max_peak 
!  REAL(8), DIMENSION(natoms) :: x_ave , y_ave  , z_ave  
!  COMPLEX(8), DIMENSION(nframes,natoms) :: rho_k_t_x  , rho_k_t_y , rho_k_t_z  
!  COMPLEX(8), DIMENSION(2 * nframes) :: correl_x , correl_y , correl_z   
!  
!  OPEN(unit=1,file=input_file,status="old",iostat = status_open ) 
!  IF (status_open > 0) STOP "error opening file 1"
!  READ(1,*)
!  READ(1,*)
!  READ(1,*)
!  READ(1,*) 
!  READ(1,*)
!  READ(1,*) box_lower_limit , box_upper_limit
!  READ(1,*)
!  READ(1,*)
!  READ(1,*)
!  
!  L = box_upper_limit - box_lower_limit
!  rho = natoms / (L * L * L)
!  q_max_peak = 2.0
!  
!  
!  OPEN(unit=7,file="Fs_fft.out",status="replace",iostat = status_open ) 
!  IF (status_open > 0) STOP "error opening file 7"
!  
!  !getting all atoms positions
!  DO i = 1,nframes
!    DO j = 1,natoms       
!        READ(1,*) aux1 , aux2 , rxyz(j,1,i) , rxyz(j,2,i) , rxyz(j,3,i)    
!    END DO      
!    
!    
!    IF (i < nframes) THEN 
!      DO j = 1,9
!        READ(1,*)
!      END DO
!    END IF
!    
!  END DO 
!  !End getting all atoms positions
!  
!  
!  !calculating rho(k,t,i)
!  DO i = 1,nframes  
!    DO j = 1,natoms  
!      rho_k_t_x(i,j) = exp(-1.0 * (0.0, 1.0) * q_max_peak * rxyz(j,1,i) ) 
!      rho_k_t_y(i,j) = exp(-1.0 * (0.0, 1.0) * q_max_peak * rxyz(j,2,i) ) 
!      rho_k_t_z(i,j) = exp(-1.0 * (0.0, 1.0) * q_max_peak * rxyz(j,3,i) ) 
!    END DO 
!  END DO
!  !End calculating rho(k,t,i)
!  
!  Fsx = 0.0
!  Fsy = 0.0
!  Fsz = 0.0
!  
!  DO j = 1,natoms
!    correl_x = 0.0
!    correl_y = 0.0
!    correl_z = 0.0
!    
!    DO i = 1,nframes        
!      correl_x(i) = rho_k_t_x(i,j)
!      correl_y(i) = rho_k_t_y(i,j) 
!      correl_z(i) = rho_k_t_z(i,j)
!    END DO  
!    
!  
!    !fourier transforming rho(k,t) in rho(k,frequency)
!    CALL fft(correl_x)
!    CALL fft(correl_y)
!    CALL fft(correl_z)
!    !End fourier transforming rho(k,t) in rho(k,frequency)
!   
!    !calculating correlations in frequency
!     correl_x = correl_x * conjg(correl_x)
!     correl_y = correl_y * conjg(correl_y)
!     correl_z = correl_z * conjg(correl_z)
!     !End calculating correlations in frequency
!   
!     !calculating correlations in time
!     CALL ifft(correl_x)
!     CALL ifft(correl_y)
!     CALL ifft(correl_z)
!     !End calculating correlations in time
!   
!     DO i = 1,nframes    !normalizing results    
!       Fsx(i) = Fsx(i) + (1.0 / (nframes - i + 1) ) * correl_x(i) 
!       Fsy(i) = Fsy(i) + (1.0 / (nframes - i + 1) ) * correl_y(i) 
!       Fsz(i) = Fsz(i) + (1.0 / (nframes - i + 1) ) * correl_z(i) 
!     END DO  
!  END DO
!  
!  DO j = 1,nframes
!    WRITE(7,*) "t[",j,"]" , Fsx(j) / natoms , Fsy(j) / natoms , Fsz(j) / natoms 
!  END DO
!
!    
!  CLOSE(1) 
!  CLOSE(7)
!END SUBROUTINE self_intermediate_scattering_function_FFT


SUBROUTINE self_intermediate_scattering_function 
  IMPLICIT NONE
  INTEGER :: i , j , k , kk   
  COMPLEX(8) :: a , b , c
  
  OPEN(unit=5,file="Fs.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file Fs.out"
  
  WRITE(5,'(3x,a5,13x,a3,23x,a3,23x,a3)') "time","Fsx","Fsy","Fsz"     
  
  !this is a correlation function for each particle (averaged)
  
  DO k = 1,nframes,Fs_Tgap  !delta t, Tgap for selecting time origins
    a = (0.0,0.0)  
    b = (0.0,0.0)  
    c = (0.0,0.0)
     
    DO kk = 1,nframes - k + 1 !this is tmax, it is the number of origins used for calculating correlation
      DO j = 1,natoms      
        a =  a + exp(-1.0 * (0.0, 1.0) * Fs_q * ( rxyz(j,1,k + kk - 1)  - rxyz(j,1,kk) ) )
        b =  b + exp(-1.0 * (0.0, 1.0) * Fs_q * ( rxyz(j,2,k + kk - 1)  - rxyz(j,2,kk) ) )
        c =  c + exp(-1.0 * (0.0, 1.0) * Fs_q * ( rxyz(j,3,k + kk - 1)  - rxyz(j,3,kk) ) )     
      END DO
    END DO
    
    !averaging in tmax and in natoms
    Fsx(k) =  abs(a) / (nframes - k + 1) / natoms  
    Fsy(k) =  abs(b) / (nframes - k + 1) / natoms
    Fsz(k) =  abs(c) / (nframes - k + 1) / natoms
    WRITE(5,*) (k - 1) * dump_interval * delta_t , Fsx(k) , Fsy(k) , Fsz(k)
  END DO

  CLOSE(5)
END SUBROUTINE self_intermediate_scattering_function

SUBROUTINE static_factor_1 
  IMPLICIT NONE
  INTEGER :: i , j , k , kk , m , n , count , pos
  COMPLEX(8) :: a , b , c 
  REAL(8) :: dq, dx, dy, dz   
  
  dq = 2 * pi / L 
  
  OPEN(unit=3,file="S(q).out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file S(q).out"
  
  OPEN(unit=4,file="Ls.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file Ls.out"
  
  WRITE(4,*) "time (simu. units)   ","Ls"
  
  DO i = 1,nframes
      
    WRITE(3,*) "timestep    ",i - 1
    
    
    !calculating S(q)
    DO j = 1,number_of_bins_sq 
      q(j) = j * dq
      
      
      a = (0.0,0.0)  
      b = (0.0,0.0)  
      c = (0.0,0.0)
      DO k = 1,natoms
        DO kk = 1,natoms 
          dx = rxyz(kk,1,i) - rxyz(k,1,i)
          dx = dx - L*nint(dx/L)
          dy = rxyz(kk,2,i) - rxyz(k,2,i) 
          dy = dy - L*nint(dy/L)
          dz = rxyz(kk,3,i) - rxyz(k,3,i)  
          dz = dz - L*nint(dz/L)  
          a =  a + exp( (0.0, 1.0) * q(j) * dx )
          b =  b + exp( (0.0, 1.0) * q(j) * dy )
          c =  c + exp( (0.0, 1.0) * q(j) * dz ) 
        END DO
      END DO
      
      S_qx(j) =  abs(a) / natoms
      S_qy(j) =  abs(b) / natoms
      S_qz(j) =  abs(c) / natoms
      S_q(j) = ( S_qx(j) + S_qy(j) + S_qz(j) ) / 3.0
      WRITE(3,*) q(j) , S_qx(j) , S_qy(j) , S_qz(j) , S_q(j) 
      
      S_qx_ave(j) = S_qx_ave(j) + S_qx(j)
      S_qy_ave(j) = S_qy_ave(j) + S_qy(j)
      S_qz_ave(j) = S_qz_ave(j) + S_qz(j)
      S_q_ave(j) = S_q_ave(j) + S_q(j)
            
    END DO  
    !end calculating S(q)
    
    CALL ordenar_cres(S_q,pos)
    Ls(i) = 2 * pi / q(pos) / radius ! Ls is dimensionless
    WRITE(3,*) "Ls = ", Ls(i)
    
    WRITE(4,*) real(i) * delta_t * dump_interval,Ls(i)
    
    WRITE(3,*)
    WRITE(3,*)
    
    WRITE(*,*) "S(q) of frame" ,i - 1 , "Done" 
    
  END DO  
  
  WRITE(3,*) "S(q) averaged  "
  DO j = 1, number_of_bins_sq
    S_qx_ave(j) = S_qx_ave(j) / nframes
    S_qy_ave(j) = S_qy_ave(j) / nframes
    S_qz_ave(j) = S_qz_ave(j) / nframes
    S_q_ave(j) = S_q_ave(j) / nframes
    WRITE(3, *) q(j), S_qx_ave(j), S_qy_ave(j), S_qz_ave(j), S_q_ave(j)
  END DO

  CLOSE(3)
  CLOSE(4)
END SUBROUTINE static_factor_1 

SUBROUTINE fitting_gr
IMPLICIT NONE
INTEGER :: i , j , k ,pos
REAL(8), DIMENSION(number_of_bins_gr) :: x ,y , y2 
REAL(8):: h , aux

IF (radial_distribution_flag == 0) THEN
  CALL radial_distri
END IF

OPEN(unit=2,file="gr_fitted.out",status="replace",iostat = status_open ) 
IF (status_open > 0) STOP "error opening file gr_fitted.out"


DO j = 1,number_of_bins_gr 
  x(j) = distance_rdf(j)
END DO

h = (x(number_of_bins_gr - 1) - x(1)) / n_interval

DO j = 1,n_interval+1
  distance_rdf_fitted(j) = x(1)+(j-1)*h
END DO 

DO i = 1,nframes 
  
  DO j = 1,number_of_bins_gr     
    y(j) = rdf(j,i)
    IF (rdf(j,i) > 0.0) THEN
      pos = j
      EXIT
    END IF    
  END DO
  
  DO j = pos+1,number_of_bins_gr     
    y(j) = rdf(j,i)    
  END DO
  
  CALL  spline(x,y,10e31_dp,10e31_dp,y2) 
  
  DO j = 1,n_interval+1
    rdf_fitted(j,i) = splint(x,y,y2,distance_rdf_fitted(j))
    IF ( ( distance_rdf_fitted(j) < distance_rdf(pos) ) .or. (rdf_fitted(j,i) < 0.0)) THEN 
      rdf_fitted(j,i) = 0.0  
    END IF    
  END DO 
 
END DO  

DO i = 1,nframes
    
  IF (i == 1 .or. i == nframes .or. (mod(i-1,gr_Tgap) == 0) ) THEN
  !writing gr for one frame
  WRITE(2,*) " timestep ", i - 1
  WRITE(2,*) " distance                   g(r)"  
  
  DO j = 1,n_interval+1 
    WRITE(2,*) distance_rdf_fitted(j) , rdf_fitted(j,i)
  END DO
  WRITE(2,*)
  WRITE(2,*)
  !end writing gr for one frame
  END IF 
  
END DO    
  
CLOSE(2)
END SUBROUTINE fitting_gr   

SUBROUTINE static_factor_2 
  IMPLICIT NONE
  INTEGER :: i , j , k , m , n , count , pos
  REAL(8), DIMENSION(number_of_bins_gr) :: y
  REAL(8) :: dq , aux
  
  dq = 2 * pi / L 
  
  IF (radial_distribution_flag == 0) THEN
    CALL radial_distri
  END IF
  
  OPEN(unit=3,file="S(q)_2.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file S(q)_2.out"
  
  OPEN(unit=4,file="Ls.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file Ls.out"
  
   WRITE(4,*) "  time (simu. units)   ","     Ls", "                       Peak height"
  
  DO i = 1,nframes
    
    WRITE(3,*) "timestep    ",i - 1
        
    DO j = 1,number_of_bins_sq 
      q(j) = j * dq
      
      !making vector to integrate r²g(r)sin(qr)/qr
      DO k = 1,number_of_bins_gr 
     
        y(k)  =  distance_rdf(k) * (rdf(k,i) - 1.0 ) * (sin( q(j) * distance_rdf(k) ) ) 
        
      END DO    
     
      S_q(j) = 1 + 4.0 * pi * rho * trapezios(y,distance_rdf,number_of_bins_gr) / q(j)
      S_q_ave(j) = S_q_ave(j) + S_q(j)
     
      WRITE(3,*) q(j) , S_q(j) 
      
    END DO     
    
    CALL ordenar_cres(S_q,pos)
    Ls(i) = 2 * pi / q(pos) / radius ! Ls is dimensionless
    WRITE(3,*) "Ls = ", Ls(i)
    
    WRITE(4,*) real(i-1) * delta_t * dump_interval,Ls(i), S_q(number_of_bins_sq)
    
    WRITE(3,*)
    WRITE(3,*) 
    
    WRITE(*,*) "S(q) of frame" ,i - 1 , "Done" 
    
  END DO 
  
  WRITE(3,*) "S(q) averaged  "
  DO j = 1,number_of_bins_sq 
    S_q_ave(j) = S_q_ave(j) / nframes
    WRITE(3,*) q(j) , S_q_ave(j)
  END DO  
  
  CLOSE(3)
  CLOSE(4)
END SUBROUTINE static_factor_2

SUBROUTINE static_factor_3
  IMPLICIT NONE
  INTEGER :: i , j , k 
  REAL(8), DIMENSION(number_of_bins_gr) :: y
  REAL(8) :: dq , aux
  
  dq = 2 * pi / L
  
  DO j = 1,number_of_bins_sq 
    q(j) = j * dq   
  END DO
  
  IF (radial_distribution_flag == 0) THEN
    CALL radial_distri
  END IF
  
  OPEN(unit=3,file="S(q)_3.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening file S(q)_3.out"

      
    DO k = 1,number_of_bins_gr     
      y(k)  =  rdf_ave(k)       
    END DO  
        
    CALL realft_dp(y,1)  
    
    S_q = 1.0 + rho * y
    
!    WRITE(3,*) "timestep    ",i - 1
    
    DO j=1,number_of_bins_sq
      WRITE(3,*) q(j) , S_q(j) 
    END DO     

    WRITE(3,*)
    WRITE(3,*) 
    
!    S_q_ave = S_q_ave + S_q 
!    
!
!  
!  WRITE(3,*) "S(q) averaged"
!  DO j = 1,number_of_bins_sq 
!    S_q_ave(j) = S_q_ave(j) / ave_number_Sq
!    WRITE(3,*) q(j) , S_q_ave(j)
!  END DO
!  
  CLOSE(3)
!    
!  IMPLICIT NONE
!  INTEGER :: i , j , k
!  REAL(8) :: dq
!  COMPLEX(8) , DIMENSION (number_of_bins_sq,3,nframes) :: rho_q , rho_q_negative
!  COMPLEX(8) , DIMENSION (number_of_bins_sq) :: aux_x, aux_y, aux_z , aux_x_negative, aux_y_negative, aux_z_negative  
!  
!  rho_q = 0.0
!  rho_q_negative = 0.0
!  
!  dq = 2 * pi / L !cubic box all the directions have the same dq
!  
!  OPEN(unit=3,file="S(q)_3.out",status="replace",iostat = status_open ) 
!  IF (status_open > 0) STOP "error opening file S(q)_3.out"
!  
!  DO i = 1,nframes,Sq_Tgap! loop for q
!    WRITE(3,*) "timestep    ",i - 1
!    
!    
!    DO j = 1,number_of_bins_sq 
!      q(j) = j * dq
!      
!      DO k = 1, natoms 
!        rho_q(j,1,i) = rho_q(j,1,i) + exp(-(0.0, 1.0) * q(j) * rxyz(k,1,i))
!        rho_q_negative(j,1,i) = rho_q_negative(j,1,i) + exp((0.0, 1.0) * q(j) * rxyz(k,1,i))
!        
!        rho_q(j,2,i) = rho_q(j,2,i) + exp(-(0.0, 1.0) * q(j) * rxyz(k,2,i))
!        rho_q_negative(j,2,i) = rho_q_negative(j,2,i) + exp((0.0, 1.0) * q(j) * rxyz(k,2,i))
!        
!        rho_q(j,3,i) = rho_q(j,3,i) + exp(-(0.0, 1.0) * q(j) * rxyz(k,3,i))
!        rho_q_negative(j,3,i) = rho_q_negative(j,3,i) + exp((0.0, 1.0) * q(j) * rxyz(k,1,3))
!      END DO
!      
!      aux_x(j) = rho_q(j,1,i)
!      aux_x_negative(j) = rho_q_negative(j,1,i) 
!      aux_y(j) = rho_q(j,2,i)
!      aux_y_negative(j) = rho_q_negative(j,2,i) 
!      aux_z(j) = rho_q(j,3,i)
!      aux_z_negative(j) = rho_q_negative(j,3,i) 
!          
!    END DO 
!    
!    S_qx = 0.0
!    S_qy = 0.0
!    S_qz = 0.0
!    
!    CALL correl_complex(aux_x, aux_x_negative, S_qx)
!    S_qx = S_qx  / natoms
!    CALL correl_complex(aux_y, aux_y_negative, S_qy)
!    S_qy = S_qy  / natoms
!    CALL correl_complex(aux_z, aux_z_negative, S_qz)
!    S_qz = S_qz  / natoms
!    
!    DO j = 1,number_of_bins_sq 
!      WRITE(3,*) q(j) , S_qx(j) , S_qy(j), S_qz(j) , (S_qx(j) + S_qy(j) + S_qz(j)) / 3.0   
!      S_qx_ave(j) = S_qx_ave(j) + S_qx(j)
!      S_qy_ave(j) = S_qy_ave(j) + S_qy(j)
!      S_qz_ave(j) = S_qz_ave(j) + S_qz(j)
!    END DO  
!    
!    WRITE(3,*)
!    WRITE(3,*)
!  END DO
!  
!  WRITE(3,*) "S(q) averaged"
!  
!  DO j = 1,number_of_bins_sq 
!    S_qx_ave(j) = S_qx_ave(j) / nframes 
!    S_qy_ave(j) = S_qy_ave(j) / nframes
!    S_qz_ave(j) = S_qz_ave(j) / nframes  
!    S_q_ave(j) = (S_qx_ave(j) + S_qy_ave(j) + S_qz_ave(j)) / 3.0
!    WRITE(3,*) q(j) , S_qx_ave(j) , S_qy_ave(j) , S_qz_ave(j) , S_q_ave(j)
!  END DO
!  
!  CLOSE(3)
END SUBROUTINE static_factor_3

SUBROUTINE radial_distri
  IMPLICIT NONE
  INTEGER :: i , j , k , z , pico1, pico2, n_rdfs
  REAL(8) :: dr , dx , dy , dz
  REAL(8), DIMENSION(:), ALLOCATABLE :: y
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: y_rdf,rdf_ave_per_bead
  REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: rdf_types
  
  rdf_ave = 0.0 
  n_rdfs = 10 ! for max 4 atom types, 4+3+2+1 = 10 combinations  
  ALLOCATE(y(number_of_bins_gr),y_rdf(number_of_bins_gr,n_rdfs))
  ALLOCATE(rdf_types(number_of_bins_gr,nframes,n_rdfs))
  ALLOCATE(rdf_ave_per_bead(number_of_bins_gr,n_rdfs))
  rdf_ave_per_bead = 0.0
  rdf_types = 0.0
  
  radial_distribution_flag = 1
  
  dr = 0.5 * L / number_of_bins_gr
  
  OPEN(unit=2,file="rdf.out",status="replace",iostat = status_open) 
  IF (status_open > 0) STOP "error opening rdf file"
  
  DO i = 1,nframes
    
    !calculating gr for one frame
    DO j = 1,natoms-1 !loop for number of atoms
      DO k = j+1,natoms   
        dx = rxyz(j,1,i) - rxyz(k,1,i)
        dx = dx - L*nint(dx/L)
        dy = rxyz(j,2,i) - rxyz(k,2,i)
        dy = dy - L*nint(dy/L)
        dz = rxyz(j,3,i) - rxyz(k,3,i)   
        dz = dz - L*nint(dz/L)
        vetor_length = sqrt(dx*dx+dy*dy+dz*dz)
        IF (vetor_length < 0.5 * L) THEN
          ! total rdf  
          rdf(ceiling(vetor_length / dr),i) = rdf(ceiling(vetor_length / dr),i) + 2
          
          ! rdf for beads
          IF (particle_type(j,i) == 1 .and. particle_type(k,i) == 1) THEN
            rdf_types(ceiling(vetor_length / dr),i,1) = rdf_types(ceiling(vetor_length / dr),i,1) + 2
          ELSE IF ((particle_type(j,i) == 1 .and. particle_type(k,i) == 2) .or. (particle_type(j,i) == 2 .and. &
              particle_type(k,i) == 1)) THEN
            rdf_types(ceiling(vetor_length / dr),i,2) = rdf_types(ceiling(vetor_length / dr),i,2) + 1          
          ELSE IF ((particle_type(j,i) == 1 .and. particle_type(k,i) == 3) .or. (particle_type(j,i) == 3 .and. &
              particle_type(k,i) == 1)) THEN    
            rdf_types(ceiling(vetor_length / dr),i,3) = rdf_types(ceiling(vetor_length / dr),i,3) + 1            
          ELSE IF ((particle_type(j,i) == 1 .and. particle_type(k,i) == 4) .or. (particle_type(j,i) == 4 .and. &
              particle_type(k,i) == 1)) THEN              
            rdf_types(ceiling(vetor_length / dr),i,4) = rdf_types(ceiling(vetor_length / dr),i,4) + 1               
          ELSE IF ((particle_type(j,i) == 2 .and. particle_type(k,i) == 2)) THEN
            rdf_types(ceiling(vetor_length / dr),i,5) = rdf_types(ceiling(vetor_length / dr),i,5) + 2               
          ELSE IF ((particle_type(j,i) == 2 .and. particle_type(k,i) == 3) .or. (particle_type(j,i) == 3 .and. &
              particle_type(k,i) == 2)) THEN              
            rdf_types(ceiling(vetor_length / dr),i,6) = rdf_types(ceiling(vetor_length / dr),i,6) + 1               
          ELSE IF ((particle_type(j,i) == 2 .and. particle_type(k,i) == 4) .or. (particle_type(j,i) == 4 .and. &
              particle_type(k,i) == 2)) THEN              
            rdf_types(ceiling(vetor_length / dr),i,7) = rdf_types(ceiling(vetor_length / dr),i,7) + 1               
          ELSE IF ((particle_type(j,i) == 3 .and. particle_type(k,i) == 3)) THEN
            rdf_types(ceiling(vetor_length / dr),i,8) = rdf_types(ceiling(vetor_length / dr),i,8) + 2               
          ELSE IF ((particle_type(j,i) == 3 .and. particle_type(k,i) == 4) .or. (particle_type(j,i) == 4 .and. &
              particle_type(k,i) == 3)) THEN              
            rdf_types(ceiling(vetor_length / dr),i,9) = rdf_types(ceiling(vetor_length / dr),i,9) + 1               
          ELSE IF ((particle_type(j,i) == 4 .and. particle_type(k,i) == 4)) THEN
            rdf_types(ceiling(vetor_length / dr),i,10) = rdf_types(ceiling(vetor_length / dr),i,10) + 2               
          ELSE                
          END IF    
          
        END IF
      END DO  
    END DO 
    
    DO z = 1,number_of_bins_gr
      distance_rdf(z) = (z-0.5)*dr   
      rdf(z,i) =  rdf(z,i) / ( (natoms - 1)  * 4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho)
      IF (rho_rdf(1)/= 0.0) rdf_types(z,i,1) =  rdf_types(z,i,1) / ( (n_atoms_rdf(1) - 1)  * &
      4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(1))       
      IF (rho_rdf(2)/= 0.0) rdf_types(z,i,2) =  rdf_types(z,i,2) * ((n_atoms_rdf(1) - 1) + (n_atoms_rdf(2) - 1)) / &
      ( ( (n_atoms_rdf(1) - 1) * (n_atoms_rdf(2) - 1) )  * 4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(2)) 
      IF (rho_rdf(3)/= 0.0) rdf_types(z,i,3) =  rdf_types(z,i,3) * ((n_atoms_rdf(1) - 1) + (n_atoms_rdf(3) - 1)) / &
      ( ( (n_atoms_rdf(1) - 1) * (n_atoms_rdf(3) - 1) )  * 4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(3))        
      IF (rho_rdf(4)/= 0.0) rdf_types(z,i,4) =  rdf_types(z,i,4) * ((n_atoms_rdf(1) - 1) + (n_atoms_rdf(4) - 1)) / &
      ( ( (n_atoms_rdf(1) - 1) * (n_atoms_rdf(4) - 1) )  * 4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(4))      
      IF (rho_rdf(5)/= 0.0) rdf_types(z,i,5) =  rdf_types(z,i,5) / ( (n_atoms_rdf(2) - 1)  * &
      4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(5))       
      IF (rho_rdf(6)/= 0.0) rdf_types(z,i,6) =  rdf_types(z,i,6) * ((n_atoms_rdf(2) - 1) + (n_atoms_rdf(3) - 1)) / &
      ( ( (n_atoms_rdf(2) - 1) * (n_atoms_rdf(3) - 1) )  * 4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(6))      
      IF (rho_rdf(7)/= 0.0) rdf_types(z,i,7) =  rdf_types(z,i,7) * ((n_atoms_rdf(2) - 1) + (n_atoms_rdf(4) - 1)) / &
      ( ( (n_atoms_rdf(2) - 1) * (n_atoms_rdf(4) - 1) )  * 4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(7))       
      IF (rho_rdf(8)/= 0.0) rdf_types(z,i,8) =  rdf_types(z,i,8) / ( (n_atoms_rdf(3) - 1)  * &
      4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(8))       
      IF (rho_rdf(9)/= 0.0) rdf_types(z,i,9) =  rdf_types(z,i,9) * ((n_atoms_rdf(3) - 1) + (n_atoms_rdf(4) - 1)) / &
      ( ( (n_atoms_rdf(3) - 1) * (n_atoms_rdf(4) - 1) )  * 4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(9))       
      IF (rho_rdf(10)/= 0.0) rdf_types(z,i,10) =  rdf_types(z,i,10) / ( (n_atoms_rdf(4) - 1)  * &
      4.0 / 3.0 * pi * ((z) ** 3 - (z - 1) ** 3) * (dr ** 3) * rho_rdf(10))       
    END DO
    !end calculating gr for one frame
    
    IF (i == 1 .or. i == nframes .or. (mod(i-1,gr_Tgap) == 0) ) THEN
    y = 0.0    
    y_rdf = 0.0
    !writing gr for one frame
    WRITE(2,*) " timestep ", i - 1
  WRITE(2,'(1x,a8,2x,a10,2x,a8,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,&
  a8,2x,a6,2x,a8,2x,a6)') "distance","g(r)_total","Nc_total","g(r)_1-1","Nc_1-1","g(r)_1-2","Nc_1-2",  &
  "g(r)_1-3","Nc_1-3","g(r)_1-4","Nc_1-4","g(r)_2-2","Nc_2-2","g(r)_2-3","Nc_2-3","g(r)_2-4","Nc_2-4","g(r)_3-3","Nc_3-3", &
  "g(r)_3-4","Nc_3-4","g(r)_4-4","Nc_4-4"
    DO j = 1,number_of_bins_gr
      y(j) =  rdf(j,i) * distance_rdf(j) * distance_rdf(j)  
      y_rdf(j,:) = rdf_types(j,i,:) * distance_rdf(j) * distance_rdf(j) 
      WRITE(2,*) distance_rdf(j) , rdf(j,i), 4.0 * pi * rho * trapezios(y,distance_rdf,number_of_bins_gr), &
                 (rdf_types(j,i,z), 4.0 * pi * rho_rdf(z) * trapezios(y_rdf(:,z),distance_rdf,number_of_bins_gr) , &
                 z=1,10)
    END DO
    WRITE(2,*)
    WRITE(2,*)
    !end writing gr for one frame
    END IF    
    
    WRITE(*,*) "RDF calculation, timestep  ",i
    
  END DO  
  
  DO i = 1,nframes
   DO j = 1,number_of_bins_gr
      rdf_ave(j) = rdf_ave(j) + rdf(j,i)
      rdf_ave_per_bead(j,:) = rdf_ave_per_bead(j,:) + rdf_types(j,i,:)
   END DO
  END DO
  
  WRITE(2,*) "RDF(q) averaged  "
  WRITE(2,'(1x,a8,2x,a10,2x,a8,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,a8,2x,a6,2x,&
  a8,2x,a6,2x,a8,2x,a6)') "distance","g(r)_total","Nc_total","g(r)_1-1","Nc_1-1","g(r)_1-2","Nc_1-2",  &
  "g(r)_1-3","Nc_1-3","g(r)_1-4","Nc_1-4","g(r)_2-2","Nc_2-2","g(r)_2-3","Nc_2-3","g(r)_2-4","Nc_2-4","g(r)_3-3","Nc_3-3", &
  "g(r)_3-4","Nc_3-4","g(r)_4-4","Nc_4-4"
  y = 0.0
  y_rdf = 0.0
  DO j = 1,number_of_bins_gr
    rdf_ave(j) = rdf_ave(j)/nframes
    rdf_ave_per_bead(j,:) = rdf_ave_per_bead(j,:)/nframes
    y(j) =  rdf_ave(j) * distance_rdf(j) * distance_rdf(j)
    y_rdf(j,:) = rdf_ave_per_bead(j,:) *  distance_rdf(j) * distance_rdf(j)
    WRITE(2,*) distance_rdf(j) , rdf_ave(j), 4.0 * pi * rho * trapezios(y,distance_rdf,number_of_bins_gr), & 
               (rdf_ave_per_bead(j,z), 4.0 * pi * rho_rdf(z) * trapezios(y_rdf(:,z),distance_rdf,number_of_bins_gr), &  
               z=1,10)
  END DO
  
  pico1 = minloc(rdf_ave,DIM=1) !posição do primeiro mínimo 

  DEALLOCATE(y,rdf_types,y_rdf)
  CLOSE(2)
  WRITE(*,*); WRITE(*,*) 
END SUBROUTINE radial_distri 


SUBROUTINE alloc
    
  IMPLICIT NONE
  
  ALLOCATE(rdf_ave(number_of_bins_gr))
  ALLOCATE(rxyz(natoms,3,nframes),vxyz_frames(natoms,3,nframes),Nc(natoms,nframes),index_part(natoms))
  ALLOCATE(S_qx(number_of_bins_sq),S_qy(number_of_bins_sq),S_qz(number_of_bins_sq))
  ALLOCATE(S_q(number_of_bins_sq) , bond_length(nframes))
  ALLOCATE(q(number_of_bins_sq),rdf(number_of_bins_gr,nframes),distance_rdf(number_of_bins_gr))
  ALLOCATE(Fsx(nframes),Fsy(nframes),Fsz(nframes),MSD(nframes),Dt(nframes),Pc(13,nframes),group_size(13,nframes))
  ALLOCATE(S_qx_ave(number_of_bins_sq),S_qy_ave(number_of_bins_sq),S_qz_ave(number_of_bins_sq),S_q_ave(number_of_bins_sq))
  ALLOCATE(Ls(nframes),rdf_smooth(number_of_bins_gr - num_poin_appro_rdf + 1,nframes))
  ALLOCATE(rdf_fitted(n_interval + 1,nframes),distance_rdf_fitted(n_interval + 1)) !It dont take the final 
  ALLOCATE(rx(n_molecules,n_parti_per_molecule,nframes))    
  ALLOCATE(ry(n_molecules,n_parti_per_molecule,nframes)) 
  ALLOCATE(rz(n_molecules,n_parti_per_molecule,nframes)) 
  ALLOCATE(rx2(n_molecules_orien,n_parti_per_molecule_orien,nframes))    
  ALLOCATE(ry2(n_molecules_orien,n_parti_per_molecule_orien,nframes)) 
  ALLOCATE(rz2(n_molecules_orien,n_parti_per_molecule_orien,nframes))  
  !ALLOCATE(index(300000),particle_type(natoms,nframes))  para  quando os indices das moleculas forem maiores que o numero total de moleculas da caixa (mais de um componente na simulação) 
  ALLOCATE(index(n_molecules),particle_type(natoms,nframes)) 
  ALLOCATE(n_atoms_rdf(4),rho_rdf(10))
  n_atoms_rdf = 0
  rho_rdf = 0.0
  particle_type = 0
  rx = 0.0
  ry = 0.0
  rz = 0.0
  rx2 = 0.0
  ry2 = 0.0
  rz2 = 0.0  
  index = 0
  distance_rdf_fitted = 0.0
  rdf_fitted = 0.0
  Ls = 0.0
  bond_length = 0.0
  S_qx_ave = 0.0 
  S_qy_ave = 0.0
  S_qz_ave = 0.0 
  S_q_ave  = 0.0
  group_size = 0.0
  Pc = 0.0
  Dt = 0.0
  MSD = 0.0
  Fsx = 0.0
  Fsx = 0.0
  Fsx = 0.0
  q = 0.0
  S_qx = 0.0
  S_qy = 0.0
  S_qz = 0.0 
  rdf = 0.0
  rdf_smooth = 0.0
  distance_rdf = 0.0
  Nc = 0.0

END SUBROUTINE  alloc  

SUBROUTINE dealloc
    
  IMPLICIT NONE
  
  DEALLOCATE(S_qx , S_qy , S_qz , q , distance_rdf , rdf,S_qx_ave , S_qy_ave , S_qz_ave , S_q_ave)
  DEALLOCATE(S_q,Nc,rxyz,Fsx,Fsy,Fsz,MSD,Dt,Pc,group_size,index_part,bond_length,Ls,rx,ry,rz,index,particle_type)
  DEALLOCATE(n_atoms_rdf,rho_rdf,rx2,ry2,rz2)

END SUBROUTINE  dealloc


SUBROUTINE ordenar_cres(x,pos) !ordena números positivos não repetidos em ordem crescente, quero posição maximo
    REAL(8), DIMENSION(:), INTENT(INOUT)  :: x 
    integer, INTENT(INOUT)  :: pos 
    INTEGER :: i,j,N
    REAL ::  aux,aux2 
    N=size(x)
    
    aux = 0.0
    aux2 = 0.0
    
    do i = N,1, -1
        do j = N,1,-1
            if (x(i) > x(j)) then
              aux = x(i)
              x(i) = x(j)
              x(j) = aux
              IF (aux > aux2) THEN
                pos = i
                aux2 = aux
              END IF
            end if
         end do
    end do
    
END SUBROUTINE ordenar_cres

SUBROUTINE ordenar_decres(x,pos) !ordena números positivos não repetidos em ordem decrescente , quero posição mínimo
    REAL(8), DIMENSION(:), INTENT(INOUT)  :: x 
    INTEGER, INTENT(INOUT)  :: pos 
    INTEGER :: i,j,N
    REAL ::  aux,aux2 
    N=size(x)
    
    aux = 10000000
    aux2 = 10000000
    
    DO i = 1,N-1!N,1, -1
        DO j = 1,N!N,1,-1
            IF (x(i) < x(j)) then
              aux = x(i)
              x(i) = x(j)
              x(j) = aux
              IF (aux < aux2) THEN
                pos = i
                aux2 = aux
              END IF
            END IF
        END DO
    END DO
    
END SUBROUTINE ordenar_decres

SUBROUTINE read_data
  
  READ(*,*);  READ(*,*);  READ(*,*) input_file , delta_t , dump_interval , radius , cutoff , Kb , Temp
  READ(*,*);  READ(*,*);  READ(*,*); READ(*,*); READ(*,*); READ(*,*) Pc_flag
  READ(*,*);  READ(*,*);  READ(*,*) radial_distribution_flag , number_of_bins_gr , gr_Tgap 
  READ(*,*);  READ(*,*);  READ(*,*) density_profile_flag , number_of_bins_dens_profile , dens_pro_Tgap, axis_dens
  READ(*,*);  READ(*,*);  READ(*,*) orien_profile_flag, number_of_bins_orien_profile, orien_pro_Tgap, axis_orien, &
  n_molecules_orien, n_parti_per_molecule_orien, n_components
  READ(*,*);  READ(*,*);  READ(*,*) static_factor_flag , static_factor_option , number_of_bins_sq 
  READ(*,*);  READ(*,*);  READ(*,*); READ(*,*); READ(*,*); READ(*,*) Fs_flag , Fs_Tgap , Fs_q
  READ(*,*);  READ(*,*);  READ(*,*) Fc_flag , Fc_Tgap , Fc_q
  READ(*,*);  READ(*,*);  READ(*,*) MSD_flag , MSD_Tgap ,MSD_by_contact_number_flag
  READ(*,*);  READ(*,*);  READ(*,*) Nc_by_contact_number_flag 
  READ(*,*);  READ(*,*);  READ(*,*) cluster_flag, dist_cluster , n_molecules , n_parti_per_molecule, max_contacts
  READ(*,*);  READ(*,*) mm(1),mm(2),mm(3),mm(4)
  READ(*,*);  READ(*,*) total_mass, n_bins_angle 
  READ(*,*);  READ(*,*);  READ(*,*); READ(*,*); READ(*,*); READ(*,*) flag_rheology,stress_file, delta_t_stress , &
  dump_interval_stress, tp
  READ(*,*);  READ(*,*);  READ(*,*) stress_acf_flag, stress_acf_Tgap , stress_acf_option, correl_leng
  READ(*,*);  READ(*,*);  READ(*,*) G_primes_flag,  G_primes_option, strain_ampli
  
END SUBROUTINE read_data  


SUBROUTINE get_particle_positions
  INTEGER :: control = 0 , i , j , k , atom_id , atom_type , mol_id , c
  REAL(8) :: x , y , z , v_x, v_y, v_z , a1, a2, b1, b2, c1, c2
  INTEGER, DIMENSION(:), ALLOCATABLE :: aux_vec
  INTEGER, DIMENSION(n_molecules_orien) :: aux_n_par_per_mol
  INTEGER, DIMENSION(n_molecules_orien) :: aux_n_mol
  LOGICAL :: verify 
  
  aux_n_par_per_mol = 0
  aux_n_mol = 0
  nframes = 0   
  
  OPEN(unit=1,file=input_file,status="old",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening input_file from lammps"
  
  DO WHILE (control == 0)
         nframes = nframes + 1 !getting number of lines
         READ (1 , * , iostat= control) 
  END DO

  REWIND 1
  
  READ(1,*)
  READ(1,*)
  READ(1,*)
  READ(1,*) natoms 
  READ(1,*)
  READ(1,*) a1, a2
  READ(1,*) b1, b2
  READ(1,*) c1, c2
  READ(1,*)
  
  Lx = abs(2.0 * a2)
  Ly = abs(2.0 * b2)
  Lz = abs(2.0 * c2) 
  box_lower_limit = a1
  box_upper_limit = a2
  nframes = nframes / (natoms + 9)
  L = box_upper_limit - box_lower_limit
  rho = natoms / (L * L * L) 
  
  CALL alloc
  ALLOCATE(aux_vec(natoms)) 
  aux_vec = 0  
  
  !getting all atoms positions
  DO i = 1,nframes 
    index = 0    
    aux_n_par_per_mol = 0
    aux_n_mol = 0  
    c = 0

    DO j = 1,natoms    
        !READ(1,*) atom_id, atom_type , x , y , z !, v_x, v_y, v_z        
        READ(1,*) atom_id , mol_id , atom_type , x , y , z !, v_x, v_y, v_z

        ! Positions to calculate rdf, Sq, Fs, Fc, MSD
        rxyz(j,1,i) = x
        rxyz(j,2,i) = y
        rxyz(j,3,i) = z
        particle_type(j,i) = atom_type ! particle_type is a function of the frame        
        !vxyz_frames(j,1,i) = v_x
        !vxyz_frames(j,2,i) = v_y
        !vxyz_frames(j,3,i) = v_z
        
        IF (cluster_flag == 1) THEN
           ! write(*,*) i,j,x,y,z, mol_id
          ! Positions to calculate cluster properties, angle probabilities
          index(mol_id) = index(mol_id) + 1
          rx(mol_id,index(mol_id),i) = x 
          ry(mol_id,index(mol_id),i) = y
          rz(mol_id,index(mol_id),i) = z
        END IF
        
        IF (orien_profile_flag == 1) THEN
          verify = .false.  
          kloop1: DO k = 1,n_molecules_orien
                    IF (aux_n_mol(k) == mol_id) THEN
                      verify = .true.  
                      c = k
                      EXIT kloop1
                    END IF    
                  END DO kloop1  
          
          IF (.not. verify) THEN
            kloop2: DO k = 1,n_molecules_orien
                      IF (aux_n_mol(k) == 0) THEN
                        c = k  
                        aux_n_mol(c) = mol_id 
                        EXIT kloop2
                      END IF
                    END DO kloop2                           
          END IF    
          
    
          aux_n_par_per_mol(c) = aux_n_par_per_mol(c) + 1
          rx2(c,aux_n_par_per_mol(c),i) = x 
          ry2(c,aux_n_par_per_mol(c),i) = y
          rz2(c,aux_n_par_per_mol(c),i) = z
          
                !write(*,*) atom_id, mol_id,  c, aux_n_mol(c), aux_n_par_per_mol(c) 
        END IF        
        
    END DO      
  
    IF (i < nframes) THEN 
      DO j = 1,9
        READ(1,*)
      END DO
    END IF
  
  END DO     
  !end getting all atoms positions
  
  
  
  DO j = 1, natoms
    aux_vec(j) = particle_type(j,1) ! particle types of first frame  
  END DO    
  n_types = MAXVAL(aux_vec)
 
  !Number of particles of each type in the first frame, it is constant in time
  DO j = 1,natoms
    SELECT CASE (aux_vec(j))
      CASE (1)
        n_atoms_rdf(1) = n_atoms_rdf(1) + 1
      CASE (2)
        n_atoms_rdf(2) = n_atoms_rdf(2) + 1
      CASE (3)
        n_atoms_rdf(3) = n_atoms_rdf(3) + 1 
      CASE (4)
        n_atoms_rdf(4) = n_atoms_rdf(4) + 1      
    END SELECT        
  END DO      

  ! Calculating rho for each bead-bead interaction
  DO i = 1,nframes
    rho_rdf(1) = REAL(n_atoms_rdf(1),8)*(1.0/L/L/L) 
    rho_rdf(2) = REAL(n_atoms_rdf(1),8)*(1.0/L/L/L)  + REAL(n_atoms_rdf(2),8)*(1.0/L/L/L)  
    rho_rdf(3) = REAL(n_atoms_rdf(1),8)*(1.0/L/L/L)  + REAL(n_atoms_rdf(3),8)*(1.0/L/L/L) 
    rho_rdf(4) = REAL(n_atoms_rdf(1),8)*(1.0/L/L/L)  + REAL(n_atoms_rdf(4),8)*(1.0/L/L/L) 
    rho_rdf(5) = REAL(n_atoms_rdf(2),8)*(1.0/L/L/L) 
    rho_rdf(6) = REAL(n_atoms_rdf(2),8)*(1.0/L/L/L)  + REAL(n_atoms_rdf(3),8)*(1.0/L/L/L)  
    rho_rdf(7) = REAL(n_atoms_rdf(2),8)*(1.0/L/L/L)  + REAL(n_atoms_rdf(4),8)*(1.0/L/L/L) 
    rho_rdf(8) = REAL(n_atoms_rdf(3),8)*(1.0/L/L/L) 
    rho_rdf(9) = REAL(n_atoms_rdf(3),8)*(1.0/L/L/L)  + REAL(n_atoms_rdf(4),8)*(1.0/L/L/L)  
    rho_rdf(10) = REAL(n_atoms_rdf(4),8)*(1.0/L/L/L) 
  END DO
  
  CLOSE(1)
  DEALLOCATE(aux_vec)
END SUBROUTINE get_particle_positions

SUBROUTINE density_profile
INTEGER :: i, j, k, status_open, T, bin
REAL(8) :: delta_dens_pro, position
REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: dens_pro ! 3 Dimensions: Nº comp, Nº bins, Nº frames
REAL(8), DIMENSION(:), ALLOCATABLE :: discrete_space
REAL(8), DIMENSION(:,:), ALLOCATABLE :: average

  density_profile_flag = 1

  ALLOCATE(dens_pro(number_of_bins_dens_profile-1,n_components,nframes))
  ALLOCATE(discrete_space(number_of_bins_dens_profile-1))
  ALLOCATE(average(number_of_bins_dens_profile-1,n_components))
  
    
  OPEN(unit=1,file="density_profile.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error creating density_profile.out"
  
  dens_pro = 0.
  average = 0.
  SELECT CASE (axis_dens)
    CASE (1)
    delta_dens_pro = Lx / REAL(number_of_bins_dens_profile-1,8)
    
    DO i = 1, number_of_bins_dens_profile-1
      discrete_space(i) = (i-1)*delta_dens_pro+delta_dens_pro
    END DO
    
    WRITE (*,*) "Calculating Linear Density Profile"
    DO i = 1, nframes
      WRITE(*,*) "Frame = ", i  
      DO j = 1, natoms 
          position = rxyz(j,1,i)
          
          IF  (position > 0) THEN
            T = floor(position / (Lx / 2.))

            IF (mod(T,2) == 1) THEN
	      position = position - (T + 1) * (Lx / 2.)
	    ELSE
	      position = position - T * (Lx / 2.)
	    END IF

	  ELSE
	    T = ceiling(position / (Lx / 2.))

	    IF (mod(T,2) == -1) THEN
	      position = position - (T - 1) * (Lx / 2.)
	    ELSE
	      position = position - T * (Lx / 2.)
	    END IF

	  END IF

	  IF  (position > 0) THEN
            bin = floor(number_of_bins_dens_profile / 2.) + floor(position / delta_dens_pro)
	  ELSE
	    bin = floor(number_of_bins_dens_profile / 2.) + ceiling(position / delta_dens_pro)
	  END IF
          
          IF (n_components == 1) THEN
	    dens_pro(bin,1,i) = dens_pro(bin,1,i)+1
          ELSE
            dens_pro(bin,particle_type(j,i),i) = dens_pro(bin,particle_type(j,i),i)+1  
          END IF  
          
      END DO
    END DO
    
    
    !Writing output file  of linear density  
    DO i = 1, nframes, dens_pro_Tgap
      WRITE(1,*) "Frame = ", i  
      WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)"
      DO j = 1, number_of_bins_dens_profile - 1
        WRITE(1,*)  discrete_space(j), dens_pro(j,:,i) / 1. / 1. / delta_dens_pro ! dens_pro(j,:,i) / Ly / Lz / DELTAX
      END DO    
      WRITE(1,*)
    END DO
    
    !Calculating and writing averaged values
    DO i = 1, nframes
      DO j = 1, number_of_bins_dens_profile - 1 
        average(j,:) = average(j,:) + dens_pro(j,:,i)   
      END DO
    END DO    
    
    DO j = 1, number_of_bins_dens_profile - 1 
      average(j,:) = average(j,:) / nframes
    END DO 
    
    WRITE(1,*) "Averaged values"
    WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)"
    DO j = 1, number_of_bins_dens_profile - 1
      WRITE(1,*) discrete_space(j), average(j,:) / 1. / 1. / delta_dens_pro ! dens_pro(j,:,i) / Ly (1.0) / Lz (1.0) / DELTAX
    END DO      
    
    CASE (2)
    delta_dens_pro = Ly / REAL(number_of_bins_dens_profile-1,8)
    
    DO i = 1, number_of_bins_dens_profile-1
      discrete_space(i) = (i-1)*delta_dens_pro+delta_dens_pro
    END DO
    
    WRITE (*,*) "Calculating Linear Density Profile"
    DO i = 1, nframes
      WRITE(*,*) "Frame = ", i  
      DO j = 1, natoms 
          position = rxyz(j,2,i)
          
          IF  (position > 0) THEN
            T = floor(position / (Ly / 2.))

            IF (mod(T,2) == 1) THEN
	      position = position - (T + 1) * (Ly / 2.)
	    ELSE
	      position = position - T * (Ly / 2.)
	    END IF

	  ELSE
	    T = ceiling(position / (Ly / 2.))

	    IF (mod(T,2) == -1) THEN
	      position = position - (T - 1) * (Ly / 2.)
	    ELSE
	      position = position - T * (Ly / 2.)
	    END IF

	  END IF

	  IF  (position > 0) THEN
            bin = floor(number_of_bins_dens_profile / 2.) + floor(position / delta_dens_pro)
	  ELSE
	    bin = floor(number_of_bins_dens_profile / 2.) + ceiling(position / delta_dens_pro)
	  END IF

          IF (n_components == 1) THEN
	    dens_pro(bin,1,i) = dens_pro(bin,1,i)+1
          ELSE
            dens_pro(bin,particle_type(j,i),i) = dens_pro(bin,particle_type(j,i),i)+1  
          END IF 
          
      END DO
    END DO
    
    
    !Writing output file  of linear density  
    DO i = 1, nframes, dens_pro_Tgap
      WRITE(1,*) "Frame = ", i  
      WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)"
      DO j = 1, number_of_bins_dens_profile - 1
        WRITE(1,*)  discrete_space(j), dens_pro(j,:,i) / 1. / 1. / delta_dens_pro ! dens_pro(j,:,i) / Lx / Lz / DELTAY
      END DO    
      WRITE(1,*)
    END DO
    
    !Calculating and writing averaged values
    DO i = 1, nframes
      DO j = 1, number_of_bins_dens_profile - 1 
        average(j,:) = average(j,:) + dens_pro(j,:,i)   
      END DO
    END DO    
    
    DO j = 1, number_of_bins_dens_profile - 1 
      average(j,:) = average(j,:) / nframes
    END DO 
    
    WRITE(1,*) "Averaged values"
    WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)"
    DO j = 1, number_of_bins_dens_profile - 1
      WRITE(1,*) discrete_space(j), average(j,:) / 1. / 1. / delta_dens_pro ! dens_pro(j,:,i) / Lx (1.0) / Lz (1.0) / DELTAY
    END DO      
    
    CASE (3)
    delta_dens_pro = Lz / REAL(number_of_bins_dens_profile-1,8)
    
    DO i = 1, number_of_bins_dens_profile-1
      discrete_space(i) = (i-1)*delta_dens_pro+delta_dens_pro
    END DO
    
    WRITE (*,*) "Calculating Linear Density Profile"
    DO i = 1, nframes
      WRITE(*,*) "Frame = ", i  
      DO j = 1, natoms 
          
         ! IF (particle_type(j,i) == 2) THEN
            position = rxyz(j,3,i)
         ! ELSE 
         !     CYCLE
         ! END IF
          
          
          IF  (position > 0) THEN
            T = floor(position / (Lz / 2.))

            IF (mod(T,2) == 1) THEN
	      position = position - (T + 1) * (Lz / 2.)
	    ELSE
	      position = position - T * (Lz / 2.)
	    END IF

	  ELSE
	    T = ceiling(position / (Lz / 2.))

	    IF (mod(T,2) == -1) THEN
	      position = position - (T - 1) * (Lz / 2.)
	    ELSE
	      position = position - T * (Lz / 2.)
	    END IF

	  END IF

	  IF  (position > 0) THEN
            bin = floor(number_of_bins_dens_profile / 2.) + floor(position / delta_dens_pro)
	  ELSE
	    bin = floor(number_of_bins_dens_profile / 2.) + ceiling(position / delta_dens_pro)
	  END IF

          IF (n_components == 1) THEN
	    dens_pro(bin,1,i) = dens_pro(bin,1,i)+1
          ELSE
            dens_pro(bin,particle_type(j,i),i) = dens_pro(bin,particle_type(j,i),i)+1  
          END IF 
          
      END DO
    END DO
    
    
    !Writing output file  of linear density  
    DO i = 1, nframes, dens_pro_Tgap
      WRITE(1,*) "Frame = ", i  
      WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)"
      DO j = 1, number_of_bins_dens_profile - 1
        WRITE(1,*)  discrete_space(j), dens_pro(j,:,i) / 1. / 1. / delta_dens_pro ! dens_pro(j,:,i) / Lx / Ly / DELTAZ
      END DO    
      WRITE(1,*)
    END DO
    
    !Calculating and writing averaged values
    DO i = 1, nframes
      DO j = 1, number_of_bins_dens_profile - 1 
        average(j,:) = average(j,:) + dens_pro(j,:,i)   
      END DO
    END DO    
    
    DO j = 1, number_of_bins_dens_profile - 1 
      average(j,:) = average(j,:) / nframes
    END DO 
    
    WRITE(1,*) "Averaged values"
    WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)"
    DO j = 1, number_of_bins_dens_profile - 1
      WRITE(1,*) discrete_space(j), average(j,:) / 1. / 1. / delta_dens_pro ! dens_pro(j,:,i) / Lx (1.0) / Ly (1.0) / DELTAZ
    END DO      
      
  END SELECT     
  
  DEALLOCATE(dens_pro,discrete_space,average)
  CLOSE(1)
END SUBROUTINE density_profile    


SUBROUTINE orientation_profile
INTEGER :: i, j, k, status_open, T, bin
REAL(8) :: delta_orien_pro, position, angle
REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: orien_pro, dens_pro ! 3 Dimensions: Nº comp, Nº bins, Nº frames
REAL(8), DIMENSION(:), ALLOCATABLE :: discrete_space
REAL(8), DIMENSION(:,:), ALLOCATABLE :: average_dens, average_orien
REAL(8), DIMENSION(n_parti_per_molecule_orien) :: x, y, z

  orien_profile_flag = 1
  ALLOCATE(orien_pro(number_of_bins_orien_profile-1,n_components,nframes))
  ALLOCATE(dens_pro(number_of_bins_orien_profile-1,n_components,nframes))  
  ALLOCATE(discrete_space(number_of_bins_orien_profile-1))
  ALLOCATE(average_orien(number_of_bins_orien_profile-1,n_components))
  ALLOCATE(average_dens(number_of_bins_orien_profile-1,n_components))  
    
  OPEN(unit=1,file="orientation_profile.out", status="replace", iostat = status_open ) 
  IF (status_open > 0) STOP "error creating orientation_profile.out"
  
  orien_pro = 0.
  dens_pro = 0.  
  average_dens = 0.
  average_orien = 0.  
  SELECT CASE (axis_orien)
    CASE (1)
    delta_orien_pro = Lx / REAL(number_of_bins_orien_profile-1,8)
    
    DO i = 1, number_of_bins_orien_profile-1
      discrete_space(i) = (i-1)*delta_orien_pro+delta_orien_pro
    END DO
    
    WRITE (*,*) "Calculating Orientation Profile"
    DO i = 1, nframes
      WRITE(*,*) "Frame = ", i  
      DO j = 1, n_molecules_orien 
          x = 0.
          y = 0.
          z = 0.
          position = 0.
          angle = 0.  
          
          x = rx2(j,:,i)
          y = ry2(j,:,i)
          z = rz2(j,:,i)
          
          CALL CoM_angle_orien_pro(x,y,z,position,angle)
          
          !POSITION is the x coodinate of Center of Mass Molecule, all particles should have the same mass
          !ANGLE it the Angle between unit vector from first to last particle in the molecule and interface
          
          IF  (position > 0) THEN
            T = floor(position / (Lx / 2.))

            IF (mod(T,2) == 1) THEN
	      position = position - (T + 1) * (Lx / 2.)
	    ELSE
	      position = position - T * (Lx / 2.)
	    END IF

	  ELSE
	    T = ceiling(position / (Lx / 2.))

	    IF (mod(T,2) == -1) THEN
	      position = position - (T - 1) * (Lx / 2.)
	    ELSE
	      position = position - T * (Lx / 2.)
	    END IF

	  END IF

	  IF  (position > 0) THEN
            bin = floor(number_of_bins_orien_profile / 2.) + floor(position / delta_orien_pro)
	  ELSE
	    bin = floor(number_of_bins_orien_profile / 2.) + ceiling(position / delta_orien_pro)
	  END IF

          IF (n_components == 1) THEN
	    orien_pro(bin,1,i) = orien_pro(bin,1,i)+angle
	    dens_pro(bin,1,i) = dens_pro(bin,1,i)+1          
          ELSE
	    orien_pro(bin,particle_type(j,i),i) = orien_pro(bin,particle_type(j,i),i)+angle
	    dens_pro(bin,particle_type(j,i),i) = dens_pro(bin,particle_type(j,i),i)+1               
          END IF    
      END DO
    END DO
    
    
    !Writing output file  of orientation profile 
    DO i = 1, nframes, orien_pro_Tgap
      WRITE(1,*) "Frame = ", i  
      WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)", &
      "     Orientation Profile of components (1 to Max Nº of components)"
      DO j = 1, number_of_bins_orien_profile - 1 
        IF (n_components == 1) THEN  
          IF (dens_pro(j,1,i) /= 0.) THEN  
            WRITE(1,*)  discrete_space(j), dens_pro(j,1,i) / 1. / 1. / delta_orien_pro, orien_pro(j,1,i) / dens_pro(j,1,i) ! dens_pro(j,:,i) / Ly / Lz / DELTAX
          ELSE
            WRITE(1,*)  discrete_space(j), dens_pro(j,1,i) / 1. / 1. / delta_orien_pro, 0. ! dens_pro(j,:,i) / Ly / Lz / DELTAX
          END IF
        ELSE
          STOP "por enquanto só serve pra um componente"   
        END IF
      END DO    
      WRITE(1,*)
    END DO
    
    !Calculating and writing averaged values
    DO i = 1, nframes
      DO j = 1, number_of_bins_orien_profile - 1 
        average_orien(j,:) = average_orien(j,:) + orien_pro(j,:,i)   
        average_dens(j,:) = average_dens(j,:) + dens_pro(j,:,i) 
        
      END DO
    END DO    

    DO j = 1, number_of_bins_orien_profile - 1 
      average_orien(j,:) = average_orien(j,:) / nframes
      average_dens(j,:) = average_dens(j,:) / nframes  
    END DO 
    
    WRITE(1,*) "Averaged values"
    WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)", &
    "     Orientation Profile of components (1 to Max Nº of components)"    
    DO j = 1, number_of_bins_orien_profile - 1
        IF (n_components == 1) THEN  
          IF (average_dens(j,1) /= 0.) THEN  
            WRITE(1,*)  discrete_space(j), average_dens(j,1) / 1. / 1. / delta_orien_pro, average_orien(j,1) / average_dens(j,1) ! dens_pro(j,:,i) / Ly / Lz / DELTAX
          ELSE
            WRITE(1,*)  discrete_space(j), average_dens(j,1) / 1. / 1. / delta_orien_pro, 0. ! dens_pro(j,:,i) / Ly / Lz / DELTAX
          END IF
        ELSE
          STOP "por enquanto só serve pra um componente"   
        END IF               
    END DO      
    
    CASE (2)
    delta_orien_pro = Ly / REAL(number_of_bins_orien_profile-1,8)
    
    DO i = 1, number_of_bins_orien_profile-1
      discrete_space(i) = (i-1)*delta_orien_pro+delta_orien_pro
    END DO
    
    WRITE (*,*) "Calculating Orientation Profile"
    DO i = 1, nframes
      WRITE(*,*) "Frame = ", i  
      DO j = 1, n_molecules_orien 
          x = 0.
          y = 0.
          z = 0.
          position = 0.
          angle = 0.  
          
          x = rx2(j,:,i)
          y = ry2(j,:,i)
          z = rz2(j,:,i)
          
          CALL CoM_angle_orien_pro(x,y,z,position,angle)
          
          !POSITION is the y coodinate of Center of Mass Molecule, all particles should have the same mass
          !ANGLE it the Angle between unit vector from first to last particle in the molecule and interface
          
          IF  (position > 0) THEN
            T = floor(position / (Ly / 2.))

            IF (mod(T,2) == 1) THEN
	      position = position - (T + 1) * (Ly / 2.)
	    ELSE
	      position = position - T * (Ly / 2.)
	    END IF

	  ELSE
	    T = ceiling(position / (Ly / 2.))

	    IF (mod(T,2) == -1) THEN
	      position = position - (T - 1) * (Ly / 2.)
	    ELSE
	      position = position - T * (Ly / 2.)
	    END IF

	  END IF

	  IF  (position > 0) THEN
            bin = floor(number_of_bins_orien_profile / 2.) + floor(position / delta_orien_pro)
	  ELSE
	    bin = floor(number_of_bins_orien_profile / 2.) + ceiling(position / delta_orien_pro)
	  END IF

          IF (n_components == 1) THEN
	    orien_pro(bin,1,i) = orien_pro(bin,1,i)+angle
	    dens_pro(bin,1,i) = dens_pro(bin,1,i)+1          
          ELSE
	    orien_pro(bin,particle_type(j,i),i) = orien_pro(bin,particle_type(j,i),i)+angle
	    dens_pro(bin,particle_type(j,i),i) = dens_pro(bin,particle_type(j,i),i)+1               
          END IF    
      END DO
    END DO
    
    
    !Writing output file  of orientation profile 
    DO i = 1, nframes, orien_pro_Tgap
      WRITE(1,*) "Frame = ", i  
      WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)", &
      "     Orientation Profile of components (1 to Max Nº of components)"
      DO j = 1, number_of_bins_orien_profile - 1 
        IF (n_components == 1) THEN  
          IF (dens_pro(j,1,i) /= 0.) THEN  
            WRITE(1,*)  discrete_space(j), dens_pro(j,1,i) / 1. / 1. / delta_orien_pro, orien_pro(j,1,i) / dens_pro(j,1,i) ! dens_pro(j,:,i) / Ly / Lz / DELTAY
          ELSE
            WRITE(1,*)  discrete_space(j), dens_pro(j,1,i) / 1. / 1. / delta_orien_pro, 0. ! dens_pro(j,:,i) / Ly / Lz / DELTAY
          END IF
        ELSE
          STOP "por enquanto só serve pra um componente"   
        END IF
      END DO    
      WRITE(1,*)
    END DO
    
    !Calculating and writing averaged values
    DO i = 1, nframes
      DO j = 1, number_of_bins_orien_profile - 1 
        average_orien(j,:) = average_orien(j,:) + orien_pro(j,:,i)   
        average_dens(j,:) = average_dens(j,:) + dens_pro(j,:,i) 
        
      END DO
    END DO    

    DO j = 1, number_of_bins_orien_profile - 1 
      average_orien(j,:) = average_orien(j,:) / nframes
      average_dens(j,:) = average_dens(j,:) / nframes  
    END DO 
    
    WRITE(1,*) "Averaged values"
    WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)", &
    "     Orientation Profile of components (1 to Max Nº of components)"    
    DO j = 1, number_of_bins_orien_profile - 1
        IF (n_components == 1) THEN  
          IF (average_dens(j,1) /= 0.) THEN  
            WRITE(1,*)  discrete_space(j), average_dens(j,1) / 1. / 1. / delta_orien_pro, average_orien(j,1) / average_dens(j,1) ! dens_pro(j,:,i) / Ly / Lz / DELTAY
          ELSE
            WRITE(1,*)  discrete_space(j), average_dens(j,1) / 1. / 1. / delta_orien_pro, 0. ! dens_pro(j,:,i) / Ly / Lz / DELTAY
          END IF
        ELSE
          STOP "por enquanto só serve pra um componente"   
        END IF               
    END DO  
        
        
    CASE (3)
    delta_orien_pro = Lz / REAL(number_of_bins_orien_profile-1,8)
    
    DO i = 1, number_of_bins_orien_profile-1
      discrete_space(i) = (i-1)*delta_orien_pro+delta_orien_pro
    END DO
    
    WRITE (*,*) "Calculating Orientation Profile"
    DO i = 1, nframes
      WRITE(*,*) "Frame = ", i  
      DO j = 1, n_molecules_orien 
          x = 0.
          y = 0.
          z = 0.
          position = 0.
          angle = 0.  
          
          x = rx2(j,:,i)
          y = ry2(j,:,i)
          z = rz2(j,:,i)
          
          CALL CoM_angle_orien_pro(x,y,z,position,angle)
          
          !POSITION is the z coodinate of Center of Mass Molecule, all particles should have the same mass
          !ANGLE it the Angle between unit vector from first to last particle in the molecule and interface
          
          IF  (position > 0) THEN
            T = floor(position / (Lz / 2.))

            IF (mod(T,2) == 1) THEN
	      position = position - (T + 1) * (Lz / 2.)
	    ELSE
	      position = position - T * (Lz / 2.)
	    END IF

	  ELSE
	    T = ceiling(position / (Lz / 2.))

	    IF (mod(T,2) == -1) THEN
	      position = position - (T - 1) * (Lz / 2.)
	    ELSE
	      position = position - T * (Lz / 2.)
	    END IF

	  END IF

	  IF  (position > 0) THEN
            bin = floor(number_of_bins_orien_profile / 2.) + floor(position / delta_orien_pro)
	  ELSE
	    bin = floor(number_of_bins_orien_profile / 2.) + ceiling(position / delta_orien_pro)
	  END IF

          IF (n_components == 1) THEN
	    orien_pro(bin,1,i) = orien_pro(bin,1,i)+angle
	    dens_pro(bin,1,i) = dens_pro(bin,1,i)+1          
          ELSE
	    orien_pro(bin,particle_type(j,i),i) = orien_pro(bin,particle_type(j,i),i)+angle
	    dens_pro(bin,particle_type(j,i),i) = dens_pro(bin,particle_type(j,i),i)+1               
          END IF    
      END DO
    END DO
    
    
    !Writing output file  of orientation profile 
    DO i = 1, nframes, orien_pro_Tgap
      WRITE(1,*) "Frame = ", i  
      WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)", &
      "     Orientation Profile of components (1 to Max Nº of components)"
      DO j = 1, number_of_bins_orien_profile - 1 
        IF (n_components == 1) THEN  
          IF (dens_pro(j,1,i) /= 0.) THEN  
            WRITE(1,*)  discrete_space(j), dens_pro(j,1,i) / 1. / 1. / delta_orien_pro, orien_pro(j,1,i) / dens_pro(j,1,i) ! dens_pro(j,:,i) / Ly / Lz / DELTAZ
          ELSE
            WRITE(1,*)  discrete_space(j), dens_pro(j,1,i) / 1. / 1. / delta_orien_pro, 0. ! dens_pro(j,:,i) / Ly / Lz / DELTAZ
          END IF
        ELSE
          STOP "por enquanto só serve pra um componente"   
        END IF
      END DO    
      WRITE(1,*)
    END DO
    
    !Calculating and writing averaged values
    DO i = 1, nframes
      DO j = 1, number_of_bins_orien_profile - 1 
        average_orien(j,:) = average_orien(j,:) + orien_pro(j,:,i)   
        average_dens(j,:) = average_dens(j,:) + dens_pro(j,:,i) 
        
      END DO
    END DO    

    DO j = 1, number_of_bins_orien_profile - 1 
      average_orien(j,:) = average_orien(j,:) / nframes
      average_dens(j,:) = average_dens(j,:) / nframes  
    END DO 
    
    WRITE(1,*) "Averaged values"
    WRITE(1,*) "r (Dist. units)", "     Linear Density of components (1 to Max Nº of components)", &
    "     Orientation Profile of components (1 to Max Nº of components)"    
    DO j = 1, number_of_bins_orien_profile - 1
        IF (n_components == 1) THEN  
          IF (average_dens(j,1) /= 0.) THEN  
            WRITE(1,*)  discrete_space(j), average_dens(j,1) / 1. / 1. / delta_orien_pro, average_orien(j,1) / average_dens(j,1) ! dens_pro(j,:,i) / Ly / Lz / DELTAZ
          ELSE
            WRITE(1,*)  discrete_space(j), average_dens(j,1) / 1. / 1. / delta_orien_pro, 0. ! dens_pro(j,:,i) / Ly / Lz / DELTAZ
          END IF
        ELSE
          STOP "por enquanto só serve pra um componente"   
        END IF               
    END DO         
      
  END SELECT     
  
  DEALLOCATE(orien_pro,dens_pro,discrete_space,average_dens,average_orien)
  CLOSE(1)
END SUBROUTINE orientation_profile   

SUBROUTINE  CoM_angle_orien_pro(x,y,z,position,angle)
IMPLICIT NONE    
REAL(8), DIMENSION(:), INTENT(INOUT) :: x,y,z
REAL(8), INTENT(OUT) :: position, angle
INTEGER :: i
REAL(8) :: sumx, sumy, sumz
REAL(8), DIMENSION(3) :: mole_direction, normal_to_interface !vector representing molecule direction

!All particles should have the same mass for CoM calculation

sumx = 0.
sumy = 0.
sumz = 0.

DO i = 1, size(x)
  sumx = sumx + x(i)
  sumy = sumy + y(i)  
  sumz = sumz + z(i)  
END DO    

SELECT CASE (axis_orien)
  CASE (1)
    position = sumx / size(x)  
  CASE(2)
    position = sumy / size(y) 
  CASE(3)
    position = sumz / size(z)        
END SELECT

! Calculating angle between vector from the first to the last particle in the molecule and interface 
SELECT CASE (axis_orien)
  CASE (1)
      mole_direction(1) = x(size(x)) - x(1)
      mole_direction(2) = y(size(y)) - y(1)
      mole_direction(3) = z(size(z)) - z(1)      
    !write(*,*) mole_direction(1),mole_direction(2),mole_direction(3)
      normal_to_interface(1) = 1.  
      normal_to_interface(2) = 0.
      normal_to_interface(3) = 0. 
      
      !Calculating cosine of the angle (alfa) between molecular direction and normal to interface
      angle = mole_direction(1) * normal_to_interface(1) + mole_direction(2) * normal_to_interface(2) + &
      mole_direction(3) * normal_to_interface(3)
      !write(*,*) angle
      angle = angle / sqrt(mole_direction(1)**2. + mole_direction(2)**2. + mole_direction(3)**2.)
      !write(*,*) angle
      angle = angle / sqrt(normal_to_interface(1)**2. + normal_to_interface(2)**2. + normal_to_interface(3)**2.)
      !write(*,*) angle
      !Calculating cosine of the angle (teta) between molecular direction and interface
      angle = 1 - (angle ** 2.)  !this is sin(alfa)² = 1 - cos(alfa)²
      angle = sqrt(angle)      !this is sin(alfa) = sqrt(sin(alfa)²)
      
      !alfa and teta are complementary so cos(teta) = sin(alfa)
      angle = angle            !this  is cos(teta) = sin(alfa)
      
      angle = acos(angle)     !this is teta = arcos(cos(teta)) this angle is in rad
      
      angle = angle * 180.0 / 3.14159265359 ! this angle is in degrees      
        !write(*,*) angle; write(*,*);write(*,*)
  CASE(2)
      mole_direction(1) = x(size(x)) - x(1)
      mole_direction(2) = y(size(y)) - y(1)
      mole_direction(3) = z(size(z)) - z(1)      
    
      normal_to_interface(1) = 0.  
      normal_to_interface(2) = 1.
      normal_to_interface(3) = 0. 
      
      !Calculating cosine of the angle (alfa) between molecular direction and normal to interface
      angle = mole_direction(1) * normal_to_interface(1) + mole_direction(2) * normal_to_interface(2) + &
      mole_direction(3) * normal_to_interface(3)
      angle = angle / sqrt(mole_direction(1)**2. + mole_direction(2)**2. + mole_direction(3)**2.)
      angle = angle / sqrt(normal_to_interface(1)**2. + normal_to_interface(2)**2. + normal_to_interface(3)**2.)
      
      !Calculating cosine of the angle (teta) between molecular direction and interface
      angle = 1 - (angle ** 2.)  !this is sin(alfa)² = 1 - cos(alfa)²
      angle = sqrt(angle)      !this is sin(alfa) = sqrt(sin(alfa)²)
      
      !alfa and teta are complementary so cos(teta) = sin(alfa)
      angle = angle            !this  is cos(teta) = sin(alfa)
      
      angle = acos(angle)     !this is teta = arcos(cos(teta)) this angle is in rad
      
      angle = angle * 180.0 / 3.14159265359 ! this angle is in degrees 
      
  CASE(3)
      mole_direction(1) = x(size(x)) - x(1)
      mole_direction(2) = y(size(y)) - y(1)
      mole_direction(3) = z(size(z)) - z(1)      
    
      normal_to_interface(1) = 0.  
      normal_to_interface(2) = 0.
      normal_to_interface(3) = 1. 
      
      !Calculating cosine of the angle (alfa) between molecular direction and normal to interface
      angle = mole_direction(1) * normal_to_interface(1) + mole_direction(2) * normal_to_interface(2) + &
      mole_direction(3) * normal_to_interface(3)
      angle = angle / sqrt(mole_direction(1)**2. + mole_direction(2)**2. + mole_direction(3)**2.)
      angle = angle / sqrt(normal_to_interface(1)**2. + normal_to_interface(2)**2. + normal_to_interface(3)**2.)
      
      !Calculating cosine of the angle (teta) between molecular direction and interface
      angle = 1 - (angle ** 2.)  !this is sin(alfa)² = 1 - cos(alfa)²
      angle = sqrt(angle)      !this is sin(alfa) = sqrt(sin(alfa)²)
      
      !alfa and teta are complementary so cos(teta) = sin(alfa)
      angle = angle            !this  is cos(teta) = sin(alfa)
      
      angle = acos(angle)     !this is teta = arcos(cos(teta)) this angle is in rad
      
      angle = angle * 180.0 / 3.14159265359 ! this angle is in degrees 
      
END SELECT  


END SUBROUTINE CoM_angle_orien_pro

END MODULE functions
