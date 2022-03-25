!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!
! Maximum number of clusters = n_molecules
!     
! File:   cluster_analysis.f90
! Author: fellipe
!
! Created on 21 de Fevereiro de 2019, 08:49
!

MODULE cluster_analysis
USE global_variables
USE fitting_planes
USE statistical_functions

CONTAINS

SUBROUTINE clustering
  IMPLICIT NONE    
  INTEGER :: k,i,j,m,n, count_links = 0, status_open, ii 
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: n_contacts_per_molecule
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: node_next
  INTEGER, DIMENSION(nframes) :: n_links
  REAL(8) :: dist, dx , dy , dz , dist_or, proba_incre, a, b, c, sum_sq_erro_fit_plane, alpha
  REAL(8), DIMENSION(:), ALLOCATABLE :: dados_molx,dados_moly,dados_molz  
  REAL(8), DIMENSION(:,:), ALLOCATABLE :: plane_equations 
  INTEGER, DIMENSION(:), ALLOCATABLE :: count_bin_angle
  INTEGER :: index_proba, counter_total, counter_d
  REAL(8) :: d, average , d_ave = 0.0 , aux2, aux3 
  REAL(8), DIMENSION(:), ALLOCATABLE :: proba_angle,proba_ave  
  
  !Maximun of max_contacts particles attached to one particle, so there are max_contacts maximum links
  ALLOCATE(node_next(n_molecules,max_contacts,nframes))
  ALLOCATE(n_contacts_per_molecule(n_molecules,nframes))
  ALLOCATE(mass(n_types))
  ALLOCATE(dados_molx(n_parti_per_molecule))
  ALLOCATE(dados_moly(n_parti_per_molecule))
  ALLOCATE(dados_molz(n_parti_per_molecule)) 
  ALLOCATE(plane_equations(n_molecules,4))  
  ALLOCATE(proba_angle(n_bins_angle),count_bin_angle(n_bins_angle),proba_ave(n_bins_angle))
  proba_ave = 0.0  
  
  DO i = 1,n_types
      mass(i) = mm(i)
  END DO 
  
  proba_incre = pi/2.0/REAL(n_bins_angle,8)
  
  OPEN(unit=5,file="proba_included_angle.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error writing proba_included_angle.out" 
  WRITE (5,*) "Incremento bin do angulo = ",  proba_incre
  WRITE (5,*);WRITE (5,*);  
  
  OPEN(unit=7,file="plane_equations.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error writing plane_equations.out" 
  
  !dist_or = dist_cluster
  !do ii = 1,200
  !dist_cluster = dist_or + real((ii-1),8)*0.2
      
  node_next = 0
  n_contacts_per_molecule = 0
  n_links = 0               ! Number of links per frame
  
  !Exibiting node_next
  OPEN(unit=2,file="neighbor_nodes.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening clustering.out"
  
  !beginning algorithm create nodes
  DO k = 1,nframes      
  WRITE(2,*) "Frame ",k  
  
  WRITE(7,*) "Frame = ",k
  WRITE(7,'(a9,5x,a2,24x,a2,24x,a1,26x,8a)') "Molecule ","a ","b ","c", "erro fit"
  ! Cálculo das equações dos planos de cada molécula no frame   
  
!  !Calculating plane equations
!  DO i = 1,n_molecules    
!    
!    DO j = 1, n_parti_per_molecule-1
!      dados_molx(j+1) = rx(i,j+1,k) - L*nint((rx(i,j+1,k)-rx(i,j,k))/L)
!      dados_moly(j+1) = ry(i,j+1,k) - L*nint((ry(i,j+1,k)-ry(i,j,k))/L)
!      dados_molz(j+1) = rz(i,j+1,k) - L*nint((rz(i,j+1,k)-rz(i,j,k))/L)
!    END DO
!    
!    ! A equação do plano é a*x+b*y+c*z+1 = 0
!    CALL fit_planes(dados_molx,dados_moly,dados_molz,a,b,c,sum_sq_erro_fit_plane) 
!    plane_equations(i,1) = a
!    plane_equations(i,2) = b
!    plane_equations(i,3) = c
!    plane_equations(i,4) = sum_sq_erro_fit_plane
!    
!    WRITE(7,*) i,a,b,c,sum_sq_erro_fit_plane
!      
!  END DO
!  WRITE(7,*); WRITE(7,*);   
!  
!  counter_total = 0
!  count_bin_angle = 0
!  counter_d = 0
!  d = 0.0    
  
  DO m = 1,n_molecules !loop for number of atoms
    DO n = 1,n_molecules        
      !Looking for some contact between two particles belonging to these two molecules
      !I'm doing it because molecules can stack perpendicular as well
          
      loop: DO i = 1, n_parti_per_molecule
              DO j = 1, n_parti_per_molecule
                dx = rx(m,i,k) - rx(n,j,k)
                dx = dx - L*nint(dx/L)
                dy = ry(m,i,k) - ry(n,j,k)
                dy = dy - L*nint(dy/L)
                dz = rz(m,i,k) - rz(n,j,k)  
                dz = dz - L*nint(dz/L)                
                dist =  dx*dx+dy*dy+dz*dz  
          
                IF ((dist < dist_cluster * dist_cluster) .AND. (m /= n)) THEN
!                  !Angle calculations
!                  alpha = plane_equations(m,1)*plane_equations(n,1)
!                  alpha = alpha + plane_equations(m,2)*plane_equations(n,2)
!                  alpha = alpha + plane_equations(m,3)*plane_equations(n,3)
!                  alpha = abs(alpha)
!           
!                  aux2 = plane_equations(m,1)*plane_equations(m,1)
!                  aux2 = aux2 + plane_equations(m,2)*plane_equations(m,2)
!                  aux2 = aux2 + plane_equations(m,3)*plane_equations(m,3)
!                  aux2 = sqrt(aux2)
!            
!                  aux3 = plane_equations(n,1)*plane_equations(n,1)
!                  aux3 = aux3 + plane_equations(n,2)*plane_equations(n,2)
!                  aux3 = aux3 + plane_equations(n,3)*plane_equations(n,3)
!                  aux3 = sqrt(aux3)  
!                    
!                  alpha = alpha / aux2 / aux3 ! cosseno do ângulo agudo entre planos
!                  IF (alpha > 1.0) alpha = 1.0
!                  alpha = acos(alpha)         ! ângulo agudo entre planos  em radianos 
!                  
!                  !Cálculo da distância média, apenas as moléculas paralelas entre si
!                  IF ((alpha) < proba_incre) THEN
!                    average = (plane_equations(m,1)/plane_equations(n,1))
!                    average = average + (plane_equations(m,2)/plane_equations(n,2))
!                    average = average + (plane_equations(m,3)/plane_equations(n,3))
!                    average = average / 3.0
!          
!                    d = d + abs(1.0*average-1.0)/aux2
!                    counter_d = counter_d + 1
!            
!                  END IF   
!
!                  DO index_proba = 1,n_bins_angle
!                    IF( (alpha >= (proba_incre*(index_proba-1)) ).AND.(alpha < (proba_incre*(index_proba)) ) )THEN
!                        count_bin_angle(index_proba) = count_bin_angle(index_proba) + 1      
!                    END IF                    
!                  END DO    
!                  counter_total = counter_total + 1                  
                    
                  !cluster calculations  
                  n_contacts_per_molecule(m,k) = n_contacts_per_molecule(m,k) + 1             
                  IF (n_contacts_per_molecule(m,k) > max_contacts) STOP "There are more than max_contacts particles around one"
                  node_next(m,n_contacts_per_molecule(m,k),k) = n
                  n_links(k) = n_links(k) + 1
                  EXIT loop
                END IF              
              END DO     
            END DO loop         
    END DO 
  END DO
  
!  IF ((counter_total/2) /= 0) THEN
!    DO index_proba = 1,n_bins_angle
!      proba_angle(index_proba) = REAL(count_bin_angle(index_proba),8)/REAL(counter_total,8)                  
!    END DO 
!    d = d/REAL(counter_d,8)
!    d_ave = d_ave + d
!  END IF  
!  
!  WRITE(5,*) "Frame ", k-1  
!  WRITE(5,*) "Number of contacts", counter_total/2
!  WRITE(5,*) "Number of parallel contacts", counter_d/2
!  WRITE(5,*) "Cutoff Distance", dist_cluster
!  WRITE(5,*) "Average separation between molecule planes", d
!  WRITE(5,'(a9,5x,a20,7x,a20,9x,a11)') "bin angle", 'lower limit(degrees)', 'upper limit(degrees)' ,"Probability"  
!  
!  DO index_proba = 1,n_bins_angle
!    WRITE(5,*) index_proba, proba_incre*(index_proba-1)*180/pi, proba_incre*(index_proba)*180/pi , &
!    proba_angle(index_proba) 
!    proba_ave(index_proba) = proba_ave(index_proba) + proba_angle(index_proba)
!  END DO 
!  WRITE(5,*);WRITE(5,*);WRITE(5,*)      
  
  WRITE(2,*) "Number of de Links ", n_links(k) / 2  ! Links counted twice in the loop  
  DO i = 1, n_molecules
    WRITE(2,*) "node",i
    WRITE(2,*) "neighbors : ", (node_next(i,j,k) , j=1,n_contacts_per_molecule(i,k)) 
    WRITE(2,*)
  END DO 
  WRITE(2,*); WRITE(2,*);  WRITE(2,*)
  
  END DO ! End algorithm create nodes 
  
!  WRITE(5,*) "Average Probability"
!  WRITE(5,*) "Average total distance",d_ave/REAL(nframes,8)
!  WRITE(5,'(a9,5x,a20,7x,a20,9x,a11)') "bin angle", 'low limit(degrees)', 'up limit(degrees)' ,"Probability"  
!  DO index_proba = 1,n_bins_angle
!    WRITE(5,*) index_proba, proba_incre*(index_proba-1)*180/pi, proba_incre*(index_proba)*180/pi , &
!    proba_ave(index_proba)/ REAL(nframes,8)
!  END DO  
  
  CLOSE(2)
  CLOSE(5)
  CLOSE(7)
  !write(*,*) dist_cluster, n_links(nframes)/2
  
  CALL HK(node_next,n_contacts_per_molecule,n_links)
  
  !end do
  
  DEALLOCATE(node_next, n_contacts_per_molecule)  
END SUBROUTINE 

SUBROUTINE HK(node_next,n_contacts_per_molecule,n_links)
  IMPLICIT NONE
  INTEGER :: i, j, k, status_open
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: n_contacts_per_molecule
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: node_next
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nodeL
  INTEGER, DIMENSION(:), ALLOCATABLE :: labels
  INTEGER, DIMENSION(n_molecules,n_molecules) :: clusters  ! Maximum of n_molecules nº of clusters, extreme case = 1 cluster with all molecules together
  INTEGER, DIMENSION(n_molecules) :: id, n_mol_per_cluster  ! Maximum of  n_molecules nº of clusters
  INTEGER, DIMENSION(nframes):: n_cluster  
  INTEGER, DIMENSION(nframes) :: n_links, nclus
  INTEGER, DIMENSION(n_molecules,2) :: cluster_distribution ! Number of cluster with i molecules and the second column saves the i value for statistic
  REAL(8), DIMENSION(n_molecules) :: cluster_probability ! Probability of finding a cluster with i molecules
  INTEGER :: N, tol = 1000
  LOGICAL :: condition,condition1,condition2
  REAL(8), DIMENSION(nframes) :: mean_size
  REAL(8) :: total_size,start,finish
  ALLOCATE(nodeL(n_molecules,nframes), labels(max_contacts))
  
  WRITE(*,*) "Running clustering routine"
  WRITE(*,*) "Remember : The maximum number of cluster is 1000"
  
  !HK algorithm implementation
  nodeL = 0     !node label
  n_cluster = 0
  nclus = 0
  mean_size = 0.0
  ! node = molecule
  ! nodeL starts corresponding to molecules but finishes corresponding to cluster label  
  
  !Exibiting nodeL
  OPEN(unit=1,file="clustering.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening clustering.out"

  OPEN(unit=3,file="cluster_statistics.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening cluster_statistics.out"  

  OPEN(unit=4,file="cluster_time_evolution.out",status="replace",iostat = status_open ) 
  IF (status_open > 0) STOP "error opening cluster_time_evolution.out"  
  WRITE(4,"(1x,a5,12x,a4,2x,a11,2x,a8,2x,a20)") "Frame", "time", "Nº clusters", "Nº links", "Average cluster size"
  
  !WRITE(*,*); WRITE(*,*); WRITE(*,*) "Elapsed time for clustering routine"
  !CALL cpu_time(start)
  
  DO k = 1, nframes
      WRITE(*,*) "Running frame ",k
      
  DO i = 1, n_molecules
      !If this node doesn't have any neighbor, it should be a single cluster
      IF (n_contacts_per_molecule(i,k) == 0 ) THEN
        n_cluster(k) = n_cluster(k) + 1    
        nodeL(i,k) = n_cluster(k)
      ELSE 
        !Testing if all neighbors of this node don't have labels as well  
        labels = 0  
        DO j = 1, n_contacts_per_molecule(i,k)
          labels(j) = nodeL(node_next(i,j,k),k)  
        END DO
        
        IF (MAXVAL(labels) == 0) THEN
          n_cluster(k) = n_cluster(k) + 1   
          nodeL(i,k) = n_cluster(k)
        ELSE
          nodeL(i,k) = MINVAL(labels, MASK = labels > 0)
          DO j = 1, n_contacts_per_molecule(i,k)
            nodeL(node_next(i,j,k),k) = nodeL(i,k)  
          END DO          
        END IF
          
      END IF     
  END DO 

  N = 1   
  condition = .TRUE.
  DO WHILE ( (N < tol) .AND. (condition .EQV. .TRUE.) )
  
  condition = .FALSE.
  !Checking if all connected particles have the same cluster label
         DO i = 1,n_molecules
           IF (n_contacts_per_molecule(i,k) /= 0) THEN   
             labels = 0   
             DO j = 1, n_contacts_per_molecule(i,k)
               labels(j) = nodeL(node_next(i,j,k),k)  
             END DO 
           
             nodeL(i,k) = MINVAL(labels, MASK = labels > 0)
             DO j = 1, n_contacts_per_molecule(i,k)
               nodeL(node_next(i,j,k),k) = nodeL(i,k)  
             END DO            
           END IF
         END DO   
         
  loop2: DO i = 1, n_molecules  
           labels = 0   
           DO j = 1, n_contacts_per_molecule(i,k)
             labels(j) = nodeL(node_next(i,j,k),k)  
           END DO 
           DO j = 1, n_contacts_per_molecule(i,k) - 1
             IF (labels(j) /= labels(j+1)) THEN
               N = N + 1  
               condition = .TRUE.
               EXIT loop2 
             END IF
           END DO           
         END DO loop2 
           
  END DO    ! End do while  
  
  id = 0
  !Finding number of clusters
  DO i = 1, n_molecules
      condition1 = .FALSE.
      
      !If i = n_molecules this molecules(node) cannot form aggregate forth and condition1 = .FALSE. and it will be the opportunity
      ! to be a single cluster (if condition2 = .TRUE.)
      IF (i < n_molecules) THEN
        DO j = i+1, n_molecules   
          IF (nodeL(i,k) == nodeL(j,k) ) THEN
            condition1 = .TRUE.
          END IF
        END DO
      END IF
      
      condition2 = .TRUE. 
      DO j = i, 1, -1   
        IF ((nodeL(i,k) == nodeL(j,k) ) .AND. (i /= j) ) THEN
            condition2 = .FALSE.
        END IF
      END DO  
      
      IF ((condition1 .EQV. .TRUE. ) .AND. (condition2 .EQV. .TRUE. )) THEN
        nclus(k) = nclus(k) + 1   !Cluster of many particles
        id(nclus(k)) = nodeL(i,k)     
      END IF    
      
      IF ((condition1 .EQV. .FALSE. ) .AND. (condition2 .EQV. .TRUE. )) THEN
        nclus(k) = nclus(k) + 1  !Cluster of one particle
        id(nclus(k)) = nodeL(i,k)
      END IF
          
  END DO 

  !CALL cpu_time(finish)
  !PRINT '("Time for clustering = ",f12.9," seconds.")', finish-start 

  
  n_mol_per_cluster = 0
  clusters = 0
  DO i = 1, nclus(k)
    DO j = 1, n_molecules  
      IF ( id(i) == nodeL(j,k)) THEN  
        n_mol_per_cluster(i) = n_mol_per_cluster(i) + 1  
        clusters(i,j) = j         
      END IF    
    END DO
  END DO
  
  !Writing dump file
  WRITE(1,*) "Frame ",k  
  WRITE(1,*) "Numbers of Links ", n_links(k)/2  ! Links counted twice in the loop    
  WRITE(1,*) "Nº clusters", nclus(k)
  WRITE(1,*) "Nº iterations for convergence",N 
  WRITE(1,*) 
  DO i = 1, nclus(k)
    WRITE(1,*) "Cluster : ",i
    WRITE(1,"(1x,a11,i4,a4)") "Molecules (",n_mol_per_cluster(i),") : "
    DO j = 1,n_molecules
      IF (clusters(i,j) > 0) THEN  
        WRITE(1,*) clusters(i,j) 
      END IF
    END DO  
    WRITE(1,*)  
  END DO  
  WRITE(1,*);WRITE(1,*);WRITE(1,*)    
  
  !Cluster statistics  
  
  ! Cluster Distribution
  cluster_distribution = 0
  DO i = 1, n_molecules
    DO j = 1,nclus(k)  
      IF (n_mol_per_cluster(j) == i) THEN
        cluster_distribution(i,1) =  cluster_distribution(i,1) + 1  ! -> ns(t)
        cluster_distribution(i,2) =  i                              ! -> s
      END IF
    END DO
  END DO
  
  total_size = 0.0
  cluster_probability = 0.0
  
  DO i = 1, n_molecules
    total_size = total_size + REAL(cluster_distribution(i,1) * cluster_distribution(i,2),8)  
  END DO    
  
  !Probabilities
  DO i = 1, n_molecules
    cluster_probability(i) =  REAL(cluster_distribution(i,1) * cluster_distribution(i,2),8) / total_size   
    mean_size(k) = mean_size(k) + cluster_probability(i) *  REAL(cluster_distribution(i,2),8)
  END DO
  
  WRITE(3,*) "Frame ",k
  WRITE(3,"(4x,a25,1x,a19,1x,a12)") "Nº molecules per Cluster ", "Number of clusters ", " Probability"
  ! Calculating Histogram
  DO i = 1, n_molecules
    IF (cluster_distribution(i,1) /= 0) THEN
      WRITE(3,"(1x,i4,22x,i4,19x,F10.8)") i, cluster_distribution(i,1), cluster_probability(i)
    END IF    
  END DO
  WRITE(3,*)
  WRITE(3,*) "Mean size of clusters in this frame ", mean_size(k)
  WRITE(3,*);WRITE(3,*)
  
  !Writing time evolution file
  WRITE(4,"(1x,i4,2x,F15.3,2x,i4,5x,i5,2x,F10.3)") k-1, (k-1)*delta_t*dump_interval, nclus(k), n_links(k)/2, mean_size(k)
 
  END DO    ! End nframes's do 
  
  CLOSE(1)
  CLOSE(3)
  CLOSE(4)
  DEALLOCATE(nodeL,labels)
END SUBROUTINE HK      


END MODULE cluster_analysis
