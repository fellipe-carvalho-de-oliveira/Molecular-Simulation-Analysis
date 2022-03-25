!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   cluster_analysis.f90
! Author: fellipe
!
! Created on 21 de Fevereiro de 2019, 08:49
!

MODULE cluster_analysis
USE global_variables

CONTAINS

SUBROUTINE clustering
  IMPLICIT NONE    
  INTEGER :: k,i,j,m,n,n_links, count_links = 0
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: node, node_s, link, link_s, n_contacts_per_molecule
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: node_next, links_of_node
  REAL(8) :: dist, dx , dy , dz 
  
  !Maximun of max_contacts particles attached to one particle, so there are max_contacts maximum links
  ALLOCATE(node(n_molecules,nframes), node_s(n_molecules,nframes), node_next(n_molecules,max_contacts,nframes))
  ALLOCATE(n_contacts_per_molecule(n_molecules,nframes))
  
  node = 0
  node_s = 0
  link = 0
  link_s = 0
  node_next = 0
  links_of_node = 0
  n_contacts_per_molecule = 0
  
  !beginning algorithm create nodes,links vector
  DO k = 1,nframes
  
  DO i = 1,natoms
    node(i,k) = i   !index for all nodes (molecules)
    node_s(i,k) = 1 !All nodes occupied
  END DO    
  
  n_links = 0
  DO m = 1,n_molecules-1 !loop for number of atoms
    DO n = m+1,n_molecules   
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
          
                IF (dist < dist_cluster * dist_cluster) THEN
                  n_contacts_per_molecule(m,k) = n_contacts_per_molecule(m,k) + 1             
                  IF (n_contacts_per_molecule(m,k) > max_contacts) STOP "There are more than max_contacts particles around one"
                
                  node_next(m,n_contacts_per_molecule(m,k),k) = n
                  n_links = n_links + 1
                  EXIT loop
                END IF              
              END DO     
            END DO loop         
    END DO 
  END DO
  
  ! Treating the Links
  ALLOCATE(link(n_links,nframes), link_s(n_links,nframes), links_of_node(n_links,max_contacts,nframes))

  DO i = 1,n_links
    link(i,k) = i   !index for all links
    link_s(i,k) = 1 !All links occupied
  END DO   
  
  count_links = 0
  DO i = 1, n_molecules
    DO j = 1, n_contacts_per_molecule(i,k)
      count_links = count_links +1  
      links_of_node(i,j,k) = count_links     
    END DO
  END DO
  
  END DO ! End algorithm create nodes and links vectors    
  
  CALL HK(node, node_s, link, link_s, node_next, links_of_node,n_contacts_per_molecule)
  
  DEALLOCATE(node, node_s, link, link_s, node_next, links_of_node,n_contacts_per_molecule)  
END SUBROUTINE 

SUBROUTINE HK(node, node_s, link, link_s, node_next, links_of_node,n_contacts_per_molecule)
  IMPLICIT NONE
  INTEGER :: n_clusters, i, j, k
  INTEGER, DIMENSION(:,:), ALLOCATABLE, INTENT(INOUT) :: node, node_s, link, link_s, n_contacts_per_molecule
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE, INTENT(INOUT) :: node_next, links_of_node
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nodeL, linkL, nodeLP
  INTEGER, DIMENSION(:), ALLOCATABLE :: labels
  INTEGER :: n_cluster , N, tol = 10000
  LOGICAL :: condition = .TRUE.
  ALLOCATE(nodeL(n_molecules,nframes) , linkL(n_molecules,nframes), labels(max_contacts))
  
  !HK algorithm implementation
  nodeL = 0     !node label
  linkL = 0     !link label
  
  ! node = molecule
  ! nodeL starts corresponding to molecules but finishes corresponding to cluster label  
  
  DO k = 1, nframes
  N = 0    
  DO WHILE ( (N < tol) .OR. (condition .EQV. .TRUE.) )
  
  n_cluster = 0
  DO i = 1, n_molecules
      !If this node doesn't have any neighbor, it should be a cluster alone
      IF (n_contacts_per_molecule(i,k) == 0 ) THEN
        n_cluster = n_cluster + 1    
        nodeL(i,k) = n_cluster
      ELSE 
        !Testing if all neighbors of this node don't have labels as well  
        labels = 0  
        DO j = 1, n_contacts_per_molecule(i,k)
          labels(j) = nodeL(node_next(i,j,k),k)  
        END DO
        
        IF (MAXVAL(labels) == 0) THEN
          n_cluster = n_cluster + 1   
          nodeL(i,k) = n_cluster
        ELSE
          nodeL(i,k) = MINVAL(labels, MASK = labels > 0)
          DO j = 1, n_contacts_per_molecule(i,k)
            nodeL(node_next(i,j,k),k) = nodeL(i,k)  
          END DO          
        END IF
          
      END IF    
  END DO 
  
  condition = .FALSE.
  !Checking if all connected particles have the same cluster label
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
  END DO    ! End nframes's do 
  
  !Exibiting nodeL
  DO k = 1, nframes
    WRITE(*,*) "Frame ",i   
    DO i = 1, n_molecules
      WRITE(*,*) nodeL(i,k)   
    END DO
    WRITE(*,*);WRITE(*,*) 
    READ(*,*)
  END DO
  
  DEALLOCATE(nodeL,linkL,nodeLP,labels)
END SUBROUTINE HK    


END MODULE cluster_analysis
