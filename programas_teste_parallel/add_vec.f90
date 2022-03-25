PROGRAM VEC_ADD_SECTIONS
!$ use omp_lib
INTEGER I,tid,NTHREADS
integer ,PARAMETER :: N=100
REAL A(N), B(N), C(N)



!Some initializations
DO I = 1, N
A(I) = I * 1.0
B(I) = A(I)
ENDDO
!$OMP PARALLEL SHARED(A,B,C), PRIVATE(I)
!$OMP SECTIONS
!$OMP SECTION
DO I = 1, N/2
C(I) = A(I) + B(I)
tid = OMP_GET_THREAD_NUM()
NTHREADS = OMP_GET_NUM_THREADS()
write(*,*) i,c(i),tid, nthreads
ENDDO
!$OMP SECTION
DO I = 1+N/2, N
C(I) = A(I) + B(I)
tid = OMP_GET_THREAD_NUM()
NTHREADS = OMP_GET_NUM_THREADS()
write(*,*) i,c(i),tid, nthreads
ENDDO
!$OMP END SECTIONS  nowait  
!$OMP END PARALLEL

!esse nowait parece deixar mais rápido, ele cria o vetor c mais rápido pois um processador não espera o outro, porém pra mostrar ele mistura
!justamente por um processador não esperar o outro



END
