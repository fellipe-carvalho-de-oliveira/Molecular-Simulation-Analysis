PROGRAM VECTOR_ADD
!$ use omp_lib
INTEGER  I,NTHREADS,tid
integer,PARAMETER ::N=1000
integer, PARAMETER ::CHUNK=4
REAL A(N), B(N), C(N)
!Some initializations
DO I = 1, N
A(I) = I * 1.0
B(I) = A(I)
ENDDO
!$OMP PARALLEL DO &
!$OMP SHARED(A,B,C) PRIVATE(I) &
!$OMP SCHEDULE(dynamic,CHUNK)
DO I = 1, N
C(I) = A(I) + B(I)
tid = OMP_GET_THREAD_NUM()
NTHREADS = OMP_GET_NUM_THREADS()
write(*,*)  i, c(i),tid, nthreads
ENDDO


!$OMP END PARALLEL DO


END
