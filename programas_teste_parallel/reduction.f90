PROGRAM DOT_PRODUCT
!$ use omp_lib
INTEGER I,NTHREADS,tid
real:: start,finish
integer, PARAMETER :: N=100
integer, PARAMETER :: CHUNK=10
REAL A(N), B(N), RESULT

!result = 1

call cpu_time(start)

!Some initializations
DO I = 1, N
A(I) = I * 1.0
B(I) = I * 2.0
ENDDO
RESULT= 0.0
!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) PRIVATE(I) &
!$OMP SCHEDULE(STATIC,CHUNK) &
!$OMP REDUCTION(+:RESULT)
DO I = 1, N
RESULT = RESULT + (A(I) * B(I)) 
!tid = OMP_GET_THREAD_NUM()
!NTHREADS = OMP_GET_NUM_THREADS()
!write(*,*)  i,tid, nthreads
ENDDO
!$OMP END PARALLEL DO
PRINT *, "Final Result= ", RESULT

call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
END
