!
! To change this license header, choose License Headers in Project Properties.
! To change this template file, choose Tools | Templates
! and open the template in the editor.
!

!     
! File:   main.f90
! Author: fellipe
!
! Created on 1 de Março de 2019, 18:08
!

program main
!    !Exemplo1
!    implicit none
!    real:: start,finish
!    integer, parameter :: n=1000
!    integer :: i,a(n),j,k,b=0
!    
!    call cpu_time(start)
!    
!    !$omp parallel 
!    do i=1,n  
!        do j=1,n
!          do k=1,n    
!            b=b+1
!          enddo
!        enddo
!    enddo
!    !$omp end parallel 
!    write(*,*) "olhe bem nos meus olhos",b
    
!    !Example 2 
!    !Example of how to use thread numbers
!    !$ use omp_lib
!    implicit none
!    integer :: tid
!    real:: start,finish
!    integer :: i,j,k,b=0,n=1000
!    call cpu_time(start)
!    
!    !$omp parallel
!    tid=2 !Dummy thread number for the case OMP is not enabled
!    !$ tid=omp_get_thread_num()
!    write (6,*) 'Doing some stuff in thread ',tid
!    do i=1,n  
!        do j=1,n
!          do k=1,n    
!            b=b+1
!          enddo
!        enddo
!    enddo
!    !$omp end parallel

!    !Example 3 do loop shared
!    !Example of a parallel loop
!    !$ use omp_lib
!    implicit none
!    real:: start,finish
!    integer,parameter :: n=1000000
!    integer :: i,a(n),t(n),tid
    
!    call cpu_time(start)
    
!    !$omp parallel shared(a,t) private(i,tid)
!    tid=0 !Dummy thread number, for the case OMP is not enabled
!    !$ tid=omp_get_thread_num()
!    !$omp do
!    do i=1   ,n  
!        a(i)=2*i
!        t(i)=tid !Record which thread did which iteration
!    enddo
!    !$omp end do
!    !$omp end parallel
!    !Show what was produced
!    write(6,*) 'i,a(i),thread number'
!    do i=1   ,n  
!        write(6,*) i,a(i),t(i)
!    enddo

!     	!Example 4 sections
!     	!Example of sections
!     	!$ use omp_lib
!     	implicit none
!	real:: start,finish
!     	integer,parameter :: n=7
!	integer :: i,a(n),ta(n),b(n),tb(n),tid
!	call cpu_time(start)

!	!$omp parallel shared(a,ta,b,tb) private(i,tid)
!	!$omp sections
!	!$omp section
!	tid  =0 !Dummy thread number, for the case OMP is not enabled
!	!$ tid=omp_get_thread_num()
!	do i=1,n
!		a(i)=2*i
!		ta(i)=tid !Record which thread did which iteration
!	enddo
!	!$omp section
!	tid  =0 !Dummy thread number, for the case OMP is not enabled
!	!$ tid=omp_get_thread_num()
!	do i=1,n
!		b(i)=3*i + a(i) ! as vezes dá tempo de ele somar, se o processador de cima for mais rápido em gerar a(i)
!		tb(i)=tid !Record which thread did which iteration
!	enddo
!	!$omp end sections
!	!$omp end parallel
!	!Show what was produced
!	write(6,*) 'i,a(i),thread that made a(i),b(i), thread that made b(i)'
!	do i=1,n
!		write(6,*) i,a(i),ta(i),b(i),tb(i) 
!	enddo

!Example $ race condition    
!Example of a race condition
!$ use omp_lib
implicit none
integer :: i,n,a
real:: start,finish

call cpu_time(start)
n=1000000
a=0
!$omp parallel shared(a,n) private(i)
!$omp do
do i=1,n	!When one thread writes to a, the value of a!it had read might not be current anymore  
a=a+1
enddo
!$omp end do
!$omp end parallel
!Show what was produced
write(6,*) 'a=',a

    
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
    
end
