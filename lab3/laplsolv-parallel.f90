program laplsolv
    use omp_lib

    !-----------------------------------------------------------------------
    ! Serial program for solving the heat conduction problem 
    ! on a square using the Jacobi method. 
    ! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
    ! Modified by Berkant Savas (besav@math.liu.se) April 2006
    ! Paralleled by Hao-Hsiang Liao May 2018
    ! Modified by Santi Pagola May 2018
    !-----------------------------------------------------------------------
    integer, parameter                        :: n=1000, maxiter=1000
    double precision,parameter                :: tol=1.0E-3
    double precision,dimension(0:n+1,0:n+1)   :: T
    double precision                          :: error
    double precision                          :: t1,t0
    integer                                   :: i,j,k
    integer                                   :: nr_threads
    integer                                   :: ratio,my_id
    integer,dimension(0:99)                   :: start_pos,stop_pos ! Initialize with a big number
    character(len=20)                         :: str
    double precision,dimension(n)             :: tmp2
    double precision,dimension(0:n-1,0:99)    :: tmp1,tmp3

    ! Set boundary conditions and initial values for the unknowns
    T=0.0D0
    T(0:n+1 , 0)     = 1.0D0
    T(0:n+1 , n+1)   = 1.0D0
    T(n+1   , 0:n+1) = 2.0D0

    !$omp parallel shared(nr_threads)
    nr_threads    = OMP_GET_NUM_THREADS()
    !$omp end parallel


    ! Determine how many stuffs should be implemented by a thread
    ratio = n/nr_threads

      ! Determine every thread's start/stop position
     do i=0,nr_threads-1
        start_pos(i) = (i*ratio)+1
        stop_pos(i) = (i+1)*ratio
     end do

  ! Solve the linear system of equations using the Jacobi method
  t0 = omp_get_wtime()
  
  do k=1,maxiter
     
     !tmp1=T(1:n,0)
     error=0.0D0

     ! We need to assign various start_pos/stop_pos's previous one array for each thread, their relationship is tmp1=T(1:n,j-1) when tmp2=T(1:n,j)
     ! tmp3 is an array that consist original stop_pos's value
     do i=0,nr_threads-1
        tmp1(0:n-1, i) = T(1:n, start_pos(i)-1)
        tmp3(0:n-1, i) = T(1:n, stop_pos(i)+1)
	 end do

     !$omp parallel private(j,tmp2,my_id) shared(T,tmp1,tmp3,k) reduction(max:error)
     my_id = OMP_GET_THREAD_NUM()

     do j=start_pos(my_id),(stop_pos(my_id)-1)
        tmp2=T(1:n,j)
        T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+tmp1(0:n-1,my_id))/4.0D0
        error=max(error,maxval(abs(tmp2-T(1:n,j))))
        tmp1(0:n-1,my_id)=tmp2
     end do
     
     tmp2=T(1:n,stop_pos(my_id))
     T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+tmp3(0:n-1,my_id)+tmp1(0:n-1,my_id))/4.0D0
     error=max(error,maxval(abs(tmp2-T(1:n,j))))
     tmp1(0:n-1,my_id)=tmp2
     
     !$omp end parallel

     if (error<tol) then
        exit
     end if
     
  end do
  
  t1 = omp_get_wtime()

  print *, 'The number of threads is',nr_threads
  write(unit=*,fmt=*) 'Time:',t1-t0,'Number of Iterations:',k
  write(unit=*,fmt=*) 'Temperature of element T(1,1)  =',T(1,1)

  ! Uncomment the next part if you want to write the whole solution
  ! to a file. Useful for plotting. 
  
  !open(unit=7,action='write',file='result.dat',status='unknown')
  !write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
  !do i=0,n+1
  !   write (unit=7,fmt=str) T(i,0:n+1)  
  !end do
  !close(unit=7)
  
end program laplsolv
