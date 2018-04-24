program tarea01
  use parametros, only: pi, dim, k1, k2, dx, Long
  use diferencias_finitas
  implicit none
  integer :: i, j, ii
  integer, parameter :: szi=3
  integer, parameter, dimension(szi) :: scl = (/ 1, 2, 5 /)
  real(Long), dimension(:), allocatable :: x, f, df, d2f
  real(Long), dimension(:), allocatable :: df_2, df_4, df_6
  real(Long), dimension(:), allocatable :: d2f_2, d2f_4, d2f_6
  real(Long), dimension(1:3) :: K_2
  real(Long), dimension(1,3) :: Kb_2
  real(Long), dimension(1:5) :: K_4
  real(Long), dimension(2,5) :: Kb_4
  real(Long), dimension(1:7) :: K_6
  real(Long), dimension(3,7) :: Kb_6
  real(Long) :: err_df_2 = 0.0_Long, err_df_4 = 0.0_Long, err_df_6 = 0.0_Long
  real(Long) :: err_d2f_2 = 0.0_Long, err_d2f_4 = 0.0_Long, err_d2f_6 = 0.0_Long

  do i=1,6
     do ii = 1,szi
        dim = scl(ii)*10**i+1
        dx = 2.0*pi/(dim-1);
        k1 = 0.5*pi/1.0; k2 = 0.75*k1
        
        allocate( x(dim) )
        allocate( f(dim) )
        allocate( df(dim) )
        allocate( d2f(dim) )

        allocate( df_2(dim) )
        allocate( df_4(dim) )
        allocate( df_6(dim) )

        allocate( d2f_2(dim) )
        allocate( d2f_4(dim) )
        allocate( d2f_6(dim) )

        df_2 = 0.0_Long; df_4 = 0.0_Long; df_6 = 0.0_Long
        d2f_2 = 0.0_Long; d2f_4 = 0.0_Long; d2f_6 = 0.0_Long
        
        call get_Kernel(Kb_2(1,:), 3, 0, dx)
        call get_Kernel(K_2, 3, -1, dx)

        call get_Kernel(Kb_4(1,:), 5, 0, dx)
        call get_Kernel(Kb_4(2,:), 5,-1, dx)
        call get_Kernel(K_4, 5, -2, dx)

        call get_Kernel(Kb_6(1,:), 7, 0, dx)
        call get_Kernel(Kb_6(2,:), 7,-1, dx)
        call get_Kernel(Kb_6(3,:), 7,-2, dx)
        call get_Kernel(K_6, 7, -3, dx)

        x(1:dim) = (/ (2*pi*(i-1)/(dim-1)-pi, i=1,dim) /)

        f = exp(sin(x)) + 0.5*cos(k1*x) - 0.8*sin(k2*x)
        df = cos(x)*exp(sin(x)) - 0.5*k1*sin(k1*x) - 0.8*k2*cos(k2*x)
        d2f = (cos(x)**2-sin(x))*exp(sin(x)) - 0.5*k1**2*cos(k1*x) + 0.8*k2**2*sin(k2*x)

        call convolution(f,dim,K_2,3,Kb_2,3,1,df_2)
        call convolution(f,dim,K_4,5,Kb_4,5,2,df_4)
        call convolution(f,dim,K_6,7,Kb_6,7,3,df_6)
        call convolution(df_2,dim,K_2,3,Kb_2,3,1,d2f_2)
        call convolution(df_4,dim,K_4,5,Kb_4,5,2,d2f_4)
        call convolution(df_6,dim,K_6,7,Kb_6,7,3,d2f_6)

        err_df_2 = maxval(abs(df-df_2))
        err_df_4 = maxval(abs(df-df_4))
        err_df_6 = maxval(abs(df-df_6))
        err_d2f_2 = maxval(abs(d2f-d2f_2))
        err_d2f_4 = maxval(abs(d2f-d2f_4))
        err_d2f_6 = maxval(abs(d2f-d2f_6))

        write(*,*) dim, dx, err_df_2, err_df_4, err_df_6, &
             maxval( (/ (abs(df((dim-1)*i/10)-df_2((dim-1)*i/10)), i=3,7) /) ), &
             maxval( (/ (abs(df((dim-1)*i/10)-df_4((dim-1)*i/10)), i=3,7) /) ), &
             maxval( (/ (abs(df((dim-1)*i/10)-df_6((dim-1)*i/10)), i=3,7) /) ), &
             err_d2f_2, err_d2f_4, err_d2f_6, &
             maxval( (/ (abs(d2f((dim-1)*i/10)-d2f_2((dim-1)*i/10)), i=3,7) /) ), &
             maxval( (/ (abs(d2f((dim-1)*i/10)-d2f_4((dim-1)*i/10)), i=3,7) /) ), &
             maxval( (/ (abs(d2f((dim-1)*i/10)-d2f_6((dim-1)*i/10)), i=3,7) /) )

        if (allocated(x)) deallocate( x )
        if (allocated(f)) deallocate( f )
        if (allocated(df)) deallocate( df )
        if (allocated(d2f)) deallocate( d2f )

        if (allocated(df_2)) deallocate( df_2 )
        if (allocated(df_4)) deallocate( df_4 )
        if (allocated(df_6)) deallocate( df_6 )

        if (allocated(d2f_2)) deallocate( d2f_2 )
        if (allocated(d2f_4)) deallocate( d2f_4 )
        if (allocated(d2f_6)) deallocate( d2f_6 )
     end do
  end do

contains

  function f_fcn(x)
    use parametros, only: k1, k2
    real(Long) :: x
    real(Long) :: f_fcn
    f_fcn = exp(sin(x)) + 0.5*cos(k1*x) - 0.8*sin(k2*x)
  end function f_fcn
  
  function df_fcn(x)
    use parametros, only: k1, k2
    real(Long) :: x
    real(Long) :: df_fcn
    df_fcn = cos(x)*exp(sin(x)) - 0.5*k1*sin(k1*x) - 0.8*k2*cos(k2*x)
  end function df_fcn
  
  function d2f_fcn(x)
    use parametros, only: k1, k2
    real(Long) :: x
    real(Long) :: d2f_fcn
    d2f_fcn = (cos(x)**2-sin(x))*exp(sin(x)) - 0.5*k1**2*cos(k1*x) + 0.8*k2**2*sin(k2*x)
  end function d2f_fcn
  
end program tarea01
