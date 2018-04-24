program diferencias
  use parametros, only: pi, dim, dx, k1, k2, Long
  use diferencias_finitas
  implicit none
  integer :: i, j, k
  integer, parameter :: szi=3
  integer, parameter, dimension(szi) :: scl = (/ 1, 2, 5 /)
  real(Long), dimension(1:3) :: K_2
  real(Long), dimension(1,3) :: Kb_2
  real(Long), dimension(1:5) :: K_4
  real(Long), dimension(2,5) :: Kb_4
  real(Long), dimension(1:7) :: K_6
  real(Long), dimension(3,7) :: Kb_6
  real(Long) :: err_df_2 = 0.0_Long, err_df_4 = 0.0_Long, err_df_6 = 0.0_Long
  real(Long) :: err_b_df_2 = 0.0_Long, err_b_df_4 = 0.0_Long, err_b_df_6 = 0.0_Long
  real(Long) :: err_d2f_2 = 0.0_Long, err_d2f_4 = 0.0_Long, err_d2f_6 = 0.0_Long
  real(Long) :: err_b_d2f_2 = 0.0_Long, err_b_d2f_4 = 0.0_Long, err_b_d2f_6 = 0.0_Long
  real(Long) :: dx_max = 1.0_Long
  
  do i=7,8
     do j=1,szi
        dim = scl(j)*10**i+1
        dx = 2.0_Long*pi/(dim-1);
        k1 = 0.5_Long*pi/dx_max; k2 = 0.75_Long*k1

        !Kernels de contorno
        call get_Kernel(Kb_2(1,:), 3, 0, dx)

        call get_Kernel(Kb_4(1,:), 5, 0, dx)
        call get_Kernel(Kb_4(2,:), 5,-1, dx)

        call get_Kernel(Kb_6(1,:), 7, 0, dx)
        call get_Kernel(Kb_6(2,:), 7,-1, dx)
        call get_Kernel(Kb_6(3,:), 7,-2, dx)
        ! Kernels de dominio
        call get_Kernel(K_2, 3, -1, dx)
        call get_Kernel(K_4, 5, -2, dx)
        call get_Kernel(K_6, 7, -3, dx)
        !call errors_kernel(dim,K_domain,Kd_cols,K_boundary,Kb_cols,Kb_rows,err_max,err_b_max,errd2_max,errd2_b_max)
        call errors_kernel(dim,K_2,3,Kb_2,3,1,err_b_df_2,err_df_2,err_b_d2f_2,err_d2f_2)
        call errors_kernel(dim,K_4,5,Kb_4,5,2,err_b_df_4,err_df_4,err_b_d2f_4,err_d2f_4)
        call errors_kernel(dim,K_6,7,Kb_6,7,3,err_b_df_6,err_df_6,err_b_d2f_6,err_d2f_6)        
        ! postprocess
        write(*,*) dim, dx, &
             err_b_df_2, err_b_df_4, err_b_df_6, &
             err_df_2, err_df_4, err_df_6, &
             err_b_d2f_2, err_b_d2f_4, err_b_d2f_6, &
             err_d2f_2, err_d2f_4, err_d2f_6
     end do
  end do

contains

  ! x: is an index in [1:dim]
  function fun(x)
    use parametros, only: k1, k2, dim, Long, pi
    real(Long) :: x
    real(Long) :: fun
    x = 2*pi*(x-1)/(dim-1)-pi
    fun = exp(sin(x)) + 0.5*cos(k1*x) - 0.8*sin(k2*x)
  end function fun
  
  function dfun(x)
    use parametros, only: k1, k2, dim, Long, pi
    real(Long) :: x
    real(Long) :: dfun
    x = 2*pi*(x-1)/(dim-1)-pi
    dfun = cos(x)*exp(sin(x)) - 0.5*k1*sin(k1*x) - 0.8*k2*cos(k2*x)
  end function dfun

  function d2fun(x)
    use parametros, only: k1, k2, dim, Long, pi
    real(Long) :: x
    real(Long) :: d2fun
    x = 2*pi*(x-1)/(dim-1)-pi
    d2fun = (cos(x)**2-sin(x))*exp(sin(x)) - 0.5*k1**2*cos(k1*x) + 0.8*k2**2*sin(k2*x)
  end function d2fun

  subroutine push_back(v,e,dim)
    use parametros, only: Long
    implicit none
    real(Long), dimension(:), intent(inout) :: v
    real(Long), intent(in) :: e
    integer, intent(in) :: dim
    integer :: i
    do i=1,dim-1
       v(i) = v(i+1)
    end do
    v(dim) = e
  end subroutine push_back

  subroutine errors_kernel(dim,K_domain,Kd_cols,K_boundary,Kb_cols,Kb_rows,err_max,err_b_max,errd2_max,errd2_b_max)
    use parametros, only: Long, pi, k1, k2
    implicit none
    real(Long), intent(in), dimension(:) :: K_domain
    real(Long), intent(in), dimension(:,:) :: K_boundary
    integer, intent(in) :: dim, Kd_cols, Kb_cols, Kb_rows
    real(Long), intent(out) :: err_max, err_b_max
    real(Long), intent(out) :: errd2_max, errd2_b_max
    real(Long), dimension(Kd_cols) :: x_tmp, f_tmp, df_tmp
    real(Long) :: err_i, errd2_i 
    real(Long) :: fun_i, dfun_i, d2fun_i
    integer :: i, j, k
    err_max = -1; err_b_max = -1
    errd2_max = -1; errd2_b_max = -1

    x_tmp(1:Kb_cols) = (/ (2*pi*(i-1)/(dim-1)-pi, i=1,Kb_cols) /)
    f_tmp = exp(sin(x_tmp)) + 0.5*cos(k1*x_tmp) - 0.8*sin(k2*x_tmp)

    ! left boundary
    do i=1,Kb_rows
       ! first derivative
       dfun_i = dot_product(K_boundary(i,:),f_tmp(1:Kb_cols))
       df_tmp(i) = dfun_i;
       !postprocess
       err_i = abs(dfun(1.0_Long*i)-dfun_i)
       if ( err_b_max < err_i ) err_b_max=err_i
    end do

    ! domain
    do i=Kb_rows+1,dim-Kb_rows
       ! first derivative
       if (i>Kb_rows+1) call push_back(f_tmp,fun(1.0_Long*(i+Kb_rows)),Kd_cols)
       dfun_i = dot_product(K_domain,f_tmp)
       ! second derivative
       if (i<Kd_cols) then
          df_tmp(i) = dfun_i;
       elseif (i==Kd_cols) then
          df_tmp(i) = dfun_i;
          do j=1,Kb_rows ! boundary
             d2fun_i = dot_product(K_boundary(j,:),df_tmp)
             errd2_i = abs(d2fun(1.0_Long*j)-d2fun_i)
             if ( errd2_b_max < errd2_i ) errd2_b_max=errd2_i
          end do
          d2fun_i = dot_product(K_domain,df_tmp)
          errd2_i = abs(d2fun(1.0_Long*(i-Kb_rows))-d2fun_i)
          if ( errd2_b_max < errd2_i ) errd2_b_max=errd2_i
       else ! domain
          call push_back(df_tmp,dfun_i,Kb_cols)
          d2fun_i = dot_product(K_domain,df_tmp)
          errd2_i = abs(d2fun(1.0_Long*(i-Kb_rows))-d2fun_i)
          if ( errd2_max < errd2_i ) errd2_max=errd2_i
       end if
       ! postprocess
       err_i = abs(dfun(1.0_Long*i)-dfun_i)
       if ( err_max < err_i ) err_max=err_i
    end do

    ! right boundary
    do i=dim-Kb_rows+1,dim
       ! first derivative
       dfun_i = dot_product(-K_boundary(dim-i+1,:),f_tmp(Kb_cols:1:-1))
       ! second derivative
       call push_back(df_tmp,dfun_i,Kb_cols)
       d2fun_i = dot_product(K_domain,df_tmp)
       errd2_i = abs(d2fun(1.0_Long*(i-Kb_rows))-d2fun_i)
       if ( errd2_max < errd2_i ) errd2_max=errd2_i
       ! postprocess
       err_i = abs(dfun(1.0_Long*i)-dfun_i)
       if ( err_b_max < err_i ) err_b_max=err_i
    end do
    ! second derivative
    do j=1,Kb_rows
       d2fun_i = dot_product(-K_boundary(Kb_rows-j+1,Kb_cols:1:-1),df_tmp)
       errd2_i = abs(d2fun(1.0_Long*(dim-Kb_rows+j))-d2fun_i)
       if ( errd2_b_max < errd2_i ) errd2_b_max=errd2_i
    end do
    
  end subroutine errors_kernel
  
end program diferencias
