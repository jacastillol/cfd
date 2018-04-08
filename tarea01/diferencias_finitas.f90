module diferencias_finitas
contains
  subroutine get_Kernel(Kernel, order, shift, dx)
    implicit none
    integer, intent(in) :: order, shift
    real*8, intent(in) :: dx
    real*8, dimension(:), intent(out) :: Kernel
    integer :: j, delta, info
    integer, dimension(order) :: pivs
    real*8, dimension(order,order) :: A
    real*8, dimension(order) :: b

    A = 0.0D0
    b = 0.0D0; b(2) = 1.0D0
    do delta=shift,order+shift-1
       do j=0,order-1
          A(j+1,delta-shift+1) = (real(delta)**j/factorial(j))*dx
       end do
    end do
    !    dgesv(N,NRHS,A,LDA,IPIV,b,LDB,info)
    call dgesv(order, 1, A, order, pivs, b, order, info)
    Kernel = b

  contains
    function factorial(n)
      implicit none
      integer :: n, i
      integer :: factorial
      factorial = 1
      do i=1,n
         factorial = factorial*i
      end do
    end function factorial
  end subroutine get_Kernel

  subroutine convolution(fun,dim,K_domain,Kd_cols,K_boundary,Kb_cols,Kb_rows,dfun)
    implicit none
    integer, intent(in) :: dim, Kd_cols, Kb_cols, Kb_rows
    real*8, dimension(:), intent(in) :: fun
    real*8, dimension(:), intent(in) :: K_domain
    real*8, dimension(:,:), intent(in) :: K_boundary
    real*8, dimension(:), intent(out) :: dfun
    integer :: i, j
    ! left boundary
    do i=1,Kb_rows
       do j=1,Kb_cols
          dfun(i) = dfun(i) + K_boundary(i,j)*fun(j)
       end do
    end do
    ! domain
    do i=Kb_rows+1,dim-Kb_rows
       do j=1,Kd_cols
          dfun(i) = dfun(i) + K_domain(j)*fun(i-((Kd_cols-1)/2)-1+j)
       end do
    end do
    ! right boundary
    do i=dim-Kb_rows+1,dim
       do j=1,Kb_cols
          dfun(i) = dfun(i) - K_boundary(dim-i+1,Kb_cols-j+1)*fun(dim-Kb_cols+j)
       end do
    end do
  end subroutine convolution

  subroutine error_fds(fun,dfun,dim,K_domain,Kd_cols,K_boundary,Kb_cols,Kb_rows)
    implicit none
    integer, intent(in) :: dim, Kd_cols, Kb_cols, Kb_rows
    real*8, dimension(:), intent(in) :: K_domain
    real*8, dimension(:,:), intent(in) :: K_boundary
    real*8 :: fun_i, dfun_i
    real*8 :: err_i, err_max = -1, err_b_max = -1
    integer :: i, j

    interface AFunc
       function fun (x)
         real*8 :: fun
         real*8 :: x
       end function fun
       function dfun (x)
         real*8 :: dfun
         real*8 :: x
       end function dfun
    end interface AFunc

    ! left boundary
    do i=1,Kb_rows
       dfun_i = 0
       do j=1,Kb_cols
          dfun_i = dfun_i + K_boundary(i,j)*fun(dx*j)
       end do
       err_i = abs(dfun(dx*i)-dfun_i)
       if ( err_b_max < err_i ) err_b_max=err_i
    end do
    ! domain
    do i=Kb_rows+1,dim-Kb_rows
       dfun_i = 0
       do j=1,Kd_cols
          dfun_i = dfun_i + K_domain(j)*fun(dx*(i-((Kd_cols-1)/2)-1+j))
       end do
       err_i = abs(dfun(dx*i)-dfun_i)
       if ( err_max < err_i ) err_max=err_i
    end do
    ! right boundary
    do i=dim-Kb_rows+1,dim
       dfun_i = 0
       do j=1,Kb_cols
          dfun_i = dfun_i - K_boundary(dim-i+1,Kb_cols-j+1)*fun(dx*(dim-Kb_cols+j))
       end do
       err_i = abs(dfun(dx*i)-dfun_i)
       if ( err_b_max < err_i ) err_b_max=err_i
    end do
  end subroutine error_fds

  subroutine diff_matrix(fun,A,dim,K_domain,Kd_cols,K_boundary,Kb_cols,Kb_rows,dfun)
    implicit none
    integer, intent(in) :: dim, Kd_cols, Kb_cols, Kb_rows
    real*8, dimension(:), intent(in) :: fun
    real*8, dimension(:), intent(in) :: K_domain
    real*8, dimension(:,:), intent(in) :: K_boundary
    real*8, dimension(:), intent(out) :: dfun
    integer :: i, j
    real*8, dimension(:,:), intent(inout) :: A

    ! left boundary
    do i=1,Kb_rows
       do j=1,Kb_cols
          A(i,j) = K_boundary(i,j)
       end do
    end do
    ! domain
    do i=Kb_rows+1,dim-Kb_rows
       do j=1,Kd_cols
          A(i,i-((Kd_cols-1)/2)-1+j) = K_domain(j)
       end do
    end do
    ! right boundary
    do i=dim-Kb_rows+1,dim
       do j=1,Kb_cols
          A(i,dim-Kb_cols+j) = - K_boundary( dim-i+1,Kb_cols-j+1 )
       end do
    end do

    write(*,'(11(f6.2,x))') A
    write(*,'(11(f9.5,x))') matmul(A,fun)
    !  DGEMM  ('N','N',  M,  N,  K, ALPHA, A,  M, B,  K,  BETA,   C,  M)
    !call dgemm('N','N',dim,  1,dim, 1.0D0, A,dim, f,dim, 0.0D0, df_,dim)
    ! write(*,'(11(f6.2,x))') y
  end subroutine diff_matrix

end module diferencias_finitas
