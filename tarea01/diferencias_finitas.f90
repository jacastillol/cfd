module diferencias_finitas
contains
  subroutine get_Kernel(Kernel, order, shift, dx)
    use parametros, only: Long
    implicit none
    integer, intent(in) :: order, shift
    real(Long), intent(in) :: dx
    real(Long), dimension(:), intent(out) :: Kernel
    integer :: j, delta, info
    integer, dimension(order) :: pivs
    real*8, dimension(order,order) :: A
    real*8, dimension(order) :: b

    A = 0.0D0
    b = 0.0D0; b(2) = 1.0D0
    do delta=shift,order+shift-1
       do j=0,order-1
          A(j+1,delta-shift+1) = real((real(delta,Long)**j/factorial(j))*dx,8)
       end do
    end do
    !    dgesv(N,NRHS,A,LDA,IPIV,b,LDB,info)
    call dgesv(order, 1, A, order, pivs, b, order, info)
    Kernel = real(b,Long)

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
    use parametros, only: Long
    implicit none
    integer, intent(in) :: dim, Kd_cols, Kb_cols, Kb_rows
    real(Long), dimension(:), intent(in) :: fun
    real(Long), dimension(:), intent(in) :: K_domain
    real(Long), dimension(:,:), intent(in) :: K_boundary
    real(Long), dimension(:), intent(out) :: dfun
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

  subroutine diff_matrix(fun,A,dim,K_domain,Kd_cols,K_boundary,Kb_cols,Kb_rows,dfun)
    use parametros, only: Long
    implicit none
    integer, intent(in) :: dim, Kd_cols, Kb_cols, Kb_rows
    real(Long), dimension(:), intent(in) :: fun
    real(Long), dimension(:), intent(in) :: K_domain
    real(Long), dimension(:,:), intent(in) :: K_boundary
    real(Long), dimension(:), intent(out) :: dfun
    integer :: i, j
    real(Long), dimension(:,:), intent(inout) :: A

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

    dfun = matmul(A,fun)
    !write(*,'(11(f9.5,x))') matmul(A,fun)
    !  DGEMM  ('N','N',  M,  N,  K, ALPHA, A,  M, B,  K,  BETA,   C,  M)
    !call dgemm('N','N',dim,  1,dim, 1.0_Long, A,dim, f,dim, 0.0_Long, df_,dim)
  end subroutine diff_matrix

  subroutine conv_difu_1D_uniform(phi,dim,K_stencil,phi0,phiL)
    use parametros, only: Long
    implicit none
    integer, intent(in) :: dim
    real(Long), dimension(3), intent(in) :: K_stencil
    real(Long), dimension(dim), intent(out) :: phi
    real(Long), dimension(dim-2) :: D
    real(Long), dimension(dim-3) :: DL, DU
    real(Long) :: phi0, phiL
    integer :: info

    ! boudaries
    phi(1) = phi0
    phi(2) = -K_stencil(1)*phi0
    phi(dim-1) = -K_stencil(3)*phiL
    phi(dim) = phiL

    ! domain
    DL = K_stencil(1)
    D = K_stencil(2)
    DU = K_stencil(3)

    !    dgtsv(N,NRHS,DL,D,DU,b,LDB,info)
    call dgtsv(dim-2, 1, DL, D, DU, phi(2:dim-1), dim-2, info)

  end subroutine conv_difu_1D_uniform

  subroutine conv_difu_1D_non_uniform(x,phi,dim,K_fun,phi0,phiL,rho,u,Gamma)
    use parametros, only: Long        
    implicit none
    integer, intent(in) :: dim
    real(Long), dimension(dim), intent(in) :: x
    real(Long), dimension(dim), intent(out) :: phi
    real(Long), dimension(3):: K_stencil
    real(Long), dimension(dim-2) :: D
    real(Long), dimension(dim-3) :: DL, DU
    real(Long), intent(in) :: rho, u, Gamma
    real(Long) :: phi0, phiL
    integer :: i, info

    interface K_fcns
       function K_fun(i,x,Gamma,rho,u)
         use parametros
         implicit none
         integer :: i
         real(Long), dimension(:) :: x
         real(Long) :: Gamma, rho, u
         real(Long), dimension(3) :: K_fun
       end function K_fun
    end interface K_fcns

    ! boudaries
    ! x=0
    i = 2;
    K_stencil = K_fun(i,x,Gamma,rho,u)
    phi(1) = phi0
    phi(2) = -K_stencil(1)*phi0
    D(i-1) = K_stencil(2)
    DU(i-1) = K_stencil(3)
    ! x=L
    i = dim-1
    K_stencil = K_fun(i,x,Gamma,rho,u)
    DL(i-2) = K_stencil(1)
    D(i-1) = K_stencil(2)
    phi(dim-1) = -K_stencil(3)*phiL
    phi(dim) = phiL
    
    ! domain
    do i=3,dim-2
       K_stencil = K_fun(i,x,Gamma,rho,u)
       DL(i-2) = K_stencil(1)
       D(i-1) = K_stencil(2)
       DU(i-1) = K_stencil(3)        
    end do
    
    !    dgtsv(N,NRHS,DL,D,DU,b,LDB,info)
    call dgtsv(dim-2, 1, DL, D, DU, phi(2:dim-1), dim-2, info)

  end subroutine conv_difu_1D_non_uniform

end module diferencias_finitas
