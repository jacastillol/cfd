program burger2d_sol
  implicit none
  !integer, parameter :: sp = selected_real_kind(6, 37)
  !integer, parameter :: dp = selected_real_kind(15, 307)
  !integer, parameter :: qp = selected_real_kind(33, 4931)
  integer, parameter :: Long = selected_real_kind(15, 307)
  integer, parameter :: dim=21
  character(20) :: filename, x1, fmt='(I4.4)'
  real(Long) :: Uteo=0.0_Long, Vteo=0.0_Long
  real(Long), dimension(dim,dim) :: U=0.0_Long, V=0.0_Long
  real(Long), dimension(dim,dim) :: Unext=0.0_Long, Vnext=0.0_Long
  real(Long), parameter :: dt=0.0001, h=0.05
  real(Long), parameter :: Re = 80
  real(Long), parameter :: Kx=(dt/Re)/(h**2.0_Long), &
       Ky=(dt/Re)/(h**2.0_Long), Kx2=dt/(2*h), Ky2=dt/(2*h)
  integer :: i, j, n, nmax = 5001
  real(Long) :: t, x, y, err_max=-1, err_ij=-1

  ! condicion inicial dominio
  do i=2,dim-1
     do j=2,dim-1
        x = h*(i-1); y = h*(j-1)
        U(i,j) = 0.75_Long - 1.0_Long/( 4*( 1+exp( Re*(-t-4*x+4*y)/32 ) ) )
        V(i,j) = 0.75_Long + 1.0_Long/( 4*( 1+exp( Re*(-t-4*x+4*y)/32 ) ) )
     end do
  end do

  open(unit=12, file="data/error.dat")
  ! evolucion
  do n=1,nmax
     t=dt*(n-1)
     ! condicion inicial
     do i=1,dim
        do j=1,dim,dim-1
           x = h*(i-1); y = h*(j-1)
           U(i,j) = 0.75_Long - 1.0_Long/( 4*( 1+exp( Re*(-t-4*x+4*y)/32 ) ) )
           V(i,j) = 0.75_Long + 1.0_Long/( 4*( 1+exp( Re*(-t-4*x+4*y)/32 ) ) )
           y = h*(i-1); x = h*(j-1)
           U(j,i) = 0.75_Long - 1.0_Long/( 4*( 1+exp( Re*(-t-4*x+4*y)/32 ) ) )
           V(j,i) = 0.75_Long + 1.0_Long/( 4*( 1+exp( Re*(-t-4*x+4*y)/32 ) ) )
        end do
     end do

     ! posprocess
     if (mod(n,100)==1) then
        write(x1,fmt) n-1
        filename = 'data/solU_t'//trim(x1)//'.dat'
        open(unit=11, file=filename)
        write(11,'(21(1x,f8.5))') ((U(i,j), j=1,dim), i=1,dim)
        close(11)
        filename = 'data/solV_t'//trim(x1)//'.dat'
        open(unit=11, file=filename)
        write(11,'(21(1x,f8.5))') ((V(i,j), j=1,dim), i=1,dim)
        close(11)
     end if
  
     do i=2,dim-1
        do j=2,dim-1
           x = h*(i-1); y = h*(j-1)
           Uteo = 0.75_Long - 1.0_Long/( 4*( 1+exp( Re*(-t-4*x+4*y)/32 ) ) )
           Vteo = 0.75_Long + 1.0_Long/( 4*( 1+exp( Re*(-t-4*x+4*y)/32 ) ) )
           ! velocidad en x
           Unext(i,j) = (Kx-Kx2*U(i,j))*U(i+1,j) &
                +(Ky-Ky2*V(i,j))*U(i,j+1) &
                +(1.0_Long-2.0_Long*Kx-2.0_Long*Ky)*U(i,j) &
                +(Kx2*U(i,j)+Kx)*U(i-1,j) &
                +(Ky2*V(i,j)+Ky)*U(i,j-1)
           ! velocidad en y
           Vnext(i,j) = (Kx-Kx2*V(i,j))*V(i+1,j) &
                +(Ky-Ky2*U(i,j))*V(i,j+1) &
                +(1.0_Long-2.0_Long*Kx-2.0_Long*Ky)*V(i,j) &
                +(Kx2*V(i,j)+Kx)*V(i-1,j) &
                +(Ky2*U(i,j)+Ky)*V(i,j-1)
           err_ij = max(abs(Vteo-Vnext(i,j)),abs(Uteo-Unext(i,j)))
           if (err_ij>err_max) err_max=err_ij
        end do
     end do
     U(2:dim-1,2:dim-1) = Unext(2:dim-1,2:dim-1)
     V(2:dim-1,2:dim-1) = Vnext(2:dim-1,2:dim-1)
     write(12,*) t, n, err_max
  end do
  close(12)
  
end program burger2d_sol
