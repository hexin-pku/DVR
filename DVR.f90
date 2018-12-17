module DVR
use const
use linalgebra, only : lin_ev
use PES, V=>V_HO
implicit none
private
    real(8), parameter :: M = 1.0_8
    integer :: N
    real(8) :: a, b
    real(8), allocatable :: x(:), H(:,:)
    
    public :: finite_DVR, infinite_DVR, half_infinite_DVR, del_DVR

contains

subroutine finite_DVR(input_a, input_b, input_N, Es, Vs)
    integer, intent(in) :: input_N
    real(8), intent(in) :: input_a, input_b
    real(8), intent(out) :: Es(input_N-1), Vs(input_N-1, input_N-1)
    real(8) :: dx, coeff, tmp_V
    integer :: i,j
    
    a = input_a
    b = input_b
    N = input_N
    coeff = 0.5/M /(b-a)**2 * 0.5*pi2
    
    if(.not. allocated(x) ) allocate(x(N-1), H(N-1,N-1))
    
    dx = (b-a) / real(N)
    do i=1,N-1
        x(i) = a + i*dx
    enddo
    
    do i=1,N-1
        H(i,i) = coeff * ( (2*N*N+1)/real(3) - 1.0 / sin(pi*i/real(N))**2 )
        do j=i+1,N-1
            H(i,j) = coeff * (-1)**(i-j) * ( 1.0 / sin(pi*(i-j)/real(2*N))**2 &
                             - 1.0 / sin(pi*(i+j)/real(2*N))**2 )
            H(j,i) = H(i,j)
        enddo
    enddo
    
    do i=1,N-1
        call V(tmp_V, x(i))
        H(i,i) = H(i,i) + tmp_V
    enddo 
    
    call lin_ev(Es, Vs, H)    
    return
end subroutine finite_DVR 

subroutine infinite_DVR(input_L, input_N, Es, Vs)
    integer, intent(in) :: input_N
    real(8), intent(in) :: input_L
    real(8), intent(out) :: Es(input_N-1), Vs(input_N-1, input_N-1)
    real(8) :: dx, coeff, tmp_V
    integer :: i,j
    
    if(input_L <= 0 ) stop 'length error'
    if(mod(N,2) .eq. 1) stop 'input_N should be even'
    print *, 'make sure L is enough large and N is even'
    
    a = -input_L/2
    b = input_L/2
    N = input_N
    dx = (b-a) / real(N)
    coeff = 0.5/M /dx**2
    
    if(.not. allocated(x) ) allocate(x(N-1), H(N-1,N-1))
    
    do i=1,N-1
        x(i) = a + i*dx
    enddo
    
    do i=1,N-1
        H(i,i) = coeff * pi2 / real(3)
        do j=i+1,N-1
            H(i,j) = coeff * (-1)**(i-j) * 2/real(i-j)**2
            H(j,i) = H(i,j)
        enddo
    enddo
    
    do i=1,N-1
        call V(tmp_V, x(i))
        H(i,i) = H(i,i) + tmp_V
    enddo 
    
    call lin_ev(Es, Vs, H)    
    return
end subroutine infinite_DVR 

subroutine half_infinite_DVR(input_L, input_N, Es, Vs)
    integer, intent(in) :: input_N
    real(8), intent(in) :: input_L
    real(8), intent(out) :: Es(input_N-1), Vs(input_N-1, input_N-1)
    real(8) :: dx, coeff, tmp_V
    integer :: i,j
    
    if(input_L <= 0 ) stop 'length error'
    
    a = 0
    b = input_L
    N = input_N
    dx = (b-a) / real(N)
    coeff = 0.5/M /dx**2
    
    if(.not. allocated(x) ) allocate(x(N-1), H(N-1,N-1))
    
    do i=1,N-1
        x(i) = a + i*dx
    enddo
    
    do i=1,N-1
        H(i,i) = coeff * ( pi2 / real(3) - 0.5*i*i )
        do j=i+1,N-1
            H(i,j) = coeff * (-1)**(i-j) * ( 2/real(i-j)**2 - 2/real(i+j)**2 )
            H(j,i) = H(i,j)
        enddo
    enddo
    
    do i=1,N-1
        call V(tmp_V, x(i))
        H(i,i) = H(i,i) + tmp_V
    enddo 
    
    call lin_ev(Es, Vs, H)    
    return
end subroutine half_infinite_DVR 

subroutine del_DVR()
    deallocate(x,H)
end subroutine del_DVR

end module DVR



program numerical_wigner
use DVR
implicit none
    integer :: Nx
    integer :: Np
    real(8) :: L, xmin, xmax, pmin, pmax, x, p, dx, tmp, beta
    real(8), allocatable :: Es(:), Vs(:,:), rho(:,:)
    integer :: i,j,k,m,ixmd,ixpd
    
    beta = 8
    Nx = 200
    NP = 200
    xmin = -10
    xmax = 10
    pmin = -10
    pmax = 10

    dx = (xmax-xmin) / real(Nx)
    
    allocate(Es(Nx-1), Vs(Nx-1,Nx-1),rho(Nx-1, Np-1))
    
    call finite_DVR(xmin, xmax, Nx, Es, Vs)
    print *, Es(1:10)
    
    do i=1,Nx-1
        ! x = xmin + (xmax-xmin) * i / real(Nx)
        do j=1,Np-1
            p = pmin + (pmax-pmin) * j /real(Np)
            rho(i,j) = 0
            do k=1, 10
                tmp = 0
                do m=0, min(i,Nx-i)
                    tmp = tmp + Vs(k,i+m)*Vs(k,i-m)*cos(2*m*dx*p)*2*dx
                enddo
                rho(i,j) = rho(i,j) + exp(-beta*Es(k)) * tmp
            enddo
        enddo
    enddo

    open(unit=222, file='dist.dat',status='replace')
    do i=1,Nx-1
    write(222,*) rho(i,:)
    enddo 
    close(unit=222)  


    deallocate(Es,Vs)
    call del_DVR()
end program numerical_wigner














