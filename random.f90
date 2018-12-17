!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- random module for sci-computations 
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module random
use const
implicit none
    logical, private :: has_initialed = .false.
    real(8), public  :: rand_u
    
contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- initial seed for random number's generation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine init_seed()
	integer :: n, ival(8), v(3), i
	integer, allocatable :: seed(:)
	call date_and_time(values=ival)
	v(1) = ival(8) + 2048*ival(7)
    v(2) = ival(6) + 64*ival(5)
    !-- value(4) isn't real(8)ly 'random'
    v(3) = ival(3) + 32*ival(2) + 32*8*ival(1)
 	call random_seed(size=n)
    allocate(seed(n))
    !-- Give the seed an implementation-dependent kick
	call random_seed()
	call random_seed(get=seed)
	do i=1, n
    	seed(i) = seed(i) + v(mod(i-1, 3) + 1)
  	enddo
  	call random_seed(put=seed)
  	deallocate(seed)
  	has_initialed = .true.
end subroutine init_seed

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- normal distribution
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine random_norm(rn_norm, iparms)
implicit none
    real(8), intent(inout) :: rn_norm
    real(8), dimension(2), intent(in), optional :: iparms
    real(8) :: mu, sigma, tmp1
    !-- check initialization
    if(has_initialed .eqv. .false.) call init_seed()
    !-- get parameters, default N(0,1)
    if(present(iparms)) then
        mu = iparms(1)
        sigma = iparms(2)
    else
        mu = 0.0
        sigma = 1.0
    endif
    !-- Box-Muller algorithm (cos side)
    call random_number(rand_u)
    tmp1 = dsqrt(-2*dlog(rand_u))
    call random_number(rand_u)
    tmp1 = tmp1*dcos(2*pi*rand_u)
    rn_norm = mu + tmp1*sigma
end subroutine random_norm

end module

