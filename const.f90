!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!-- base module for sci-computations 
!-- (Copyright) Xin He <1500011805@pku.edu.cn>)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module const
implicit none
    !-- precision of float number and length of strings
    integer, parameter :: sp=kind(0.0)
    integer, parameter :: dp=kind(0.d0)
    integer, parameter :: len0=4
    integer, parameter :: len1=20
    integer, parameter :: len2=100
    integer, parameter :: len3=500

    !-- mathematical constant numbers #XXX fixed bugs, the kind must noted at the parameters!
    real(8), parameter :: pi = 3.14159265358979323846_8
    real(8), parameter :: pi2 = pi**2
    real(8), parameter :: twopi = 2._8*pi
    real(8), parameter :: sqrtpi = dsqrt(pi)
    real(8), parameter :: eu = 2.71828182845904523536_8
    real(8), parameter :: ec = 0.57721566490153286060_8

    !-- physics constant
    real(8), parameter :: hb = 1.0545718001E-34_8
    real(8), parameter :: kb = 1.3806490351E-23_8
    
    !-- atom unit
    !-- quantity(a.u.) == au2si_quantity * quantity(si)
    real(8), parameter :: au2si_mass   = 9.10938291E-31_8     ! (kg)
    real(8), parameter :: au2si_charge = 1.602176565E-19_8    ! (C)
    real(8), parameter :: au2si_electric_const = 8.9875517873681E+9_8 ! ( kg*m^3*s^-2*C^-2 )
    real(8), parameter :: au2si_action = 1.054571726E-34_8    ! (J*s)
    real(8), parameter :: au2si_length = 5.2917721092E-11_8   ! (m)
    real(8), parameter :: au2si_momentum = 1.99285188224E-24_8 ! (kg*m/s)
    real(8), parameter :: au2si_energy = 4.35974417E-18_8     ! (J)
    real(8), parameter :: au2si_time = 2.418884326505E-17_8   ! (s)
    real(8), parameter :: au2si_temperatue = 3.157746455E+5_8 ! (K)
    real(8), parameter :: au2si_beta = 3.166815367E-6_8       ! (1/K)
    
    !-- others:
    !----- a_in_b : 1(a) == a_in_b (b) 
    !-- energy transfrom
    real(8), parameter :: au_in_wn = 219474.6313702_8
    
    !-- length transfrom
    real(8), parameter :: au_in_ai = 5.2917721092E-1_8
    
    !-- time transfrom
    real(8), parameter :: au_in_fs = 2.418884326505E-2_8
    real(8), parameter :: au_in_ps = 2.418884326505E-5_8
    
end module



