module PES
implicit none
private
    public :: V_HO, V_QUAR, V_DOUB
contains

subroutine V_HO(V,x)
    real(8), intent(out) :: V
    real(8), intent(in) :: x
    V = 0.5_8 * x*x
    return
end subroutine V_HO

subroutine V_QUAR(V,x)
    real(8), intent(out) :: V
    real(8), intent(in) :: x
    V = 0.25_8 * x**4
    return
end subroutine V_QUAR

subroutine V_DOUB(V,x)
    real(8), intent(out) :: V
    real(8), intent(in) :: x
    V = 0.5_8 * ( x**2 - 2.0 )**2
    return
end subroutine V_DOUB

end module PES
