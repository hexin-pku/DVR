module PES
implicit none
private
    public :: V_HO, V_QUAR
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
    V = 0.5_8 * x*x
    return
end subroutine V_QUAR

end module PES
