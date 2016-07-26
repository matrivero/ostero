subroutine demo_subroutine(xx,yy)

real*8, dimension(2,2), intent(in) :: xx
real*8, dimension(2,2), intent(out) :: yy

print *,xx
yy = 5*xx

end subroutine demo_subroutine