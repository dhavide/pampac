subroutine calculate_residual_f(x,res,conv)
! in: x, vector of length nt, point on a Poincare plane of
!       intersection close to a periodic point
!       dx, vector of length nt, approximate tangent to the solution curve
! out: if conv is TRUE, res is the residual ||(P^k(x)-x)||
implicit none
include 'constants.h'
logical :: error
integer :: conv
real(8), dimension(nt-1) :: x,dx,dx0,dx1,x1,res
real(8) :: period,edr
conv=1
dx0=0d0
dx=0d0

call C0_to_Ac0(x,x1,dx0,dx1,dx,period,edr,0,error)
if(error) then
   conv=0
   return
end if
res=x-x1

return
end subroutine calculate_residual_f
