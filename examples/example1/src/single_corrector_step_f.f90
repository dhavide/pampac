subroutine single_corrector_step_f(x,dx,out_vars,conv)
! in: x, vector of length nt, point on a Poincare plane of
!       intersection close to a periodic point
!       dx, vector of length nt, approximate tangent to the solution curve
! out: if conv is TRUE, then
!     x is the new point on the Poincare plane,
!     dx is unchanged
!     the period and mean energy dissipation rate of the 
!     orbit segment are stored in preiod and edr
!     dim is the Krylov subspace dimension
!     errf is the new residue and errx the size of the Newton step
!      if conv is FALSE, ignore all output
! maxdim, maxit and gmrestol are the maximal Krylov subspace dimension, the maximal
! number of restarts and the tolerance of the GMRES routine - these are constant
! troughout the continuation and you might like to specify there in a different way
implicit none
include 'n_power.for'
include 'constants.h'
logical :: gmconv,abort,error
integer :: dim,maxdim,maxit,conv
real(8), dimension(nt-1) :: x,dx,x1,dx0,dx1,res
real(8) :: period,edr,errf,gmrestol
real(8),dimension(3) :: out_vars
common /gmres_pars/ gmrestol,maxdim,maxit

out_vars=0d0
conv=1
dx0=0d0

! compute the image of x and the residual of (P^k(x)-x), which 
! is the RHS in (DP^k-I) dx=x-P^k(x)

call C0_to_Ac0(x,x1,dx0,dx1,dx,period,edr,0,error)
if(error) then
   conv=0
   return
end if
res=x-x1

! approximate up to GMRESTOL the solution of
! (DP^k-I) dx=x-P^k(x) 
call gmres(nt-1,maxdim,maxit,gmrestol,x,dx,res,period,edr,gmconv,abort,dim)
if((.not.(gmconv)).or.(abort)) then
   conv=0
   return
end if

x=x+res

out_vars(1)=period!period of the orbit
out_vars(2)=edr!energy dissipation rate averaged along the orbit
out_vars(3)=float(dim)!Krylov subspace dimension

return
end subroutine single_corrector_step_f

