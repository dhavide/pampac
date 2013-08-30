!_____________________________________________________________________________!
!____________________________________GMRES____________________________________!
subroutine gmres(NN,KK,iter,eps,xx,dxx,bb,flight,edr,conv,abort,k)
! Find an approximate solution to Ax=b using GMRES
  implicit none
! to get myrank
  include 'mpif.h'
  integer myrank,ierr
!for GMRES
  logical :: conv,abort
  integer,intent(in) :: NN,KK,iter
  real(8),intent(in) :: eps
  real(8),dimension(NN),intent(in) :: xx,dxx
  real(8),dimension(NN),intent(inout) :: bb
  real(8),dimension(NN) :: yy0,rr,xx1
  real(8),dimension(NN,KK+1) :: Q
  real(8),dimension(KK+1,KK+1) :: Omg
  real(8),dimension(KK,KK) :: H
  real(8),dimension(KK+1) :: yy
  real(8),intent(out) :: flight,edr
  real(8) :: dif,sqrtb2,r02,alpha,beta
  integer :: i,j,k,l,info

  call mpi_comm_rank(mpi_comm_world,myrank,ierr)
  if(ierr.ne.0) then
     write(*,*) 'mpi problem, exciting...'
     abort=.true.
     return
  end if

!test
!   open(33,file='krylovstatus') 


  conv=.false.
  abort=.false.
!  hess=0d0
!  df=bb

  sqrtb2=sqrt(sum(bb**2))
  yy0=0.d0
  do l=1,iter
     Q=0.d0; Omg=0.d0; Omg(1,1)=1.d0 ; H=0.d0;
     if(l==1)then
        rr=bb
     else
        call c0_to_Ac0(xx,xx1,yy0,rr,dxx,FLIGHT,EDR,1,abort)
        if(abort) return
        rr=bb-rr
     end if
     r02=sqrt(sum(rr**2))
     Q(:,1)=rr/r02
     do i=1,KK
        !___orthonormalization______________________________________________!
        !  Array 'yy' is temporarily used to store components of the i-th
        !  column of the Hessenberg matrix.
        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
        call c0_to_Ac0(xx,xx1,Q(:,i),rr,dxx,FLIGHT,EDR,1,abort)
        if(abort) return
        Q(:,i+1)=0.d0
        do j=1,i
           yy(j)=dot_product(rr,Q(:,j))

!           hess(j,i)=yy(j)

           Q(:,i+1)=Q(:,i+1)+yy(j)*Q(:,j)
        end do
        Q(:,i+1)=rr-Q(:,i+1)
        yy(i+1)=sqrt(sum(Q(:,i+1)**2))

!        hess(i+1,i)=yy(i+1)

        Q(:,i+1)=Q(:,i+1)/yy(i+1)

        !___convert_to_upper_triangular_matrix___
        yy(1:i)      =  matmul(yy(1:i),Omg(1:i,1:i))
        beta         =  1.d0/sqrt(yy(i)**2+yy(i+1)**2)
        alpha        =  yy(i  )*beta
        beta         =  yy(i+1)*beta
        yy(i)        =  yy(i)*alpha+yy(i+1)*beta
        H(1:i,i)     =  yy(1:i)
        yy(1:i)      =  Omg(1:i,i)
        Omg(1:i,i  ) =  alpha*yy(1:i)
        Omg(i+1,i  ) =  beta
        Omg(1:i,i+1) = -beta *yy(1:i)
        Omg(i+1,i+1) =  alpha
        yy=0.d0; yy(1)=r02
        yy(1:i+1)=matmul(yy(1:i+1),Omg(1:i+1,1:i+1))
        dif=abs(yy(i+1))/sqrtb2
        select case(myrank)
           case(0)
              open(30,file='0_status')
              write(30,*) l,i,dif
              close(30)
           case(1)
              open(30,file='1_status')
              write(30,*) l,i,dif
              close(30)
           case(2)
              open(30,file='2_status')
              write(30,*) l,i,dif
              close(30)
           case(3)
              open(30,file='3_status')
              write(30,*) l,i,dif
              close(30)
           end select
        
        if(dif<eps.or.i==KK)then
           k=i
           exit
        end if
     end do
!     write(33,*)'#',l,k,dif
     call dtrtrs('U','N','N',k,1,H(1:k,1:k),k,yy(1:k),k,info)
     yy0=yy0+matmul(Q(:,1:k),yy(1:k))
     if(dif<eps) then
        conv=.true.
        exit
     end if
  end do
  if(.not.(conv)) write(*,*) 'WARNING: no convergence in GMRES'
  bb=yy0
!  open(33,file='krylovstatus')
!  write(33,*) 'restarted ',l-1,'dimension=',k
!  write(33,90) bb
!  close(33) 
!90 format(e21.14)

  return
  end subroutine gmres
