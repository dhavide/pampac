      INTEGER N,N1,N2,NE,NO,NT
      PARAMETER(N=42,N1=N/2,N2=(N-1)/2)
      PARAMETER(NE=(2*N1**3+3*N1**2-5*N1)/6,NO=(N2**3+3*N2**2+2*N2)/3)
      PARAMETER(NT=NE+NO)
!  for NN=6 (i.e. 64^3 resolution) set N=20,
!  for NN=7 (i.e 128^3 resolution) set N=42,
!  for NN=8 (i.e 256^3 resolution) set N=84, the largest even
!  number less then 2^8/3
!  likewise, for NN=9 set N=170
!  the integration code suppresses the largest wave number components
!  if they are odd!
