      subroutine C0_to_Ac0(x0,x1,dx0,dx1,tan,flight,edr,switch,error)
c  if SWITCH=0: computes the iterated Poincare map on x0
c  if SWITCH.ne.0:
c  computes the action of the linearized map (DP^k-I) on a perturbation
c  plus the residue inp(tan,dx0)
c  on input:
c  x0 - vector of dimension (nt-1) whith vorticity field and vorticity
c  dx0 - vector of dimension (nt-1) whith corresponding perturbations
c  k - number of iterations of the Poincare map
c  on output:
c  if flag=0 then
c    x1 - image of x0 under the k-th iterate of the 
c        Poincare map
c    dx1 - image of dx0 under the k-th
c         iterate of the linearized Poincare map minus unit (DP^k-I)
c    flight - time of flight for the k-th iterate (if x0=x1
c             this is the period of a cycle)
c  if flag>0 an error occured and x1, dx1 and flight SHOULD BE IGNORED
c    flag=1 - intersection with the Poincare plane not found
c    flag=2 - no convergence in the Newton iterations to find the 
c             intersection point

      implicit none

c  constants.h contains the dimension nt of x and dx
      include 'constants.h'
c  n_power.for contains the truncation level nn
      include 'n_power.for'
c  parameters.h contains the time step size, discrete period, poincare
c  interesection data and gmres parameters - declarations only -
      include 'parameters.h'

      logical error
      integer flag,ns,nsc,i,switch,ichoose
      double precision x0,x1,dx0,dx1,flight,nu,tan
      double precision x0l,x1l,dx0l,dx1l,a,e,dnu,varnu,edr

c  nsc is the dimension of the variables in the integration routine
      parameter(ns=2**nn/6,nsc=(ns+1)**3)

      dimension x0(nt-1),dx0(nt-1),x1(nt-1),dx1(nt-1),a(nsc),e(nsc)
      dimension x0l(nt+1),dx0l(nt+1),x1l(nt+1),dx1l(nt+1),tan(nt-1)

c  the vorticity field A and perturbation field E are in COMMON blocks
      COMMON /VOR1/A
      COMMON /VOR3/E
c  and so are the viscosity and the step size for integration
!      COMMON /CONST/ANU,DT ! now included in 'parameters.h'
c  common block /period/ holds the discrete period
!      common /period/ k ! now included in 'parameters.h'
c  common block /varnu/ holds the variation of nu
      common /varnu/ varnu
c  common block /choose/ holds SWITCH, iff SWITCH=0 no perturbation
c  is integrated
      common /choose/ ichoose
      include 'set_parameters.f'! time step, poincare, forcing and gmres
! data are set here
      ichoose=switch

      error=.false.

c  'lift' input vectors - i.e. insert the right values for the forcing
c  component and the Poincare section
      call lift(x0,x0l,0)
      call lift(dx0,dx0l,1)

c  translate the vector x and perturbation vector dx into vorticity
c  field A and perturbation field E
      call ytoa(x0l(1),a)
      call ytoa(dx0l(1),e)
      anu=x0l(nt+1)
      varnu=dx0l(nt+1)

c  call the integration routine
      call compute(disc_per,flight,edr,flag)
      if(flag.ne.0) then
         write(*,*) 'Warning: error flag ',flag
         error=.true.
         return
      endif

c  translate the vorticity field A and perturbation field E into vector 
c  x and perturbation vector dx
      call atoy(a,x1l(1))
      call atoy(e,dx1l(1))
      x1l(nt+1)=anu
      dx1l(nt+1)=varnu

c  'strip' the resulting (perturbation) vorticity fields of the fixed
c  components
      call strip(x1l,x1,0)
      call strip(dx1l,dx1,1)
c  IMPORTANT subtract the unit to compute (DP^k-I)dx0
c  consider where to put this in the code..
      do i=1,nt-2
         dx1(i)=dx1(i)-dx0(i)
      enddo
c  the extra equation:
      dx1(nt-1)=0d0
      do i=1,nt-1
         dx1(nt-1)=dx1(nt-1)+tan(i)*dx0(i)
      enddo

      return
      end

      subroutine lift(x,y,switch)
c  x: vector of (nt-1) unknowns, i.e. non-fixed Fourier components in the
c  Poincare intersection plane INPUT plus the viccosity
c  y: vector of (nt) Fourier components, including the forced component
c  and the component fixed on the Poincare plane of intersection OUTPUT
c  switch: if switch=0, x is assumed to be a vorticity vector
c          if switch<>0, x is assumed to be a perturbation vorticity
c          vector INPUT

      implicit none
      include 'constants.h'
      include 'n_power.for'

      integer i,j,switch,ifix,nnwt,nfix,npsfix
      double precision x,y,secxtol,secftol,xsec,pmax
      
      dimension x(nt-1),y(nt+1)

      common /fixed/ ifix
      common /section/ secxtol,secftol,xsec,pmax,nnwt,nfix,npsfix

      j=0
      do i=1,nt+1
         if(i.eq.ifix) then
            j=j+1
            if(switch.eq.0) then
               y(i)=-0.375d0
            else
               y(i)=0d0
            end if
         else
            if(i.eq.nfix) then
               j=j+1
               if(switch.eq.0) then
                  y(i)=xsec
               else
                  y(i)=0d0
               end if
            else
               y(i)=x(i-j)
            end if
         end if
 1    enddo

      return
      end

      subroutine strip(x,y,switch)
c  strips the (perturbation) vorticity vector of the component in the Poincare
c  intersection plane and the fixed component
c  x: vector of (nt) Fourier components, including the forced component
c  and the component fixed on the Poincare plane of intersection INPUT
c  y: vector of (nt-1) unknowns, i.e. non-fixed Fourier components in the
c  Poincare intersection plane plus viscosity INPUT
c  switch: if switch=0, x is assumed to be a vorticity vector
c          if switch<>0, x is assumed to be a perturbation vorticity
c          vector INPUT
      implicit none
      include 'constants.h'
      include 'n_power.for'
      integer i,j,switch,ifix,nnwt,nfix,npsfix
      double precision small,x,y,secxtol,secftol,xsec,pmax

      dimension x(nt+1),y(nt-1)
      parameter(small=1d-9)

      common /fixed/ ifix 
      common /section/ secxtol,secftol,xsec,pmax,nnwt,nfix,npsfix

      j=0
      do i=1,nt+1
         if(i.eq.ifix) then
            j=j+1
         else
            if(i.eq.nfix) then
               j=j+1
               if(switch.eq.0) then
                  if(dabs(x(nfix)-xsec).gt.small) then
                     write(*,*) 'WARNING: phase point does not seem to
     c           lie in intersection plane. Badness=',dabs(x(nfix)-xsec)
                  end if
               end if
               if(switch.ne.0) then
                  if(dabs(x(nfix)).gt.small) then
                     write(*,*) 'WARNING: nonzero component projected
     c                out of perturbation field. Badness=',dabs(x(nfix))
                  end if
               end if
            else
               y(i-j)=x(i)
            end if
         end if
 1    enddo

      return
      end

      subroutine compute(k,flight,edr,iflag)
      
*PROCESS DOUBLE
C
C ****************************************************************
C *                                                              *
C *           TIME-DEVELOPMENT OF A DIFFERENCE FIELD OF          *
C *                   3-D HIGH-SYMMETRIC FLOWS                   *
C *                           (LINEAR)                           *
C *                                                              *
C ****************************************************************
C
C
C            MADE BY SHIGEO KIDA ON 2 MAY 1987
C
C     THROUGH MODIFICATION TO NS3DIFS1: THE MAIN MODIFICATION IS
C     MADE BY REDUCING THE SIZE IN MEMO2 FROM NLCPH TO NLC AND
C     ADDING MEMO3.
C
C
C     CALCULATION IS MADE IN SPECTRAL REPRESENTATIONS.
C
C     TIME-STEPPING IS DONE BY FORTH ORDER RUNGE-KUTTA-GILL FORMULA.
C
C     NON-LINEAR TERMS ARE CALCULATED BY PSEUDO-SPECTARAL METHOD.
C
C     FOURIER TRANSFORMATIONS OF DIMENSION  N  ARE USED.  IN ORDER TO
C  ELIMINATE THE ALIASING ERROR FOURIER COMPONENTS ARE TRUNCATED AT
C  WAVENUMBER  N/3.  NUMBER OF FOURIER MODES RETAINED IS THEREFORE
C  (2N/3)**3.
C
C
C  ON EXECUTING THIS PROGRAM, WE MUST CHECK THE FOLLOWINGS:
C
C  (1)  NN  IN ALL THE PARAMETER STATEMENTS (IN MAIN AND SUB PROGRAMS),
C       WHICH MUST NOT BE LESS THAN  5.
C
C  (2)  PARAMETERS (ANU, DT, RUN, ISTEP1, ISTEP2, IOUT, IEN, IHEL  AND
C                   ISTART, IFORCE ETC.) IN CSETALL".
C
C  (3)  INITIAL CONDITION IN  "INIT".
C
C  (4) FORCING TERMS IN "STEP".
C
C  (5) TYPE OF PERTURBATIONS IN "PERTRB".
C
C **************************************************************
C *                                                            *
C *   RVF: REFERENCE VORTICITY FILED                           *
C *   DVF: DIFFERENCE VORTICITY FILED                          *
C *                                                            *
C *   FILES                     QUANTITY                       *
C *    11     RVF (TO BE READ IF ISTART.NE.0),                 *
C *    12     DVF (TO BE READ IF ISTART.NE.0),                 *
C *    13     RVF (TO BE WRITTEN EVERY IOUT STEPS),            *
C *    14     DVF (TO BE WRITTEN EVERY IOUT STEPS),            *
C *    15     ENERGY AND ENSTROPHY OF RVF (EVERY IEN STEPS),   *
C *    16     ENERGY AND ENSTROPHY OF DVF (EVERY IEN STEPS),   *
C *    18     HELICITY OF RVF (EVERY IHEL STEPS),              *
C *    19     HELICITY OF DVF (EVERY IHEL STEPS),              *
C *    21     SEVERAL FOURIER COMPONENTS OF RVF                *
C *                   (EVERY TIME STEP IF IFORCE=1 OR IFIX=1)  *
C *    22     SEVERAL FOURIER COMPONENTS OF DVF                *
C *                   (EVERY TIME STEP IF IFORCE=1 OR IFIX=1)  *
c   I disabled the writing of files 21-2 
C *    23     FIXED FOURIER COMPONENTS OF RVF                  *
C *                          (EVERY TIME STEP IF IFIX=1)       *
C *    24     FIXED FOURIER COMPONENTS OF DVF                  *
C *                          (EVERY TIME STEP IF IFIX=1)       *
C *                                                            *
C **************************************************************
C
C ****************************************************************
C *                                                              *
C *     CONTENTS IN MEMO1,MEMO2,MEMO3 AND MEMO4 (IN MARCH)       *
C *                                                              *
C *   STATEMENTS IN MARCH      MEMO1    MEMO2    MEMO3    MEMO4  *
C *    CALL WU                            UF                     *
C *    CALL UFR                           U                      *
C *    CALL COPY                 U        U                      *
C *    CALL SQUARE               U        S                      *
C *    CALL SRF, CALL NLS        U        SF                     *
C *    CALL COPY                 U        SF               SF    *
C *    CALL COPY                 U        U                SF    *
C *    CALL UT                   U        T                SF    *
C *    CALL TRF, CALL NLT        U        TF               SF    *
C *    CALL SUM                  U        TF              SF+TF  *
C *                                                              *
C *    CALL WU                   U        UF'             SF+TF  *
C *    CALL UFR                  U        U'              SF+TF  *
C *    CALL UT2                  U        U'       T'     SF+TF  *
C *    CALL US2                  S'       U'       T'     SF+TF  *
C *    CALL COPY                 S'       T'              SF+TF  *
C *    CALL TRF, CALL NLT        S'       TF'             SF+TF  *
C *    CALL COPY                 S'       TF'      TF'    SF+TF  *
C *    CALL COPY                 S'       S'       TF'    SF+TF  *
C *    CALL COPY                 TF'      S'       TF'    SF+TF  *
C *    CALL SRF, CALL NLS        TF'      SF'             SF+TF  *
C *    CALL SUM               SF'+TF'     SF'             SF+TF  *
C *                                                              *
C ****************************************************************
C
c
      implicit double precision (a-h,o-z)
c
      include 'n_power.for'
      logical up
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
      PARAMETER(NVW=ND8*(NN-1))
C
C  *************************************
C  *  A: REFERENCE OF VORTICITY FIELD  *
C  *  B: WORKING AREA                  *
C  *  E: DIFFERENCE OF VORTICITY FIELD *
C  *  H: WORKING AREA                  *
C  *************************************
C
      COMMON /VOR1/A(NSC)
      COMMON /VOR2/B(NSC)
      COMMON /VOR3/E(NSC)
      COMMON /VOR4/G(NSC)
      COMMON /MEMO1/C(NLC)
      COMMON /MEMO2/D(NLC)
      COMMON /MEMO3/F(NLC)
      COMMON /MEMO4/H(NSC)
      COMMON /WORK1/W1(0:NL,0:NL)
      COMMON /WORK2/W2(0:NL,0:NL)
      COMMON /WORK3/W3(0:NL)
      COMMON /WORK4/W4(0:NL,0:NL)
      COMMON /WORK5/W5(0:NL,0:NL)
      COMMON /WORK6/W6(0:NL,0:NL)
c      COMMON /WORK7/VW(NVW),IVW(NVW)
      COMMON /WORK7/ip(0:nd8)
      common /work8/wwww(0:nd8*5/2-1)
      COMMON /CONST/ANU,DT
      COMMON /PARA/IORTH,IHEL,IOUT,IEN,ISTEP2
      COMMON /TIME/TIME,ISTEP
      COMMON /FORCE/IFORCE,IFIX
c  the common block /section/ holds some constants related to the
c  Poincare section
      common /section/ secxtol,secftol,xsec,pmax,nnwt,nfix,npsfix

      dimension abc(2*nd8)
C
c      CALL DVCFT1(A,B,ND8,-1,0,VW,IVW,ICON)
      call two2one(H,B,abc,nd8)
      ip(0)=0
      call cdft(2*nd8,1,abc,ip,wwww)
      call one2two(abc,H,B,nd8)
ccc these calls only initialise the FFT routine
ccc calls two2one and one2two SHOULD NOT TAKE A FOR INPUT!
ccc H and B are dummy arguments which are zeroed below

      CALL ZERO(W1,(NL+1)*(NL+1))
      CALL ZERO(C,NLC)
      CALL ZERO(D,NLC)
      CALL ZERO(F,NLC)
c  E, the perturbation field, is now an input argument!
c      CALL ZERO(E,NSC)
      CALL ZERO(G,NSC)
      CALL ZERO(H,NSC)
C
c  set a few constants: force=0 for no external forcing, ifix=1 for fixed
c  Fourrier components
      IFORCE=0
      IFIX=1
c  every iorth steps the divergence-free property is imposed
      IORTH=300
c  istep2 is the maximal number of steps between intesection points
      ISTEP2=int(pmax/dt)
c  every iout steps the fields are written to files fort.13 and fort.14
c  iout>istep2 disables this feature, set to small value to generate
c  output for testing purposes
      iout=istep2+1

c  every ien steps the energy and enstrophy of the fields is written to 
c  fort.15 and fort.16
      ien=istep2+1
c  every ihel steps the helicity of the fields is written to 
c  fort.18 and fort.19
      ihel=istep2+1
c  set error flag to 0
      iflag=0
C
      CALL FACT
      CALL BASE
      CALL WNUM
C
      TIME=0
      ISTEP=0
      I=0
      edr=0d0

c  check in which direction we cross the section and initalise the counter

ccccccccccccccccccccccccccccccc
c      write(*,*) time,a(npsfix),e(npsfix)
ccccccccccccccccccccccccccccccccc
      TIME=TIME+DT
      ISTEP=ISTEP+1
      I=I+1
      CALL MARCH(A,E,I-1)
      call moment(a,ens,15)
      edr=edr+ens*dt

ccccccccccccccccccccccccccccccc
c      write(*,*) time,a(npsfix),e(npsfix)
ccccccccccccccccccccccccccccccccc


      up=.false.
      if(a(npsfix).gt.xsec) up=.true.
      isec=0
      amem=a(npsfix)
c  if up=.true. we intersect "upward", otherwise "downward"
c  amem stores the previous value of the component fixed in the plane
c  isec is a counter for the number of iterations of the Poincare map

   30 TIME=TIME+DT
      ISTEP=ISTEP+1
      I=I+1

      CALL MARCH(A,E,I-1)

      IF(MOD(ISTEP,IORTH).EQ.0) CALL ORTH(A,B,NS)
      IF(MOD(ISTEP,IORTH).EQ.0) CALL ORTH(E,B,NS)
      IF(MOD(I,IEN).EQ.0) CALL MOMENT(A,ens,15)
      IF(MOD(I,IEN).EQ.0) CALL MOMENT(E,ens,16)
      IF(MOD(I,IOUT).EQ.0) CALL STORE(A,13)
      IF(MOD(I,IOUT).EQ.0) CALL STORE(E,14)
C
      if((up.and.(a(npsfix).gt.xsec).and.(amem.lt.xsec))
     c   .or.(.not.(up).and.(a(npsfix).lt.xsec).and.(amem.gt.xsec)))then

c  we have computed an iteration of the Poincare map
            isec=isec+1
            i=0
            istep=0
            if(isec.eq.k) goto 77
      endif
      amem=a(npsfix) 
      call moment(a,ens,15)
      edr=edr+ens*dt

      IF(ISTEP.LT.ISTEP2) GO TO 30
c  if we come to this point no intersection could be found within pmax
      iflag=1
      return

c  call the Newton routine for finding the intersection point
 77   dum=time
      call nwt(a,e,time,iflag)
c  project the perturbation vector onto the plane of intersection
      call project
      flight=time
      call moment(a,ens,15)
      edr=edr+(dt+time-dum)*ens
      edr=edr*2d0*anu/flight

ccccccccccccccccccccccccccccccc
c      write(*,*) time,a(npsfix)
ccccccccccccccccccccccccccccccccc

C
      return
      END

C23456789012345678901234567890123456789012345678901234567890123456789012

      subroutine project
      implicit none
      include 'n_power.for'
c  projection of the perturbation vector E onto plane of intersection
c  A(npsfix)=xsec along the vector field
      integer i,ns,nsc,nnwt,nfix,npsfix
      double precision e,ee,f,a,dtmem,secxtol,secftol,xsec,pmax,ANU,DT
      double precision proj
      parameter(NS=2**NN/6,NSC=(1+NS)**3)
      dimension f(nsc),e(0:NS,0:NS,0:NS),ee(nsc),a(nsc)
      common /section/ secxtol,secftol,xsec,pmax,nnwt,nfix,npsfix
      COMMON /VOR3/e
      common /vffix/ f
      COMMON /VOR1/A
      COMMON /CONST/ANU,DT

      EQUIVALENCE (e(0,0,0),ee(1))

      dtmem=dt
      dt=0
      call march(a,e,0)
      dt=dtmem

c  project
      proj=-ee(npsfix)/f(npsfix)
      do i=1,nsc
         ee(i)=ee(i)+proj*f(i)
      enddo
c  put fixed forcing components to zero
      e(0,0,1)=0d0
      e(1,0,0)=0d0

      return
      end

C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE NWT(y,e,TIME,iflag)
      IMPLICIT NONE
      INCLUDE 'n_power.for'
c     Newton's method for finding intersection point x(fix)=xsec
c     In /vffix/ the vector field is stored upon exiting MARCH
      INTEGER I,II,STEPS,NS,NSC,iflag,nnwt,fix,psfix
      DOUBLE PRECISION Y,F,TIME,FN,secxtol,secftol,xsec,pmax
      DOUBLE PRECISION SM,DTMEM,e,NU,DT
      PARAMETER(SM=1D-8,NS=2**NN/6,NSC=(1+NS)**3)
      DIMENSION Y(NSC),F(NSC),e(nsc)

      common /section/ secxtol,secftol,xsec,pmax,nnwt,fix,psfix
      common /vffix/ f
      COMMON /CONST/ NU,DT

      DTMEM=DT
      DO 11 I=1,NNWT
         FN=0D0
         DO 5 II=1,NSC
            FN=FN+F(II)**2
 5       CONTINUE
         FN=DSQRT(FN)
         IF(DABS(F(PSFIX)/FN).LT.SM) WRITE(*,*) 
     C        "Bad transversality f_fix/|f|=",F(PSFIX)/FN
         dt=(XSEC-Y(PSFIX))/F(PSFIX)
         IF((DABS(DT).LT.SECXTOL).OR.
     C        (DABS(Y(PSFIX)-XSEC).LT.SECFTOL)) GOTO 20
         TIME=TIME+DT

         CALL MARCH(y,e,0)
 11   CONTINUE
      iflag=2
      RETURN

 20   CONTINUE
      DT=DTMEM
      RETURN
      END

C23456789012345678901234567890123456789012345678901234567890123456789012

*PROCESS DOUBLE
C
      SUBROUTINE ZERO(A,N)
C
c
      implicit double precision (a-h,o-z)
c
      DIMENSION A(N)
C
      DO 10 J=1,N
   10 A(J)=0
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE COPY(A,B,N)
C
c
      implicit double precision (a-h,o-z)
c
      DIMENSION A(N),B(N)
C
      DO 10 J=1,N
   10 B(J)=A(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FACT
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND6M1=ND6-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1)
C
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
C
      PI=4*ATAN(1.)
      DO 10 J=1,ND8M1
      SIN4I=1./SIN(4*PI*J/N)
      SIN8I=1./SIN(8*PI*J/N)
      FSINI(J)=0.25*SIN8I
      F1(J)=0.5*(1-0.5*SIN8I)
      F2(J)=0.5*(1+0.5*SIN8I)
      F3(J)=0.5*F1(J)*SIN4I
      F4(J)=-0.5*F2(J)*SIN4I
   10 CONTINUE
      DO 20 J=0,ND4M1
      FCOSI(J)=1./COS(2*PI*J/N)
      FCOS2I(J)=0.5*FCOSI(J)
   20 CONTINUE
      DO 30 J=ND4,NL
      FCOSI(J)=0
   30 FCOS2I(J)=0
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTESSS
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NS=ND6,NSC=(NS+1)**3)
      PARAMETER(NLC=(NL+1)**3,NM=NLC*3/4)
      PARAMETER(N1=(NS+1)**2,N2=(NS+1)*(NL+1),N3=(NL+1)**2)
C
      COMMON /MEMO2/A(1)
      COMMON /MEMO3/B(1)
      COMMON /WORK4/W1(1)
      COMMON /WORK5/W2(1)
      COMMON /WORK6/W3(1)
C
      CALL COR1(A,W1,W2,W3,N3)
      CALL MIXS1(A,B,N3)
      CALL COR2(B,W3,N3)
      CALL MIXC(B,A,N3)
      CALL FFT(A,N3)
      CALL MIXS2(A,B,W1,W2,W3,N3)
      CALL COPY(B,A,N3*(NS+1))
C
      CALL COR1(A,W1,W2,W3,N2)
      CALL MIXS1(A,B,N2)
      CALL COR2(B,W3,N2)
      CALL MIXC(B,A,N2)
      CALL FFT(A,N2)
      CALL MIXS2(A,B,W1,W2,W3,N2)
      CALL COPY(B,A,N2*(NS+1))
C
      CALL COR1(A,W1,W2,W3,N1)
      CALL MIXS1(A,B,N1)
      CALL COR2(B,W3,N1)
      CALL MIXC(B,A,N1)
      CALL FFT(A,N1)
      CALL MIXS2(A,B,W1,W2,W3,N1)
      CALL COPY(B,A,NSC)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTEASS(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NS=ND6,NLP1=NL+1,NLC=NLP1**3)
      PARAMETER(N1=(NS+1)**2,N2=(NS+1)*(NL+1),N3=(NL+1)**2)
C
      COMMON /MEMO3/B(0:NL,N3)
      COMMON /WORK4/W1(1)
      COMMON /WORK5/W2(1)
      COMMON /WORK6/W3(1)
C
      DIMENSION A(0:NS,N3)
C
      CALL ZERO(B,NLP1*N1)
      DO 10 J=0,NS
      DO 10 K=1,N1
   10 B(J,K)=A(J,K)
      CALL MIXA1(B,A,N1)
      CALL COR2(A,W3,N1)
      CALL MIXC(A,B,N1)
      CALL FFT(B,N1)
      CALL MIXA2(B,A,W3,N1)
C
      CALL ZERO(B,NLP1*N2)
      DO 20 J=0,NS
      DO 20 K=1,N2
   20 B(J,K)=A(J,K)
      CALL COR1(B,W1,W2,W3,N2)
      CALL MIXS1(B,A,N2)
      CALL COR2(A,W3,N2)
      CALL MIXC(A,B,N2)
      CALL FFT(B,N2)
      CALL MIXS2(B,A,W1,W2,W3,N2)
C
      CALL ZERO(B,NLC)
      DO 30 J=0,NS
      DO 30 K=1,N3
   30 B(J,K)=A(J,K)
      CALL COR1(B,W1,W2,W3,N3)
      CALL MIXS1(B,A,N3)
      CALL COR2(A,W3,N3)
      CALL MIXC(A,B,N3)
      CALL FFT(B,N3)
      CALL MIXS2(B,A,W1,W2,W3,N3)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTESAA(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NS=ND6,NLC=(NL+1)**3,NSP1=NS+1,NSC=NSP1**3)
      PARAMETER(N1=(NS+1)**2,N2=(NS+1)*(NL+1),N3=(NL+1)**2)
C
      COMMON /MEMO3/B(1)
      COMMON /WORK4/W1(1)
      COMMON /WORK5/W2(1)
      COMMON /WORK6/W3(1)
C
      DIMENSION A(1)
C
      CALL COR1(A,W1,W2,W3,N3)
      CALL MIXS1(A,B,N3)
      CALL COR2(B,W3,N3)
      CALL MIXC(B,A,N3)
      CALL FFT(A,N3)
      CALL MIXS2(A,B,W1,W2,W3,N3)
      CALL COPY(B,A,N3*NSP1)
C
      CALL MIXA1(A,B,N2)
      CALL COR2(B,W3,N2)
      CALL MIXC(B,A,N2)
      CALL FFT(A,N2)
      CALL MIXA2(A,B,W3,N2)
      CALL COPY(B,A,N2*NSP1)
C
      CALL MIXA1(A,B,N1)
      CALL COR2(B,W3,N1)
      CALL MIXC(B,A,N1)
      CALL FFT(A,N1)
      CALL MIXA2(A,B,W3,N1)
      CALL COPY(B,A,NSC)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTQSS(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,NL=ND4+1,N2=(NL+1)**2)
C
      COMMON /CONST4/B(1)
C
      DIMENSION A(1)
C
      CALL MIXY1(A)
      CALL FTESS(A)
      CALL PROD(A,B,N2)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXY1(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,NSP1=NS+1,N2=NSP1**2)
C
      COMMON /WORK2/B(N2)
C
      DIMENSION A(N2)
C
      CALL COPY(A,B,N2)
C
      DO 10 J=NSP1+2,N2
   10 A(J)=B(J)+B(J-1)+B(J-NSP1)+B(J-NSP1-1)
C
      DO 20 J=2,ND6+1
   20 A(J)=2*(B(J)+B(J-1))
C
      DO 30 J=NSP1+1,N2,NSP1
   30 A(J)=2*(B(J)+B(J-NSP1))
C
      A(1)=4*B(1)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTESS(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1,NS=ND6)
      PARAMETER(N1=NS+1,N2=NL+1)
C
      COMMON /WORK4/W1(N2),W2(N2),W3(N2)
      COMMON /WORK5/B(0:NL,0:NL)
C
      DIMENSION A(0:NS,0:NL)
C
      CALL ZERO(B,N1*N2)
      DO 10 K=0,NS
      DO 10 J=0,NS
   10 B(J,K)=A(J,K)
      CALL COR1(B,W1,W2,W3,N1)
      CALL MIXS1(B,A,N1)
      CALL COR2(A,W3,N1)
      CALL MIXC(A,B,N1)
      CALL FFT(B,N1)
      CALL MIXS2(B,A,W1,W2,W3,N1)
C
      CALL ZERO(B,N2*N2)
      DO 20 K=0,NL
      DO 20 J=0,NS
   20 B(J,K)=A(J,K)
      CALL COR1(B,W1,W2,W3,N2)
      CALL MIXS1(B,A,N2)
      CALL COR2(A,W3,N2)
      CALL MIXC(A,B,N2)
      CALL FFT(B,N2)
      CALL MIXS2(B,A,W1,W2,W3,N2)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTOSA(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,N3=(NL+1)**2)
C
      COMMON /CONST5/C(1)
C
      DIMENSION AK(0:ND4),A(0:NL,0:NL)
C
      DO 10 J=0,ND4
   10 AK(J)=A(J,ND4)
C
      CALL FTOS(AK)
C
      CALL PROD(A,C,N3)
      CALL FTESA1(A)
      CALL MIXX2(AK)
      CALL ADIM(A)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXX2(AK)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NS=ND6,NSP1=NS+1,NSP2=NS+2,NSP3=NS+3)
      PARAMETER(N2=NSP2**2,M=NSP2*NSP1)
C
      COMMON /WORK1/ A(0:NSP1,0:NSP1)
      COMMON /WORK4/B(0:NSP1,0:NSP1)
C
      DIMENSION AK(0:ND4)
      DIMENSION AA(N2),BB(N2)
C
      EQUIVALENCE (A(0,0),AA(1)),(B(0,0),BB(1))
C
      DO 10 J=1,M-1
   10 BB(J)=0.25*(AA(J)+AA(J+1)+AA(J+NSP2)+AA(J+NSP3))
C
      DO 20 K=0,ND6,2
      DO 20 J=0,ND6
   20 A(J,K)=B(J,K)-AK(J)
C
      DO 30 K=1,ND6,2
      DO 30 J=0,ND6
   30 A(J,K)=B(J,K)+AK(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTOS(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND6M1=ND6-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1)
C
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
C
      DIMENSION A(0:ND4)
C
      DO 10 J=0,ND4
   10 A(J)=A(J)*FCOSI(J)
      CALL FTES(A,2)
      AND4M1=A(ND4M1)
      DO 20 J=0,ND4M1
   20 A(J)=A(J)+A(J+1)
      A(ND4)=A(ND4)+AND4M1
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTES(A,ID)
c
      implicit double precision (a-h,o-z)
c
      include 'n_power.for'
C
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND6M1=ND6-1,ND8M1=ND8-1)
      PARAMETER(NVW=ND8*(NN-1))
      PARAMETER(NL=ND4+1)
C
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      COMMON /WORK4/W4(0:ND4M1),W5(0:ND8M1),W6(0:ND8M1)
c      COMMON /WORK7/VW(NVW),IVW(NVW)
      COMMON /WORK7/ip(0:nd8)
      common /work8/wwww(0:nd8*5/2-1)
      dimension abc(2*nd8)
      DIMENSION A(0:ND4)
C
      IF(ID.NE.2) GO TO 4
      A(ND4)=2*A(ND4)
    4 S1=0
      S2=0
      DO 10 J=0,ND8M1
      J2=J*2
      S1=S1+A(J2)
   10 S2=S2+A(J2+1)
      W1=2*(S1+S2)-A(0)+A(ND4)
      W2=2*(S1-S2)-A(0)+A(ND4)
      W4(0)=A(0)
      W4(ND8)=A(ND4)
      DO 30 J=1,ND8M1
      J2=J*2
      S1=A(J2)
      S2=A(J2+1)-A(J2-1)
      W4(J)=S1+S2
   30 W4(ND4-J)=S1-S2
      W3=0
      DO 40 J=1,ND4M1,2
   40 W3=W3+W4(J-1)-W4(J)
      W5(0)=W4(0)
      W6(0)=W4(ND4M1)-W4(1)
      DO 50 J=1,ND8M1
      J2=J*2
      W5(J)=W4(J2)
   50 W6(J)=W4(J2-1)-W4(J2+1)
c-------------------------------------------
c      CALL DVCFT1(W5,W6,ND8,-1,1,VW,IVW,ICON)
c
      call two2one(w5,w6,abc,nd8)
      call cdft(2*nd8,1,abc,ip,wwww)
      call one2two(abc,w5,w6,nd8)
c-------------------------------------------
      DO 70 J=1,ND8M1
      S1=F1(J)*W5(J)+F2(J)*W5(ND8-J)
      S2=F3(J)*W6(J)+F4(J)*W6(ND8-J)
      A(J)=S1-S2
   70 A(ND4-J)=S1+S2
      A(0)=W1
      A(ND4)=W2
      A(ND8)=W3
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTESA1(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NS=ND6,NSP1=NS+1)
      PARAMETER(N1=NS+2,N2=NL+1)
C
      COMMON /WORK4/W1(N2),W2(N2),W3(N2)
      COMMON /WORK5/B(1)
C
      DIMENSION A(1),B1(0:NL,0:NSP1),B2(0:NSP1,0:NSP1)
C
      EQUIVALENCE (B(1),B1(0,0)),(B(1),B2(0,0))
C
      CALL COR1(A,W1,W2,W3,N2)
      CALL MIXS1(A,B,N2)
      CALL COR2(B,W3,N2)
      CALL MIXC(B,A,N2)
      CALL FFT(A,N2)
      CALL MIXS2(A,B,W1,W2,W3,N2)
      CALL COPY(B1,A,N1*N2)
C
      CALL MIXA1(A,B,N1)
      CALL COR2(B,W3,N1)
      CALL MIXC(B,A,N1)
      CALL FFT(A,N1)
      CALL MIXA2(A,B,W3,N1)
      CALL COPY(B2,A,N1*N1)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTESA2(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NS=ND6)
      PARAMETER(N1=NS+1,N2=NL+1)
C
      COMMON /WORK4/W1(N2),W2(N2),W3(N2)
      COMMON /WORK5/B(1)
C
      DIMENSION A(1),B1(0:NL,0:NS),B2(0:NS,0:NS)
C
      EQUIVALENCE (B(1),B1(0,0)),(B(1),B2(0,0))
C
      CALL COR1(A,W1,W2,W3,N2)
      CALL MIXS1(A,B,N2)
      CALL COR2(B,W3,N2)
      CALL MIXC(B,A,N2)
      CALL FFT(A,N2)
      CALL MIXS2(A,B,W1,W2,W3,N2)
      CALL COPY(B1,A,N1*N2)
C
      CALL MIXA1(A,B,N1)
      CALL COR2(B,W3,N1)
      CALL MIXC(B,A,N1)
      CALL FFT(A,N1)
      CALL MIXA2(A,B,W3,N1)
      CALL COPY(B2,A,N1*N1)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE COR1(A,W1,W2,W3,M)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
C
      DIMENSION A(0:NL,M),W1(M),W2(M),W3(M)
C
      CALL ZERO(W2,M)
      CALL ZERO(W3,M)
C
      DO 10 J=1,ND4M1,2
      JM1=J-1
      DO 10 K=1,M
      W2(K)=W2(K)+A(JM1,K)
   10 W3(K)=W3(K)+A(J,K)
      DO 20 K=1,M
      W1(K)=2*(W2(K)+W3(K))-A(0,K)+A(ND4,K)
   20 W2(K)=2*(W2(K)-W3(K))-A(0,K)+A(ND4,K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXS1(A,B,M)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
C
      DIMENSION A(0:NL,M),B(M,0:NL)
C
      DO 10 K=1,M
      B(K,0)=A(0,K)
   10 B(K,ND8)=A(ND4,K)
      DO 20 J=1,ND8M1
      J2=J*2
      J2P1=J2+1
      J2M1=J2-1
      ND4MJ=ND4-J
      DO 20 K=1,M
      S2=A(J2P1,K)-A(J2M1,K)
      B(K,J)=A(J2,K)+S2
   20 B(K,ND4MJ)=A(J2,K)-S2
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE COR2(A,W,M)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
C
      DIMENSION A(M,0:NL),W(M)
C
      CALL ZERO(W,M)
C
      DO 10 J=1,ND4M1,2
      JM1=J-1
      DO 10 K=1,M
   10 W(K)=W(K)+A(K,JM1)-A(K,J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXC(A,B,M)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
C
      DIMENSION A(M,0:NL),B(0:NL,M)
C
      DO 10 K=1,M
      B(0,K)=A(K,0)
   10 B(ND8,K)=A(K,ND4M1)-A(K,1)
      DO 20 J=1,ND8M1
      J2=J*2
      J2P1=J2+1
      J2M1=J2-1
      ND8PJ=ND8+J
      DO 20 K=1,M
      B(J,K)=A(K,J2)
   20 B(ND8PJ,K)=A(K,J2M1)-A(K,J2P1)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FFT(A,M)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'

      integer omp_get_thread_num

      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
      PARAMETER(NVW=ND8*(NN-1))
C
c      COMMON /WORK7/VW(NVW),IVW(NVW)
      COMMON /WORK7/ip(0:nd8)
      common /work8/wwww(0:nd8*5/2-1)
      dimension abc(2*nd8)

      DIMENSION A(0:NL,M)

!$OMP PARALLEL DO if(m.gt.1000)
!$OMP& schedule(static)
!$OMP& default(shared) private(k,abc)
!$OMP& firstprivate(ip,www) lastprivate(ip,www) 
      DO 10 K=1,M
         call two2one(a(0,k),a(nd8,k),abc,nd8)
         call cdft(2*nd8,1,abc,ip,wwww)
         call one2two(abc,a(0,k),a(nd8,k),nd8)
10    continue
!$OMP END PARALLEL DO

      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXS2(A,B,W1,W2,W3,M)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
C
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
C
      DIMENSION A(0:NL,M),B(M,0:NL),W1(M),W2(M),W3(M)
C
      DO 10 J=1,ND8M1
      ND8PJ=ND8+J
      ND8MJ=ND8-J
      ND4MJ=ND4-J
      DO 10 K=1,M
      S1=F1(J)*A(J,K)+F2(J)*A(ND8MJ,K)
      S2=F3(J)*A(ND8PJ,K)+F4(J)*A(ND4MJ,K)
      B(K,J)=S1-S2
   10 B(K,ND4MJ)=S1+S2
      DO 20 K=1,M
      B(K,0)=W1(K)
      B(K,ND4)=W2(K)
   20 B(K,ND8)=W3(K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE ADIM(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NS=ND6,NSP1=NS+1)
C
      COMMON /WORK4/B(0:NS,0:NS)
C
      DIMENSION A(0:NSP1,0:NSP1)
C
      DO 10 K=0,NS
      DO 10 J=0,NS
   10 B(J,K)=A(J,K)
C
      CALL COPY(B,A,NSP1**2)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXA1(A,B,M)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
C
      DIMENSION A(0:NL,M),B(M,0:NL)
C
      DO 10 K=1,M
      B(K,0)=2*A(1,K)
   10 B(K,ND8)=-2*A(ND4M1,K)
      DO 20 J=1,ND8M1
      J2=J*2
      J2P1=J2+1
      J2M1=J2-1
      ND4MJ=ND4-J
      DO 20 K=1,M
      S2=A(J2P1,K)-A(J2M1,K)
      B(K,J)=S2+A(J2,K)
   20 B(K,ND4MJ)=S2-A(J2,K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXA2(A,B,W,M)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1,NLC=(NL+1)**3)
C
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      DIMENSION A(0:NL,M),B(M,0:NL),W(M)
C
      DO 10 J=1,ND8M1
      ND8PJ=ND8+J
      ND8MJ=ND8-J
      ND4MJ=ND4-J
      DO 10 K=1,M
      S1=F2(J)*A(ND4MJ,K)-F1(J)*A(ND8PJ,K)
      S2=F4(J)*A(ND8MJ,K)-F3(J)*A(J,K)
      B(K,J)=S1+S2
   10 B(K,ND4MJ)=S2-S1
      DO 20 K=1,M
      B(K,0)=0
      B(K,ND4)=0
   20 B(K,ND8)=-0.5*W(K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE UT
C
C  *********************
C  *  B: WORKING AREA  *
C  *********************
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,NLP1=N/4+2,NLP1S=NLP1**2,NLC=NLP1**3)
C
      COMMON /MEMO2/A(1)
      COMMON /MEMO3/B(1)
C
      CALL CYCLE(B,A)
      CALL PROD(B,A,NLC)
      CALL CYCLE(A,B)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE UT2
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,NL=N/4+1)
C
      COMMON /MEMO1/A(0:NL,0:NL,0:NL)
      COMMON /MEMO2/B(0:NL,0:NL,0:NL)
      COMMON /MEMO3/C(0:NL,0:NL,0:NL)
C
      DO 10 J=0,NL
      DO 10 K=0,NL
      DO 10 L=0,NL
   10 C(J,K,L)=A(K,L,J)*B(L,J,K)+B(K,L,J)*A(L,J,K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE PROD(A,B,N)
C
c
      implicit double precision (a-h,o-z)
c
      DIMENSION A(N),B(N)
C
      DO 10 J=1,N
   10 A(J)=A(J)*B(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE CYCLE(A,B)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,NLP1=N/4+2,NLP1S=NLP1**2)
C
      DIMENSION A(NLP1,NLP1S),B(NLP1S,NLP1)
C
      DO 10 J=1,NLP1
      DO 10 K=1,NLP1S
   10 A(J,K)=B(K,J)
      RETURN
C
      END
*PROCESS DOUBLE
      SUBROUTINE SQUARE(A,N)
C
c
      implicit double precision (a-h,o-z)
c
      DIMENSION A(N)
C
      DO 10 J=1,N
   10 A(J)=A(J)*A(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE US2
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,NL=N/4+1,NLC=(NL+1)**3)
C
      COMMON /MEMO1/A(NLC)
      COMMON /MEMO2/B(NLC)
C
      CALL PROD(A,B,NLC)
      DO 10 J=1,NLC
   10 A(J)=2*A(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE NLS
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      COMMON /MEMO2/B(0:NS,0:NS,0:NS),C(0:NS,0:NS,0:NS),
     c               D(0:NS,0:NS,0:NS)
      COMMON /CONST2/E(0:NS,0:NS)
C
      CALL COPY(B,C,NSC)
      DO 10 L=0,NS
      DO 10 K=0,L
      DO 10 J=0,NS
   10 C(J,K,L)=C(J,L,K)
C
      CALL COPY(B,D,NSC)
      DO 20 K=0,NS
      DO 20 J=0,NS
   20 D(J,K,K)=0
C
      DO 30 L=0,NS
      DO 30 K=L,NS
      DO 30 J=0,NS
   30 D(J,K,L)=-D(J,L,K)
C
      DO 40 K=0,NS
      DO 40 L=0,K
      DO 40 J=0,NS
   40 B(J,K,L)=C(L,J,K)-C(K,J,L)
C
      DO 50 L=0,NS
      DO 50 K=0,L
      DO 50 J=0,NS
   50 B(J,K,L)=D(L,J,K)+D(K,J,L)
C
      CALL EXT(E,C)
      CALL PROD(B,C,NSC)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE BASE
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NS=ND6,NSC=(NS+1)**3)
C
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      COMMON /CONST1/A(0:NS,0:NS),B(0:NS,0:NS)
      COMMON /CONST2/C(0:NS,0:NS)
      COMMON /CONST3/D(0:NS,0:NS),E(0:NS,0:NS),F(0:NS,0:NS)
      COMMON /CONST4/G(0:NL,0:NL)
      COMMON /CONST5/H(0:NL,0:NL)
      COMMON /CONST6/P(0:NS,0:NS),Q(0:NS)
C
      FF=16./N**3
      DO 10 L=0,NS
      DO 10 K=0,L
      A(K,L)=(K+0.5)*0.5
   10 B(K,L)=(L+0.5)*0.5
C
      DO 20 L=0,NS
      DO 20 K=L,NS
      A(K,L)=K*0.5
   20 B(K,L)=L*0.5
C
      DO 30 L=0,NS
      DO 30 K=L,NS
   30 C(K,L)=FF*K*L
C
      DO 40 L=0,NS
      DO 40 K=0,L
   40 C(K,L)=FF*(K+0.5)*(L+0.5)
C
      DO 50 L=0,NS
      DO 50 K=L,NS
      D(K,L)=K*FF
   50 E(K,L)=L*FF
C
      DO 60 L=0,NS
      DO 60 K=0,L
      D(K,L)=(K+0.5)*FF
   60 E(K,L)=(L+0.5)*FF
C
      DO 70 L=0,NS
      DO 70 K=L,NS
   70 F(K,L)=(K*K-L*L)*FF
C
      DO 80 L=0,NS
      DO 80 K=0,L
   80 F(K,L)=((K+0.5)**2-(L+0.5)**2)*FF
C
      DO 90 J=0,NL
      DO 90 K=0,NL
   90 G(J,K)=FCOS2I(J)*FCOS2I(K)
C
      DO 100 J=0,NL
      DO 100 K=0,NL
  100 H(J,K)=4*G(J,K)
C
      DO 110 K=0,NS
      DO 110 L=0,K
  110 P(K,L)=K*K+L*L
C
      DO 120 L=0,NS
      DO 120 K=0,L
  120 P(K,L)=(K+0.5)**2+(L+0.5)**2+0.25
C
      DO 130 J=0,NS
  130 Q(J)=J*J
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE EXT(A,B)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,N1=(NS+1)**2)
C
      DIMENSION A(N1),B(0:NS,N1)
C
      DO 10 J=0,NS
      DO 10 K=1,N1
   10 B(J,K)=A(K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE NLT
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      COMMON /MEMO2/A(NSC),D(0:NS,0:NS,0:NS),E(0:NS,0:NS,0:NS)
      COMMON /MEMO3/P(NSC),Q(0:NS,0:NS,0:NS),C(0:NS,0:NS,0:NS)
      COMMON /CONST3/F(0:NS,0:NS),G(0:NS,0:NS),H(0:NS,0:NS)
C
      DIMENSION CC(NSC),DD(NSC),EE(NSC),QQ(NSC)
C
      EQUIVALENCE (C(0,0,0),CC(1)),(D(0,0,0),DD(1)),
     c             (E(0,0,0),EE(1)),(Q(0,0,0),QQ(1))
C
      CALL COPY(A,C,NSC)
      DO 10 L=0,NS
      DO 10 K=0,L
      DO 10 J=0,NS
   10 C(J,K,L)=C(J,L,K)
C
      CALL COPY(A,D,NSC)
      DO 20 K=0,NS
      DO 20 J=0,NS
   20 D(J,K,K)=0
      DO 30 L=0,NS
      DO 30 K=L,NS
      DO 30 J=0,NS
   30 D(J,K,L)=-D(J,L,K)
C
      DO 40 L=0,NS
      DO 40 K=L,NS
      DO 40 J=0,NS
   40 E(J,K,L)=C(K,L,J)
C
      DO 50 L=0,NS
      DO 50 K=0,L
      DO 50 J=0,NS
   50 E(J,K,L)=D(K,L,J)
C
      CALL EXT(F,P)
C
      DO 60 L=0,NS
      DO 60 K=L,NS
      DO 60 J=0,NS
   60 Q(J,K,L)=C(L,J,K)
C
      DO 70 K=0,NS
      DO 70 L=K,NS
      DO 70 J=0,NS
   70 Q(J,K,L)=D(L,J,K)
C
      CALL EXT(G,C)
C
      DO 80 L=0,NS
      DO 80 K=L,NS
      DO 80 J=0,NS
   80 D(J,K,L)=J
C
      DO 90 L=0,NS
      DO 90 K=0,L
      DO 90 J=0,NS
   90 D(J,K,L)=J+0.5
C
      CALL PROD(CC,QQ,NSC)
      CALL EXT(H,QQ)
C
      DO 100 J=1,NSC
  100 A(J)=A(J)*QQ(J)+(EE(J)*P(J)-CC(J))*DD(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE ASUM(A,B,N)
C
c
      implicit double precision (a-h,o-z)
c
      DIMENSION A(N),B(N)
C
      DO 10 J=1,N
   10 A(J)=A(J)+B(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE DIF(A,B,N)
C
c
      implicit double precision (a-h,o-z)
c
      DIMENSION A(N),B(N)
C
      DO 10 J=1,N
   10 A(J)=A(J)-B(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE WU(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      COMMON /MEMO2/F(NSC),D(0:NS,0:NS,0:NS)
      COMMON /MEMO3/E(0:NS,0:NS,0:NS),B(0:NS,0:NS,0:NS)
     c     ,C(0:NS,0:NS,0:NS)
      COMMON /CONST1/G(0:NS,0:NS),H(0:NS,0:NS)
C
      DIMENSION A(0:NS,0:NS,0:NS),DD(NSC)
C
      EQUIVALENCE (D(0,0,0),DD(1))
C
      CALL WUFAC
C
      CALL COPY(A,C,NSC)
      DO 10 L=0,NS
      DO 10 K=L,NS
      DO 10 J=0,NS
   10 C(J,K,L)=C(J,L,K)
C
      CALL COPY(A,B,NSC)
      DO 20 L=0,NS
      DO 20 K=0,L
      DO 20 J=0,NS
   20 B(J,K,L)=-B(J,L,K)
      DO 30 K=0,NS
      DO 30 J=0,NS
   30 B(J,K,K)=0
C
      DO 40 L=0,NS
      DO 40 K=0,L
      DO 40 J=0,NS
   40 E(J,K,L)=C(L,J,K)
C
      DO 50 L=0,NS
      DO 50 K=L,NS
      DO 50 J=0,NS
   50 E(J,K,L)=B(L,J,K)
C
      CALL EXT(G,D)
      CALL PROD(D,E,NSC)
C
      DO 60 L=0,NS
      DO 60 K=0,L
      DO 60 J=0,NS
   60 E(J,K,L)=C(K,L,J)
C
      DO 70 L=0,NS
      DO 70 K=L,NS
      DO 70 J=0,NS
   70 E(J,K,L)=B(K,L,J)
C
      CALL EXT(H,C)
      CALL PROD(C,E,NSC)
C
      CALL DIF(D,C,NSC)
C
      DO 80 J=2,NSC
   80 F(J)=DD(J)/F(J)
      F(1)=0
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE WUFAC
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6)
C
      COMMON /MEMO2/ A(0:NS,0:NS,0:NS)
C
      DO 10 K=0,NS
      DO 10 J=0,NS
   10 A(J,K,K)=J*J+2*K*K
      A(0,0,0)=1.D-50
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE SPLIT(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6,ND6M1=ND6-1)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      DIMENSION A(0:NS,0:NS,0:NS)
C
      DO 10 L=0,ND6M1
      DO 10 K=L+1,ND6
      DO 10 J=0,ND6
   10 A(J,K,L)=0.5*(A(J,K,L)+A(J,L,K))
C
      DO 20 L=1,ND6
      DO 20 K=0,L-1
      DO 20 J=0,ND6
   20 A(J,K,L)=A(J,K,L)-A(J,L,K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE UNIFY(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6,ND6M1=ND6-1)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      DIMENSION A(0:NS,0:NS,0:NS)
C
      DO 10 L=0,ND6M1
      DO 10 K=L+1,ND6
      DO 10 J=0,ND6
   10 A(J,K,L)=0.5*(A(J,K,L)-A(J,L,K))
C
      DO 20 K=1,ND6
      DO 20 L=0,K-1
      DO 20 J=0,ND6
   20 A(J,L,K)=A(J,L,K)+A(J,K,L)
C
      DO 30 K=0,ND6
      DO 30 J=0,ND6
   30 A(J,K,K)=0.5*A(J,K,K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE SRF
C
c
      implicit double precision (a-h,o-z)
c
      COMMON /MEMO2/A(1)
C
      CALL MIXZ1
      CALL FTESSS
      CALL SPLIT(A)
      CALL MIXZ2
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXZ1
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,NL=ND4+1,ND8=N/8,NLC=(NL+1)**3)
C
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      COMMON /MEMO2/B(0:NL,0:NL,0:NL)
      COMMON /CONST5/C(0:NL,0:NL)
C
      DIMENSION BB(NLC)
C
      EQUIVALENCE (B(0,0,0),BB(1))
C
      DO 40 K=0,ND4
      DO 40 L=0,K
      DO 40 J=0,ND4
      SUM=B(J,K,L)+B(J,L,K)
      DIF=(B(J,K,L)-B(J,L,K))*FCOSI(J)*C(K,L)
      B(J,K,L)=SUM+DIF
   40 B(J,L,K)=SUM-DIF
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXZ2
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6,ND6M1=ND6-1)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      COMMON /MEMO2/B(0:NS,0:NS,0:NS),C(0:NS,0:NS,0:NS)
C
      DO 120 L=1,ND6
      DO 120 K=0,L-1
      DO 120 J=0,ND6M1
  120 C(J,K,L)=0.125*(B(J,K,L)+B(J+1,K,L))
C
      DO 130 L=2,ND6M1
      LP1=L+1
      DO 130 K=0,L-2
      KP1=K+1
      DO 130 J=0,ND6M1
  130 B(J,K,L)=C(J,K,L)+C(J,KP1,L)+C(J,K,LP1)+C(J,KP1,LP1)
C
      DO 140 L=1,ND6M1
      LP1=L+1
      LM1=L-1
      DO 140 J=0,ND6M1
  140 B(J,LM1,L)=C(J,LM1,L)+C(J,LM1,LP1)+C(J,L,LP1)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE TRF
C
c
      implicit double precision (a-h,o-z)
c
      COMMON /MEMO2/ A(1)
      COMMON /WORK1/ B(1)
      COMMON /WORK2/ C(1)
      COMMON /WORK3/ D(1)
C
C
      CALL CORT1(A,B,C,D)
      CALL MIXZ1
      CALL FTESAA(A)
      CALL CORT2(A,C,D)
      CALL SPLIT(A)
      CALL MIXZ2
      CALL CORT3(A,B)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE CORT1(A,B,C,D)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,NL=ND4+1)
C
      DIMENSION A(0:NL,0:NL,0:NL),B(0:NL,0:NL),C(0:NL,0:NL),D(0:NL)
C
      DO 10 L=0,NL
   10 D(L)=2*A(L,ND4,ND4)
C
      DO 20 L=0,NL
      DO 20 J=0,NL
      B(J,L)=A(J,ND4,L)-A(J,L,ND4)
   20 C(J,L)=A(J,ND4,L)+A(J,L,ND4)
C
      CALL FTOSA(B)
      CALL FTESA2(C)
      CALL FTES(D,1)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE CORT2(A,C,D)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      DIMENSION A(0:NS,0:NS,0:NS),C(0:NS,0:NS),D(0:NS)
C
      DO 10 L=0,ND6,2
      DO 10 K=0,ND6,2
      DO 10 J=0,ND6
   10 A(J,K,L)=D(J)-A(J,K,L)-C(J,L)-C(J,K)
C
      DO 20 L=1,ND6,2
      DO 20 K=0,ND6,2
      DO 20 J=0,ND6
   20 A(J,K,L)=C(J,K)-A(J,K,L)-C(J,L)-D(J)
C
      DO 30 L=1,ND6,2
      DO 30 K=1,ND6,2
      DO 30 J=0,ND6
   30 A(J,K,L)=D(J)-A(J,K,L)+C(J,L)+C(J,K)
C
      DO 40 L=0,ND6,2
      DO 40 K=1,ND6,2
      DO 40 J=0,ND6
   40 A(J,K,L)=C(J,L)-A(J,K,L)-D(J)-C(J,K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE CORT3(A,B)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      DIMENSION A(0:NS,0:NS,0:NS),B(0:NS,0:NS)
C
      DO 10 L=1,ND6,2
      DO 10 K=0,L-1,2
      DO 10 J=0,ND6
   10 A(J,K,L)=A(J,K,L)+B(J,L)+B(J,K)
C
      DO 20 L=3,ND6,2
      DO 20 K=1,L-1,2
      DO 20 J=0,ND6
   20 A(J,K,L)=A(J,K,L)-B(J,L)+B(J,K)
C
      DO 30 L=2,ND6,2
      DO 30 K=0,L-1,2
      DO 30 J=0,ND6
   30 A(J,K,L)=A(J,K,L)+B(J,L)-B(J,K)
C
      DO 40 L=2,ND6,2
      DO 40 K=1,L-1,2
      DO 40 J=0,ND6
   40 A(J,K,L)=A(J,K,L)-B(J,L)-B(J,K)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE UFR(A)
C
c
      implicit double precision (a-h,o-z)
c
      COMMON /WORK1/ B(1)
      DIMENSION A(1)
C
      CALL CORU1(A,B)
      CALL MIXU1(A)
      CALL UNIFY(A)
      CALL FTEASS(A)
      CALL MIXU2
      CALL CORU2(A,B)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE CORU1(A,B)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,NS=ND6)
C
      DIMENSION A(0:NS,0:NS,0:NS),B(0:NS,0:NS)
C
      CALL ZERO(B,(NS+1)*(NS+1))
C
      DO 20 L=1,ND6
      DO 20 K=0,L-1
      S1=0
      S2=0
      DO 10 J=1,NS,2
      S1=S1+A(J-1,K,L)
   10 S2=S2+A(J,K,L)
      B(K,L)=S1-S2
   20 B(L,K)=S2-S1
C
      CALL FTQSS(B)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXU1(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND6M1=ND6-1,ND6M2=ND6-2)
      PARAMETER(NL=ND4+1,NS=ND6,NSC=(NS+1)**3)
C
      COMMON /MEMO3/B(0:NS,0:NS,0:NS)
      DIMENSION A(0:NS,0:NS,0:NS)
C
      DO 10 L=1,ND6
      DO 10 K=0,L-1
      DO 10 J=1,ND6
   10 B(J,K,L)=A(J,K,L)+A(J-1,K,L)
C
      DO 30 L=3,ND6
      LM1=L-1
      DO 30 K=1,L-2
      KM1=K-1
      DO 30 J=0,ND6
   30 A(J,K,L)=B(J,K,L)+B(J,K,LM1)+B(J,KM1,L)+B(J,KM1,LM1)
C
      DO 40 L=2,ND6
      LM1=L-1
      LM2=L-2
      DO 40 J=0,ND6
      A(J,0,L)=2*(B(J,0,L)+B(J,0,LM1))
   40 A(J,LM1,L)=B(J,LM1,L)+B(J,LM2,L)+B(J,LM2,LM1)
C
      DO 50 J=0,ND6
      A(J,0,1)=2*B(J,0,1)
   50 A(J,ND6M1,ND6)=B(J,ND6M2,ND6M1)
C
      DO 60 K=0,NS
      DO 60 L=0,NS
   60 A(0,K,L)=0
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MIXU2
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND4M1=ND4-1)
      PARAMETER(NL=ND4+1)
C
      COMMON /MEMO2/A(0:NL,0:NL,0:NL)
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL)
      COMMON /CONST4/B(0:NL,0:NL)
C
      DO 10 L=0,ND4
      DO 10 K=L,ND4
      DO 10 J=0,ND4
      TEMP1=A(J,K,L)+A(J,L,K)
      TEMP2=(A(J,K,L)-A(J,L,K))*FCOS2I(J)*B(K,L)
      A(J,K,L)=TEMP1+TEMP2
   10 A(J,L,K)=TEMP1-TEMP2
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE CORU2(A,B)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,NL=ND4+1,N3=(NL+1)**2)
C
      DIMENSION A(0:NL,N3),B(N3)
C
      DO 10 J=1,N3
   10 A(ND4,J)=A(ND4,J)-2*B(J)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE STORE(A,IFILE)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      COMMON /TIME/TIME,ISTEP
C
      DIMENSION A(NSC)
C
      WRITE(IFILE,*) TIME,A
C
c 90   format(e14.7,1X,<nsc>e14.7)
c      WRITE(6,1000) IFILE,TIME
C
 1000 FORMAT(' VORTICITY WAS WRITTEN IN FILE(',I2,') AT TIME=',1PE14.7)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MARCH(A,E,J)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(NL=ND4+1,NS=ND6,NLC=(NL+1)**3,NSC=(NS+1)**3)
C
      COMMON /MEMO1/C(1)
      COMMON /MEMO2/D(1)
      COMMON /MEMO3/B(1)
      COMMON /MEMO4/H(1)
      COMMON /TIME/TIME,ISTEP
      COMMON /FORCE/IFORCE,IFIX
      COMMON /PARA/IORTH,IHEL,IOUT,IEN,ISTEP2

      common /choose/ ichoose
C
      DIMENSION A(NSC),E(NSC)
C
      DO 60 I=1,4
      ID=I
      CALL WU(A)
      CALL UFR(D)
      IF(I.NE.1.OR.MOD((J+1),IHEL).NE.0) GO TO 30
      CALL HELCTY(B,D,A,18)
   30 CALL COPY(D,C,NLC)
      CALL SQUARE(D,NLC)
      CALL SRF
      CALL NLS
      CALL COPY(D,H,NSC)
      CALL COPY(C,D,NLC)
      CALL UT
      CALL TRF
      CALL NLT
      CALL ASUM(H,D,NSC)
      CALL WNUM
c skip over the next 17 lines to integrate only the
c vorticity, not the perturbation
      if(ichoose.eq.0) goto 1
      CALL WU(E)
      CALL UFR(D)
      IF(I.NE.1.OR.MOD((J+1),IHEL).NE.0) GO TO 50
      CALL HELCTY(B,D,E,19)
      IF(ISTEP.GT.ISTEP2) STOP
   50 CALL UT2
      CALL US2
      CALL COPY(B,D,NLC)
      CALL TRF
      CALL NLT
      CALL COPY(D,B,NSC)
      CALL COPY(C,D,NLC)
      CALL COPY(B,C,NSC)
      CALL SRF
      CALL NLS
      CALL ASUM(C,D,NSC)
      CALL WNUM
c resume here to integrate the vorticity field
 1    CALL STEP(A,E,ID)
      IF(IFIX.EQ.1) CALL FIX(A,E,23,24)
   60 CONTINUE
C
c      IF(IFIX.EQ.1) CALL FIX(A,E,23,24)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FIX(A,E,IFILE1,IFILE2)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      COMMON /TIME/TIME,ISTEP
C
      DIMENSION A(0:NS,0:NS,0:NS),E(0:NS,0:NS,0:NS)
C
c      WRITE(IFILE1,*) ISTEP,TIME,A(0,0,1),A(1,0,0)
C
      A(0,0,1)=-0.375
      A(1,0,0)=0.25
C
c      WRITE(IFILE2,*) ISTEP,TIME,E(0,0,1),E(1,0,0)
C
      E(0,0,1)=0.
      E(1,0,0)=0.
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE STEP(A,E,ID)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      COMMON /VOR2/B(NSC)
      COMMON /VOR4/H(NSC)
      COMMON /MEMO1/G(NSC)
      COMMON /MEMO2/D(NSC)
      COMMON /MEMO4/C(NSC)
      COMMON /CONST/ANU,DT
      COMMON /FORCE/IFORCE,IFIX
      common /section/ secxtol,secftol,xsec,npmax,nnwt,nfix,npsfix
      common /vffix/ ffix(nsc)
      common /varnu/ dnu
C
      DIMENSION A(NSC),E(NSC)
C
      F=4*ANU*DT
C
c     D in /MEMO2/ holds the factors K^2 divided by 4
      DO 10 J=1,NSC
      ffix(j)=c(j)-4d0*anu*a(j)*d(j)
      C(J)=C(J)*DT-F*A(J)*D(J)
   10 G(J)=G(J)*DT-F*E(J)*D(J)-4d0*dt*a(j)*d(j)*dnu
c  ffix is the vector field not scaled with dt, for later
c  use in NWT and PROJECT
      call akill(ffix)
C
      IF(IFORCE.EQ.0) GO TO 15
C
C ******************************************************
C
C      FORCING TERMS  F(0,0,1)=-0.375, F(1,0,0)=0.25
C
      J1=(NS+1)*(NS+1)+1
C
      C(J1)=C(J1)-0.375*DT
      C(2)=C(2)+0.25*DT
C
C ******************************************************
C
   15 GO TO (20,40,70,100),ID
C
   20 DO 30 J=1,NSC
      A(J)=A(J)+0.5*C(J)
      B(J)=C(J)
      E(J)=E(J)+0.5*G(J)
   30 H(J)=G(J)
C
      CALL AKILL(A)
      CALL AKILL(E)
C
      RETURN
C
   40 S1=1-SQRT(0.5)
      S2=2-SQRT(2.0)
      S3=3*SQRT(0.5)-2
      DO 50 J=1,NSC
      A(J)=A(J)+S1*(C(J)-B(J))
   50 E(J)=E(J)+S1*(G(J)-H(J))
      DO 60 J=1,NSC
      B(J)=S2*C(J)+S3*B(J)
   60 H(J)=S2*G(J)+S3*H(J)
C
      CALL AKILL(A)
      CALL AKILL(E)
C
      RETURN
C
   70 S1=1+SQRT(0.5)
      S2=2+SQRT(2.0)
      S3=3*SQRT(0.5)+2
      DO 80 J=1,NSC
      A(J)=A(J)+S1*(C(J)-B(J))
   80 E(J)=E(J)+S1*(G(J)-H(J))
      DO 90 J=1,NSC
      B(J)=S2*C(J)-S3*B(J)
   90 H(J)=S2*G(J)-S3*H(J)
C
      CALL AKILL(A)
      CALL AKILL(E)
C
      RETURN
C
  100 S1=1./6
      DO 110 J=1,NSC
      A(J)=A(J)+S1*(C(J)-2*B(J))
  110 E(J)=E(J)+S1*(G(J)-2*H(J))
C
      CALL AKILL(A)
      CALL AKILL(E)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE MOMENT(A,ens,IFILE)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,NSP1=NS+1,NSC=(NS+1)**3)
      PARAMETER(N2=NSP1**2)
C
      COMMON /MEMO1/B(0:NS,0:NS,0:NS)
      COMMON /MEMO2/C(NSC)
      COMMON /TIME/TIME,ISTEP
C
ccccccccccc
c      common /section/ secxtol,secftol,xsec,fmax,nnwt,nfix,npsfix
ccccccccccc
      DIMENSION A(NSC),BB(NSC)
C
      EQUIVALENCE (B(0,0,0),BB(1))
C
      DO 30 J=1,NSC
   30 BB(J)=A(J)*A(J)
C
   40 DO 50 L=0,NS-1
      DO 50 K=L+1,NS
   50 B(0,K,L)=0.5*B(0,K,L)
C
      DO 60 K=0,NS
      DO 60 J=0,NS
   60 B(J,K,K)=0.5*B(J,K,K)
C
      ENE=0
      ENS=0
C
      DO 80 J=1,NSC
      ENE=ENE+BB(J)/C(J)
   80 ENS=ENS+BB(J)
C
      ENE=6*ENE
      ENS=24*ENS
C
c      WRITE(IFILE,*) ISTEP,TIME,ENE,ENS
c      WRITE(IFILE,*) ISTEP,TIME,ENE,ENS,a(npsfix)
c      IF(IFILE.EQ.15) WRITE(6,1000) ENE,ENS,TIME
c      IF(IFILE.EQ.16) WRITE(6,2000) ENE,ENS,TIME
c      IF(IFILE.EQ.17) WRITE(6,3000) ENE,ENS,TIME
 1000 FORMAT(8X,'ENERGY=',1PE14.7,2X,'ENSTROPHY=',E14.7,
     c        ' FOR RVF AT  TIME=',E14.7)
 2000 FORMAT(8X,'ENERGY=',1PE14.7,2X,'ENSTROPHY=',E14.7,
     c        ' FOR DVF AT  TIME=',E14.7)
 3000 FORMAT(8X,'ENERGY=',1PE14.7,2X,'ENSTROPHY=',E14.7,
     c        ' FOR DIFFER  AT  TIME=',E14.7)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE WNUM
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,NSC=(NS+1)**3,N1=(NS+1)**2)
      COMMON /MEMO2/A(0:NS,0:NS,0:NS)
      COMMON /CONST6/B(N1),C(0:NS)
C
      DIMENSION AA(0:NS,N1)
C
      EQUIVALENCE (A(0,0,0),AA(0,1))
C
      DO 10 J=0,NS
      DO 10 K=1,N1
   10 AA(J,K)=B(K)+C(J)
C
      DO 20 L=0,NS
      DO 20 K=0,L
      DO 20 J=0,NS
   20 A(J,K,L)=A(J,K,L)+J
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE AKILL(A)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND6=N/6)
      PARAMETER(NS=ND6,NSC=(NS+1)**3)
C
      DIMENSION A(0:NS,0:NS,0:NS)
C
      DO 10 L=0,NS
      DO 10 K=0,L
   10 A(ND6,K,L)=0
C
      DO 20 K=0,NS
      DO 20 J=0,NS
   20 A(J,K,ND6)=0
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE ORTH(A,B,ND6)
C
C  *****************************************
C  *   B(0:ND6,0:ND6,0:ND6)  WORKING AREA  *
C  *****************************************
C
c
      implicit double precision (a-h,o-z)
c
      DIMENSION A(0:ND6,0:ND6,0:ND6),B(0:ND6,0:ND6,0:ND6)
C
      S=1./3.
      ND6P1C=(ND6+1)**3
      ND6M1=ND6-1
C
      CALL ZERO(B,ND6P1C)
      CALL UP(A,ND6,2)
C
      DO 10 J=1,ND6
      DO 10 K=1,ND6
      DO 10 L=1,K
   10 A(J,K,L)=J*A(J,K,L)
C
      DO 20 J=0,ND6M1
      J2P1=J*2+1
      DO 20 K=0,ND6M1
      DO 20 L=K,ND6M1
   20 A(J,K,L+1)=J2P1*A(J,K,L+1)
C
      DO 30 J=1,ND6
      DO 30 K=1,J
      DO 30 L=1,K
      B(J,K,L)=S*(2*A(J,K,L)+A(K,J,L)-A(L,J,K))/J
      B(K,J,L)=S*(2*A(K,J,L)+A(L,J,K)+A(J,K,L))/K
   30 B(L,J,K)=S*(2*A(L,J,K)-A(J,K,L)+A(K,J,L))/L
C
      DO 40 J=0,ND6M1
      J2P1=J*2+1
      DO 40 K=0,J
      K2P1=K*2+1
      DO 40 L=0,K
      L2P1=L*2+1
      B(J,L,K+1)=S*(2*A(J,L,K+1)-A(K,L,J+1)-A(L,K,J+1))/J2P1
      B(K,L,J+1)=S*(2*A(K,L,J+1)-A(L,K,J+1)-A(J,L,K+1))/K2P1
   40 B(L,K,J+1)=S*(2*A(L,K,J+1)-A(J,L,K+1)-A(K,L,J+1))/L2P1
C
      DO 50 K=1,ND6
      DO 50 L=1,K
   50 B(0,K,L)=A(0,K,L)
C
      CALL DOWN(B,ND6,2)
      CALL COPY(B,A,ND6P1C)
C
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE UP(A,N,ID)
c
      implicit double precision (a-h,o-z)
c
C
C ***********************************************
C
C  ID = 1   FOR VELOCITY
C                         A(J,K,L-1)   (L > K+1)
C        A(J,K,L)  <---   0            (L = K+1)
C                         A(J,K,L)     (L < K+1)
C
C  ID = 2   FOR VORTICITY
C                         A(J,K,L-1)  (L > K)
C        A(J,K,L)  <---   0           (L = K)
C                         A(J,K,L)    (L < K)
C
C ***********************************************
      DIMENSION A(0:N,0:N,0:N)
      DO 10 L=N,1,-1
      LM1=L-1
      DO 10 K=0,LM1
      DO 10 J=0,N
   10 A(J,K,L)=A(J,K,LM1)
      GO TO (20,40),ID
   20 DO 30 J=0,N
      DO 30 K=0,N-1
   30 A(J,K,K+1)=0
      RETURN
   40 DO 50 J=0,N
      DO 50 K=0,N
   50 A(J,K,K)=0
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE DOWN(A,N,ID)
c
      implicit double precision (a-h,o-z)
c
C
C ***********************************************
C
C  ID = 1   FOR VELOCITY
C                         0            (L =N, K < N)
C        A(J,K,L)  <---   A(J,K,L+1)   (L > K)
C                         A(J,K,L)     (L =< K)
C
C  ID = 2   FOR VORTICITY
C                         0            (L =N, K < N)
C        A(J,K,L)  <---   A(J,K,L+1)   (L >= K)
C                         A(J,K,L)     (L < K)
C
C ***********************************************
      DIMENSION A(0:N,0:N,0:N)
      GO TO (10,30), ID
   10 DO 20 L=1,N-1
      LP1=L+1
      DO 20 K=0,L-1
      DO 20 J=0,N
   20 A(J,K,L)=A(J,K,LP1)
      GO TO 50
   30 DO 40 L=0,N-1
      LP1=L+1
      DO 40 K=0,L
      DO 40 J=0,N
   40 A(J,K,L)=A(J,K,LP1)
   50 DO 60 J=0,N
      DO 60 K=0,N-1
   60 A(J,K,N)=0
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE WFR(A)
C
C     FOURIER TRANSFORMATION OF VORTICITY FROM FOURIER TO REAL SPACE
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND6M1=ND6-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1)
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      COMMON /WORK1/B(0:ND4,0:ND4)
      DIMENSION A(0:ND4,0:ND4,0:ND4)
      CALL ZERO(B,(ND4+1)**2)
      DO 40 J=0,ND6M1
      DO 40 K=0,ND6M1
      B(J,K)=0
      S=-1
      DO 20 L=0,K
      S=-S
   20 B(J,K)=B(J,K)+S*A(J,L,K)
      IF(K.EQ.ND6M1) GO TO 40
      DO 30 L=K+1,ND6M1
      S=-S
   30 B(J,K)=B(J,K)+S*A(J,K,L)
   40 B(J,K)=B(J,K)*2
      CALL FTQSA(B)
      DO 55 L=0,ND6
      DO 55 K=0,L
      DO 50 J=ND6,1,-1
   50 A(J,K,L)=A(J,K,L)+A(J-1,K,L)
   55 A(0,K,L)=2*A(0,K,L)
      DO 80 J=0,ND6
      A(J,ND6,ND6)=A(J,ND6M1,ND6M1)
      DO 70 K=ND6M1,1,-1
      DO 60 L=ND6,K+1,-1
   60 A(J,K,L)=A(J,K,L)+A(J,K,L-1)+A(J,K-1,L)+A(J,K-1,L-1)
   70 A(J,K,K)=A(J,K,K)+2*A(J,K-1,K)+A(J,K-1,K-1)
      DO 80 L=ND6,0,-1
   80 A(J,0,L)=0
      CALL UNIFY2(A,ND4,2)
      CALL HELSAA(A,1)
      DO 90 J=0,ND4
      DO 90 K=0,ND4
      DO 90 L=0,K
      TEMP1=A(J,K,L)-A(J,L,K)
      TEMP2=(A(J,K,L)+A(J,L,K))*FCOS2I(J)*FCOS2I(K)*FCOS2I(L)
      A(J,K,L)=0.5*(TEMP1+TEMP2)
   90 A(J,L,K)=0.5*(-TEMP1+TEMP2)
      DO 110 J=0,ND4
      DO 100 K=0,ND4M1
      A(J,K,ND4)=A(J,K,ND4)+B(J,K)
  100 A(J,ND4,K)=A(J,ND4,K)+B(J,K)
  110 A(J,ND4,ND4)=A(J,ND4,ND4)+B(J,ND4)
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE HELSAA(A,IFR)
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
C
C   **************************************************************
C
C                    ND4
C   A(P,Q,R) <---   SIGMA   A(J,K,L) EXP(4*PI*I*(P*J+Q*K+R*L)/N)
C                 J,K,L=-ND4
C
C     WHERE  A(J,K,L) = A(-J,K,L) = -A(J,-K,L) = -A(J,K,-L)
C
C     INPUT AND OUTPUT: A(0:ND4,0:ND4,0:ND4)
C
C     IF  IFR = 1,  FOURIER SPACE ---> REAL SPACE
C     IF  IFR = 2,  REAL SPACE    ---> FOURIER SPACE
C
C   *************************************************************
C
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND6M1=ND6-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1)
      PARAMETER(NVW=ND8*(NN-1))
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      COMMON /WORK2/W1(0:ND4,0:ND4)
      COMMON /WORK4/W2(0:ND4,0:ND4)
      COMMON /WORK5/W3(0:ND4,0:ND4)
      COMMON /WORK6/W4(0:ND4M1),W5(0:ND8M1),W6(0:ND8M1)
c      COMMON /WORK7/VW(NVW),IVW(NVW)
      COMMON /WORK7/ip(0:nd8)
      common /work8/wwww(0:nd8*5/2-1)
c
      dimension abc(2*nd8)
      DIMENSION A(0:ND4,0:ND4,0:ND4)
      GO TO (1,2), IFR
    1 MAX1=ND6
      MAX2=ND4
      GO TO 5
    2 MAX1=ND4
      MAX2=ND6
      DO 3 K=0,ND4
      DO 3 L=0,ND4
    3 A(ND4,K,L)=2*A(ND4,K,L)
    5 DO 20 K=0,MAX1
      DO 20 L=0,MAX1
      S1=0
      S2=0
      DO 10 J=0,ND8M1
      J2=J*2
      S1=S1+A(J2,K,L)
   10 S2=S2+A(J2+1,K,L)
      W1(K,L)=2*(S1+S2)-A(0,K,L)+A(ND4,K,L)
   20 W2(K,L)=2*(S1-S2)-A(0,K,L)+A(ND4,K,L)
      DO 70 K=0,MAX1
      DO 70 L=0,MAX1
      W4(0)=A(0,K,L)
      W4(ND8)=A(ND4,K,L)
      DO 30 J=1,ND8M1
      J2=J*2
      S1=A(J2,K,L)
      S2=A(J2+1,K,L)-A(J2-1,K,L)
      W4(J)=S1+S2
   30 W4(ND4-J)=S1-S2
      W3(K,L)=0
      DO 40 J=0,ND8M1
      J2=J*2
   40 W3(K,L)=W3(K,L)+W4(J2)-W4(J2+1)
      W5(0)=W4(0)
      W6(0)=W4(ND4M1)-W4(1)
      DO 50 J=1,ND8M1
      J2=J*2
      W5(J)=W4(J2)
   50 W6(J)=W4(J2-1)-W4(J2+1)
c-------------------------------------------
c      CALL DVCFT1(W5,W6,ND8,-1,1,VW,IVW,ICON)
c
      call two2one(w5,w6,abc,nd8)
      call cdft(2*nd8,1,abc,ip,wwww)
      call one2two(abc,w5,w6,nd8)
c-------------------------------------------
      DO 60 J=1,ND8M1
      S1=F1(J)*W5(J)+F2(J)*W5(ND8-J)
      S2=F3(J)*W6(J)+F4(J)*W6(ND8-J)
      A(J,K,L)=S1-S2
   60 A(ND4-J,K,L)=S1+S2
      A(0,K,L)=W1(K,L)
      A(ND4,K,L)=W2(K,L)
   70 A(ND8,K,L)=W3(K,L)
      DO 170 J=0,MAX2
      DO 170 L=0,MAX1
      W4(0)=2*A(J,1,L)
      W4(ND8)=-2*A(J,ND4M1,L)
      DO 130 K=1,ND8M1
      K2=K*2
      S1=A(J,K2,L)
      S2=A(J,K2+1,L)-A(J,K2-1,L)
      W4(K)=S1+S2
  130 W4(ND4-K)=-S1+S2
      W3(J,L)=0
      DO 140 K=0,ND8M1
      K2=K*2
  140 W3(J,L)=W3(J,L)+W4(K2)-W4(K2+1)
      W3(J,L)=-0.5*W3(J,L)
      W5(0)=W4(0)
      W6(0)=W4(ND4M1)-W4(1)
      DO 150 K=1,ND8M1
      K2=K*2
      W5(K)=W4(K2)
  150 W6(K)=W4(K2-1)-W4(K2+1)
c-------------------------------------------
c      CALL DVCFT1(W5,W6,ND8,-1,1,VW,IVW,ICON)
c
      call two2one(w5,w6,abc,nd8)
      call cdft(2*nd8,1,abc,ip,wwww)
      call one2two(abc,w5,w6,nd8)
c-------------------------------------------
      DO 160 K=1,ND8M1
      S1=-F1(K)*W6(K)+F2(K)*W6(ND8-K)
      S2=-F3(K)*W5(K)+F4(K)*W5(ND8-K)
      A(J,K,L)=S1+S2
  160 A(J,ND4-K,L)=-S1+S2
      A(J,0,L)=0
      A(J,ND4,L)=0
  170 A(J,ND8,L)=W3(J,L)
      DO 270 J=0,MAX2
      DO 270 K=0,MAX2
      W4(0)=2*A(J,K,1)
      W4(ND8)=-2*A(J,K,ND4M1)
      DO 230 L=1,ND8M1
      L2=L*2
      S1=A(J,K,L2)
      S2=A(J,K,L2+1)-A(J,K,L2-1)
      W4(L)=S1+S2
  230 W4(ND4-L)=-S1+S2
      W3(J,K)=0
      DO 240 L=0,ND8M1
      L2=L*2
  240 W3(J,K)=W3(J,K)+W4(L2)-W4(L2+1)
      W3(J,K)=-0.5*W3(J,K)
      W5(0)=W4(0)
      W6(0)=W4(ND4M1)-W4(1)
      DO 250 L=1,ND8M1
      L2=L*2
      W5(L)=W4(L2)
  250 W6(L)=W4(L2-1)-W4(L2+1)
c-------------------------------------------
c      CALL DVCFT1(W5,W6,ND8,-1,1,VW,IVW,ICON)
c
      call two2one(w5,w6,abc,nd8)
      call cdft(2*nd8,1,abc,ip,wwww)
      call one2two(abc,w5,w6,nd8)
c-------------------------------------------
      DO 260 L=1,ND8M1
      S1=-F1(L)*W6(L)+F2(L)*W6(ND8-L)
      S2=F4(L)*W5(ND8-L)-F3(L)*W5(L)
      A(J,K,L)=-(S1+S2)
  260 A(J,K,ND4-L)=-(-S1+S2)
      A(J,K,0)=0
      A(J,K,ND4)=0
  270 A(J,K,ND8)=-W3(J,K)
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTQSA(A)
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
C
C   **************************************************************
C
C               ND4-1
C A(P,Q) <--- I SIGMA   A(J,K) EXP(2*PI*I*(P*(2*J+1)+Q*(2*K+1))/N)
C              J,K=-ND4
C
C     WHERE  A(J,K) = A(-J-1,K) = -A(J,-K-1)
C
C     INPUT AND OUTPUT: A(0:ND4,0:ND4)
C                       A(ND4,*) = A(*,ND4) = 0  (FOR INPUT)
C
C   *************************************************************
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND6M1=ND6-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1)
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      COMMON /WORK3/AA(0:ND4)
      DIMENSION A(0:ND4,0:ND4)
      DO 20 J=0,ND4M1
      S=-1
      AA(J)=0
      DO 10 K=0,ND4M1
      S=-S
   10 AA(J)=AA(J)+S*A(J,K)
   20 AA(J)=2*AA(J)
      CALL FTQS(AA)
      DO 40 J=ND4,1,-1
      JM1=J-1
      DO 30 K=ND4,1,-1
      KM1=K-1
   30 A(J,K)=A(J,K)+A(J,KM1)+A(JM1,K)+A(JM1,KM1)
   40 A(J,0)=0
      DO 50 K=ND4,1,-1
   50 A(0,K)=2*(A(0,K)+A(0,K-1))
      A(0,0)=0
      DO 60 K=0,ND4
   60 A(ND4,K)=2*A(ND4,K)
      CALL FTESA(A)
      DO 80 J=0,ND4
      DO 70 K=0,ND4M1
   70 A(J,K)=A(J,K)*FCOS2I(J)*FCOS2I(K)
   80 A(J,ND4)=-AA(J)
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTESA(A)
c
      implicit double precision (a-h,o-z)
c
C
C   **************************************************************
C
C
C                  ND4-1
C   A(P,Q) <---   I SIGMA   A(J,K) EXP(4*PI*I*(P*J+Q*K)/N)
C                 J,K=-ND4
C
C      WHERE  A(J,K) = A(-J,K) = -A(J,-K)
C
C      INPUTAND OUTPUT: A(0:ND4,0:ND4),  A(*,ND4) = 0.
C
C   *************************************************************
C
      include 'n_power.for'
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND6M1=ND6-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1)
      PARAMETER(NVW=ND8*(NN-1))
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      COMMON /WORK2/W1(0:ND4),W2(0:ND4),W3(0:ND4),
     c               W4(0:ND4M1),W5(0:ND8M1),W6(0:ND8M1)
c      COMMON /WORK7/VW(NVW),IVW(NVW)
      COMMON /WORK7/ip(0:nd8)
      common /work8/wwww(0:nd8*5/2-1)
c
      dimension abc(2*nd8)
      DIMENSION A(0:ND4,0:ND4)
      DO 20 K=0,ND4
      S1=0
      S2=0
      DO 10 J=0,ND8M1
      J2=J*2
      S1=S1+A(J2,K)
   10 S2=S2+A(J2+1,K)
      W1(K)=2*(S1+S2)-A(0,K)+A(ND4,K)
   20 W2(K)=2*(S1-S2)-A(0,K)+A(ND4,K)
      DO 70 K=0,ND4
      W4(0)=A(0,K)
      W4(ND8)=A(ND4,K)
      DO 30 J=1,ND8M1
      J2=J*2
      S1=A(J2,K)
      S2=A(J2+1,K)-A(J2-1,K)
      W4(J)=S1+S2
   30 W4(ND4-J)=S1-S2
      W3(K)=0
      DO 40 J=0,ND8M1
      J2=J*2
   40 W3(K)=W3(K)+W4(J2)-W4(J2+1)
      W5(0)=W4(0)
      W6(0)=W4(ND4M1)-W4(1)
      DO 50 J=1,ND8M1
      J2=J*2
      W5(J)=W4(J2)
   50 W6(J)=W4(J2-1)-W4(J2+1)
c-------------------------------------------
c      CALL DVCFT1(W5,W6,ND8,-1,1,VW,IVW,ICON)
c
      call two2one(w5,w6,abc,nd8)
      call cdft(2*nd8,1,abc,ip,wwww)
      call one2two(abc,w5,w6,nd8)
c-------------------------------------------
      DO 60 J=1,ND8M1
      S1=F1(J)*W5(J)+F2(J)*W5(ND8-J)
      S2=F3(J)*W6(J)+F4(J)*W6(ND8-J)
      A(J,K)=S1-S2
   60 A(ND4-J,K)=S1+S2
      A(0,K)=W1(K)
      A(ND4,K)=W2(K)
   70 A(ND8,K)=W3(K)
      DO 170 J=0,ND4
      W4(0)=2*A(J,1)
      W4(ND8)=-2*A(J,ND4M1)
      DO 130 K=1,ND8M1
      K2=K*2
      S1=A(J,K2)
      S2=A(J,K2+1)-A(J,K2-1)
      W4(K)=S1+S2
  130 W4(ND4-K)=-S1+S2
      W3(J)=0
      DO 140 K=0,ND8M1
      K2=K*2
  140 W3(J)=W3(J)+W4(K2)-W4(K2+1)
      W3(J)=-0.5*W3(J)
      W5(0)=W4(0)
      W6(0)=W4(ND4M1)-W4(1)
      DO 150 K=1,ND8M1
      K2=K*2
      W5(K)=W4(K2)
  150 W6(K)=W4(K2-1)-W4(K2+1)
c-------------------------------------------
c      CALL DVCFT1(W5,W6,ND8,-1,1,VW,IVW,ICON)
c
      call two2one(w5,w6,abc,nd8)
      call cdft(2*nd8,1,abc,ip,wwww)
      call one2two(abc,w5,w6,nd8)
c-------------------------------------------
      DO 160 K=1,ND8M1
      S1=-F1(K)*W6(K)+F2(K)*W6(ND8-K)
      S2=F4(K)*W5(ND8-K)-F3(K)*W5(K)
      A(J,K)=S1+S2
  160 A(J,ND4-K)=-S1+S2
      A(J,0)=0
      A(J,ND4)=0
  170 A(J,ND8)=W3(J)
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE FTQS(A)
c
      implicit double precision (a-h,o-z)
c
C
C   **************************************************************
C
C              ND4-1
C   A(P) <---  SIGMA  A(J) EXP(2*PI*I*P*(2*J+1)/N)
C              J=-ND4
C
C     WHERE  A(J) = A(-J-1)
C
C     INPUT AND OUTPUT: A(0:ND4)
C                       A(ND4) = 0  (FOR INPUT)
C
C   **************************************************************
C
      include 'n_power.for'
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6,ND8=N/8)
      PARAMETER(ND4M1=ND4-1,ND6M1=ND6-1,ND8M1=ND8-1)
      PARAMETER(NL=ND4+1)
      COMMON /FAC/FCOSI(0:NL),FCOS2I(0:NL),FSINI(ND8),
     c         F1(ND8),F2(ND8),F3(ND8),F4(ND8)
      DIMENSION A(0:ND4)
      A(ND4)=A(ND4-1)
      DO 10 J=ND4-1,1,-1
   10 A(J)=A(J)+A(J-1)
      A(0)=2*A(0)
      CALL FTES(A,2)
      DO 20 J=0,ND4
   20 A(J)=A(J)*FCOS2I(J)
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE UNIFY2(A,ND4,ID)
c
      implicit double precision (a-h,o-z)
c
C
C     IF  ID = 1  (FOR VELOCITY TYPE),
C                 0.5*(B(J,K,L)+B(J,L,K))  (K >= L)
C     A(J,K,L) =
C                 0.5*(B(J,K,L)-B(J,L,K))  (K < L)
C
C        --->
C
C     A(J,K,L) = B(J,K,L)
C
C     IF  ID = 2 (FOR VORTICITY TYPE),
C                 0.5*(B(J,K,L)-B(J,L,K))  (K > L)
C     A(J,K,L) =
C                 0.5*(B(J,K,L)+B(J,L,K))  (K =< L)
C
C        --->
C
C     A(J,K,L) = B(J,K,L)
C
      DIMENSION A(0:ND4,0:ND4,0:ND4)
      GO TO (10,30),ID
   10 DO 20 J=0,ND4
      DO 20 K=1,ND4
      DO 20 L=0,K-1
      A(J,K,L)=A(J,K,L)-A(J,L,K)
   20 A(J,L,K)=2*A(J,L,K)+A(J,K,L)
      RETURN
   30 DO 40 J=0,ND4
      DO 40 K=1,ND4
      DO 40 L=0,K-1
      A(J,K,L)=A(J,K,L)+A(J,L,K)
   40 A(J,L,K)=2*A(J,L,K)-A(J,K,L)
      RETURN
      END
*PROCESS DOUBLE
C
      SUBROUTINE HELCTY(A,B,C,IFILE)
C
c
      implicit double precision (a-h,o-z)
      include 'n_power.for'
c
      PARAMETER(N=2**NN,ND4=N/4,ND6=N/6)
      PARAMETER(NS=ND6,NL=ND4+1,NSC=(NS+1)**3)
C
      COMMON /TIME/TIME,ISTEP
      COMMON /CONST/ANU,DT
C
      DIMENSION A(0:ND4,0:ND4,0:ND4),B(0:NL,0:NL,0:NL)
      DIMENSION C(0:NS,0:NS,0:NS)
C
      F=3./ND4**3
      CALL ZERO(A,NL**3)
C
      DO 10 J=0,NS
      DO 10 K=0,NS
      DO 10 L=0,NS
   10 A(J,K,L)=C(J,K,L)*F
C
      CALL WFR(A)
C
      H=0.0
C
      DO 20 J=0,ND4
      WWW=1
      IF(J.EQ.0.OR.J.EQ.ND4) WWW=0.5
      DO 20 K=0,ND4
      WW=WWW
      IF(K.EQ.0.OR.K.EQ.ND4) WW=WW*0.5
      DO 20 L=0,ND4
      W=WW
      IF(L.EQ.0.OR.L.EQ.ND4) W=W*0.5
   20 H=H+W*A(J,K,L)*B(J,K,L)
C
      TIMEMD=TIME-DT
      ISTEPM=ISTEP-1
      WRITE(IFILE,*) ISTEPM,TIMEMD,H
      IF(IFILE.EQ.18) WRITE(6,1000) H,TIMEMD
      IF(IFILE.EQ.19) WRITE(6,2000) H,TIMEMD
      IF(IFILE.EQ.20) WRITE(6,3000) H,TIMEMD
C
      RETURN
 1000 FORMAT(31X,'HELICITY= ',1PE14.7,1X,'FOR RVF AT  TIME=',1PE14.7)
 2000 FORMAT(31X,'HELICITY= ',1PE14.7,1X,'FOR DVF AT  TIME=',1PE14.7)
 3000 FORMAT(31X,'HELICITY= ',1PE14.7,1X,'FOR DIFFER  AT  TIME=',
     c        1PE14.7)
C
      END
c===========================================================
c     The following routines are added by G.K. 2002/12/26
c-----------------------------------------------------------
      subroutine two2one(ar,ai,b,n)
      implicit real*8(a-h,o-z)
      dimension ar(n),ai(n),b(2*n)
c
      do 10 i=1,n
        b(2*i-1)=ar(i)
        b(2*i)  =ai(i)
10    continue
c
      return
      end
c
c-----------------------------------------------------------
      subroutine one2two(b,ar,ai,n)
      implicit real*8(a-h,o-z)
      dimension ar(n),ai(n),b(2*n)
c
      do 10 i=1,n
        ar(i)=b(2*i-1)
        ai(i)=b(2*i)
10    continue
c
      return
      end
c===========================================================
c Fast Fourier/Cosine/Sine Transform
c     dimension   :one
c     data length :power of 2
c     decimation  :frequency
c     radix       :4, 2
c     data        :inplace
c     table       :use
c subroutines
c     cdft: Complex Discrete Fourier Transform
c     rdft: Real Discrete Fourier Transform
c     ddct: Discrete Cosine Transform
c     ddst: Discrete Sine Transform
c     dfct: Cosine Transform of RDFT (Real Symmetric DFT)
c     dfst: Sine Transform of RDFT (Real Anti-symmetric DFT)
c
c
c -------- Complex DFT (Discrete Fourier Transform) --------
c     [definition]
c         <case1>
c             X(k) = sum_j=0^n-1 x(j)*exp(2*pi*i*j*k/n), 0<=k<n
c         <case2>
c             X(k) = sum_j=0^n-1 x(j)*exp(-2*pi*i*j*k/n), 0<=k<n
c         (notes: sum_j=0^n-1 is a summation from j=0 to n-1)
c     [usage]
c         <case1>
c             ip(0) = 0  c first time only
c             call cdft(2*n, 1, a, ip, w)
c         <case2>
c             ip(0) = 0  c first time only
c             call cdft(2*n, -1, a, ip, w)
c     [parameters]
c         2*n          :data length (integer)
c                       n >= 1, n = power of 2
c         a(0:2*n-1)   :input/output data (real*8)
c                       input data
c                           a(2*j) = Re(x(j)), 
c                           a(2*j+1) = Im(x(j)), 0<=j<n
c                       output data
c                           a(2*k) = Re(X(k)), 
c                           a(2*k+1) = Im(X(k)), 0<=k<n
c         ip(0:*)      :work area for bit reversal (integer)
c                       length of ip >= 2+sqrt(n)
c                       strictly, 
c                       length of ip >= 
c                           2+2**(int(log(n+0.5)/log(2.0))/2).
c                       ip(0),ip(1) are pointers of the cos/sin table.
c         w(0:n/2-1)   :cos/sin table (real*8)
c                       w(),ip() are initialized if ip(0) = 0.
c     [remark]
c         Inverse of 
c             call cdft(2*n, -1, a, ip, w)
c         is 
c             call cdft(2*n, 1, a, ip, w)
c             do j = 0, 2 * n - 1
c                 a(j) = a(j) / n
c             end do
c         .
c
c
c -------- Real DFT / Inverse of Real DFT --------
c     [definition]
c         <case1> RDFT
c             R(k) = sum_j=0^n-1 a(j)*cos(2*pi*j*k/n), 0<=k<=n/2
c             I(k) = sum_j=0^n-1 a(j)*sin(2*pi*j*k/n), 0<k<n/2
c         <case2> IRDFT (excluding scale)
c             a(k) = (R(0) + R(n/2)*cos(pi*k))/2 + 
c                    sum_j=1^n/2-1 R(j)*cos(2*pi*j*k/n) + 
c                    sum_j=1^n/2-1 I(j)*sin(2*pi*j*k/n), 0<=k<n
c     [usage]
c         <case1>
c             ip(0) = 0  c first time only
c             call rdft(n, 1, a, ip, w)
c         <case2>
c             ip(0) = 0  c first time only
c             call rdft(n, -1, a, ip, w)
c     [parameters]
c         n            :data length (integer)
c                       n >= 2, n = power of 2
c         a(0:n-1)     :input/output data (real*8)
c                       <case1>
c                           output data
c                               a(2*k) = R(k), 0<=k<n/2
c                               a(2*k+1) = I(k), 0<k<n/2
c                               a(1) = R(n/2)
c                       <case2>
c                           input data
c                               a(2*j) = R(j), 0<=j<n/2
c                               a(2*j+1) = I(j), 0<j<n/2
c                               a(1) = R(n/2)
c         ip(0:*)      :work area for bit reversal (integer)
c                       length of ip >= 2+sqrt(n/2)
c                       strictly, 
c                       length of ip >= 
c                           2+2**(int(log(n/2+0.5)/log(2.0))/2).
c                       ip(0),ip(1) are pointers of the cos/sin table.
c         w(0:n/2-1)   :cos/sin table (real*8)
c                       w(),ip() are initialized if ip(0) = 0.
c     [remark]
c         Inverse of 
c             call rdft(n, 1, a, ip, w)
c         is 
c             call rdft(n, -1, a, ip, w)
c             do j = 0, n - 1
c                 a(j) = a(j) * 2 / n
c             end do
c         .
c
c
c -------- DCT (Discrete Cosine Transform) / Inverse of DCT --------
c     [definition]
c         <case1> IDCT (excluding scale)
c             C(k) = sum_j=0^n-1 a(j)*cos(pi*j*(k+1/2)/n), 0<=k<n
c         <case2> DCT
c             C(k) = sum_j=0^n-1 a(j)*cos(pi*(j+1/2)*k/n), 0<=k<n
c     [usage]
c         <case1>
c             ip(0) = 0  c first time only
c             call ddct(n, 1, a, ip, w)
c         <case2>
c             ip(0) = 0  c first time only
c             call ddct(n, -1, a, ip, w)
c     [parameters]
c         n            :data length (integer)
c                       n >= 2, n = power of 2
c         a(0:n-1)     :input/output data (real*8)
c                       output data
c                           a(k) = C(k), 0<=k<n
c         ip(0:*)      :work area for bit reversal (integer)
c                       length of ip >= 2+sqrt(n/2)
c                       strictly, 
c                       length of ip >= 
c                           2+2**(int(log(n/2+0.5)/log(2.0))/2).
c                       ip(0),ip(1) are pointers of the cos/sin table.
c         w(0:n*5/4-1) :cos/sin table (real*8)
c                       w(),ip() are initialized if ip(0) = 0.
c     [remark]
c         Inverse of 
c             call ddct(n, -1, a, ip, w)
c         is 
c             a(0) = a(0) / 2
c             call ddct(n, 1, a, ip, w)
c             do j = 0, n - 1
c                 a(j) = a(j) * 2 / n
c             end do
c         .
c
c
c -------- DST (Discrete Sine Transform) / Inverse of DST --------
c     [definition]
c         <case1> IDST (excluding scale)
c             S(k) = sum_j=1^n A(j)*sin(pi*j*(k+1/2)/n), 0<=k<n
c         <case2> DST
c             S(k) = sum_j=0^n-1 a(j)*sin(pi*(j+1/2)*k/n), 0<k<=n
c     [usage]
c         <case1>
c             ip(0) = 0  c first time only
c             call ddst(n, 1, a, ip, w)
c         <case2>
c             ip(0) = 0  c first time only
c             call ddst(n, -1, a, ip, w)
c     [parameters]
c         n            :data length (integer)
c                       n >= 2, n = power of 2
c         a(0:n-1)     :input/output data (real*8)
c                       <case1>
c                           input data
c                               a(j) = A(j), 0<j<n
c                               a(0) = A(n)
c                           output data
c                               a(k) = S(k), 0<=k<n
c                       <case2>
c                           output data
c                               a(k) = S(k), 0<k<n
c                               a(0) = S(n)
c         ip(0:*)      :work area for bit reversal (integer)
c                       length of ip >= 2+sqrt(n/2)
c                       strictly, 
c                       length of ip >= 
c                           2+2**(int(log(n/2+0.5)/log(2.0))/2).
c                       ip(0),ip(1) are pointers of the cos/sin table.
c         w(0:n*5/4-1) :cos/sin table (real*8)
c                       w(),ip() are initialized if ip(0) = 0.
c     [remark]
c         Inverse of 
c             call ddst(n, -1, a, ip, w)
c         is 
c             a(0) = a(0) / 2
c             call ddst(n, 1, a, ip, w)
c             do j = 0, n - 1
c                 a(j) = a(j) * 2 / n
c             end do
c         .
c
c
c -------- Cosine Transform of RDFT (Real Symmetric DFT) --------
c     [definition]
c         C(k) = sum_j=0^n a(j)*cos(pi*j*k/n), 0<=k<=n
c     [usage]
c         ip(0) = 0  c first time only
c         call dfct(n, a, t, ip, w)
c     [parameters]
c         n            :data length - 1 (integer)
c                       n >= 2, n = power of 2
c         a(0:n)       :input/output data (real*8)
c                       output data
c                           a(k) = C(k), 0<=k<=n
c         t(0:n/2)     :work area (real*8)
c         ip(0:*)      :work area for bit reversal (integer)
c                       length of ip >= 2+sqrt(n/4)
c                       strictly, 
c                       length of ip >= 
c                           2+2**(int(log(n/4+0.5)/log(2.0))/2).
c                       ip(0),ip(1) are pointers of the cos/sin table.
c         w(0:n*5/8-1) :cos/sin table (real*8)
c                       w(),ip() are initialized if ip(0) = 0.
c     [remark]
c         Inverse of 
c             a(0) = a(0) / 2
c             a(n) = a(n) / 2
c             call dfct(n, a, t, ip, w)
c         is 
c             a(0) = a(0) / 2
c             a(n) = a(n) / 2
c             call dfct(n, a, t, ip, w)
c             do j = 0, n
c                 a(j) = a(j) * 2 / n
c             end do
c         .
c
c
c -------- Sine Transform of RDFT (Real Anti-symmetric DFT) --------
c     [definition]
c         S(k) = sum_j=1^n-1 a(j)*sin(pi*j*k/n), 0<k<n
c     [usage]
c         ip(0) = 0  c first time only
c         call dfst(n, a, t, ip, w)
c     [parameters]
c         n            :data length + 1 (integer)
c                       n >= 2, n = power of 2
c         a(0:n-1)     :input/output data (real*8)
c                       output data
c                           a(k) = S(k), 0<k<n
c                       (a(0) is used for work area)
c         t(0:n/2-1)   :work area (real*8)
c         ip(0:*)      :work area for bit reversal (integer)
c                       length of ip >= 2+sqrt(n/4)
c                       strictly, 
c                       length of ip >= 
c                           2+2**(int(log(n/4+0.5)/log(2.0))/2).
c                       ip(0),ip(1) are pointers of the cos/sin table.
c         w(0:n*5/8-1) :cos/sin table (real*8)
c                       w(),ip() are initialized if ip(0) = 0.
c     [remark]
c         Inverse of 
c             call dfst(n, a, t, ip, w)
c         is 
c             call dfst(n, a, t, ip, w)
c             do j = 1, n - 1
c                 a(j) = a(j) * 2 / n
c             end do
c         .
c
c
c Appendix :
c     The cos/sin table is recalculated when the larger table required.
c     w() and ip() are compatible with all routines.
c
c
      subroutine cdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *)
      real*8 a(0 : n - 1), w(0 : *)
      if (n .gt. 4 * ip(0)) then
          call makewt(n / 4, ip, w)
      end if
      if (n .gt. 4) then
          if (isgn .ge. 0) then
              call bitrv2(n, ip(2), a)
              call cftfsub(n, a, w)
          else
              call bitrv2conj(n, ip(2), a)
              call cftbsub(n, a, w)
          end if
      else if (n .eq. 4) then
          call cftfsub(n, a, w)
      end if
      end
c
      subroutine rdft(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), nw, nc
      real*8 a(0 : n - 1), w(0 : *), xi
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 4 * nc) then
          nc = n / 4
          call makect(nc, ip, w(nw))
      end if
      if (isgn .ge. 0) then
          if (n .gt. 4) then
              call bitrv2(n, ip(2), a)
              call cftfsub(n, a, w)
              call rftfsub(n, a, nc, w(nw))
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
          xi = a(0) - a(1)
          a(0) = a(0) + a(1)
          a(1) = xi
      else
          a(1) = 0.5d0 * (a(0) - a(1))
          a(0) = a(0) - a(1)
          if (n .gt. 4) then
              call rftbsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
              call cftbsub(n, a, w)
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
      end if
      end
c
      subroutine ddct(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), j, nw, nc
      real*8 a(0 : n - 1), w(0 : *), xr
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. nc) then
          nc = n
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          xr = a(n - 1)
          do j = n - 2, 2, -2
              a(j + 1) = a(j) - a(j - 1)
              a(j) = a(j) + a(j - 1)
          end do
          a(1) = a(0) - xr
          a(0) = a(0) + xr
          if (n .gt. 4) then
              call rftbsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
              call cftbsub(n, a, w)
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
      end if
      call dctsub(n, a, nc, w(nw))
      if (isgn .ge. 0) then
          if (n .gt. 4) then
              call bitrv2(n, ip(2), a)
              call cftfsub(n, a, w)
              call rftfsub(n, a, nc, w(nw))
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
          xr = a(0) - a(1)
          a(0) = a(0) + a(1)
          do j = 2, n - 2, 2
              a(j - 1) = a(j) - a(j + 1)
              a(j) = a(j) + a(j + 1)
          end do
          a(n - 1) = xr
      end if
      end
c
      subroutine ddst(n, isgn, a, ip, w)
      integer n, isgn, ip(0 : *), j, nw, nc
      real*8 a(0 : n - 1), w(0 : *), xr
      nw = ip(0)
      if (n .gt. 4 * nw) then
          nw = n / 4
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. nc) then
          nc = n
          call makect(nc, ip, w(nw))
      end if
      if (isgn .lt. 0) then
          xr = a(n - 1)
          do j = n - 2, 2, -2
              a(j + 1) = -a(j) - a(j - 1)
              a(j) = a(j) - a(j - 1)
          end do
          a(1) = a(0) + xr
          a(0) = a(0) - xr
          if (n .gt. 4) then
              call rftbsub(n, a, nc, w(nw))
              call bitrv2(n, ip(2), a)
              call cftbsub(n, a, w)
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
      end if
      call dstsub(n, a, nc, w(nw))
      if (isgn .ge. 0) then
          if (n .gt. 4) then
              call bitrv2(n, ip(2), a)
              call cftfsub(n, a, w)
              call rftfsub(n, a, nc, w(nw))
          else if (n .eq. 4) then
              call cftfsub(n, a, w)
          end if
          xr = a(0) - a(1)
          a(0) = a(0) + a(1)
          do j = 2, n - 2, 2
              a(j - 1) = -a(j) - a(j + 1)
              a(j) = a(j) - a(j + 1)
          end do
          a(n - 1) = -xr
      end if
      end
c
      subroutine dfct(n, a, t, ip, w)
      integer n, ip(0 : *), j, k, l, m, mh, nw, nc
      real*8 a(0 : n), t(0 : n / 2), w(0 : *), xr, xi, yr, yi
      nw = ip(0)
      if (n .gt. 8 * nw) then
          nw = n / 8
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 2 * nc) then
          nc = n / 2
          call makect(nc, ip, w(nw))
      end if
      m = n / 2
      yi = a(m)
      xi = a(0) + a(n)
      a(0) = a(0) - a(n)
      t(0) = xi - yi
      t(m) = xi + yi
      if (n .gt. 2) then
          mh = m / 2
          do j = 1, mh - 1
              k = m - j
              xr = a(j) - a(n - j)
              xi = a(j) + a(n - j)
              yr = a(k) - a(n - k)
              yi = a(k) + a(n - k)
              a(j) = xr
              a(k) = yr
              t(j) = xi - yi
              t(k) = xi + yi
          end do
          t(mh) = a(mh) + a(n - mh)
          a(mh) = a(mh) - a(n - mh)
          call dctsub(m, a, nc, w(nw))
          if (m .gt. 4) then
              call bitrv2(m, ip(2), a)
              call cftfsub(m, a, w)
              call rftfsub(m, a, nc, w(nw))
          else if (m .eq. 4) then
              call cftfsub(m, a, w)
          end if
          a(n - 1) = a(0) - a(1)
          a(1) = a(0) + a(1)
          do j = m - 2, 2, -2
              a(2 * j + 1) = a(j) + a(j + 1)
              a(2 * j - 1) = a(j) - a(j + 1)
          end do
          l = 2
          m = mh
          do while (m .ge. 2)
              call dctsub(m, t, nc, w(nw))
              if (m .gt. 4) then
                  call bitrv2(m, ip(2), t)
                  call cftfsub(m, t, w)
                  call rftfsub(m, t, nc, w(nw))
              else if (m .eq. 4) then
                  call cftfsub(m, t, w)
              end if
              a(n - l) = t(0) - t(1)
              a(l) = t(0) + t(1)
              k = 0
              do j = 2, m - 2, 2
                  k = k + 4 * l
                  a(k - l) = t(j) - t(j + 1)
                  a(k + l) = t(j) + t(j + 1)
              end do
              l = 2 * l
              mh = m / 2
              do j = 0, mh - 1
                  k = m - j
                  t(j) = t(m + k) - t(m + j)
                  t(k) = t(m + k) + t(m + j)
              end do
              t(mh) = t(m + mh)
              m = mh
          end do
          a(l) = t(0)
          a(n) = t(2) - t(1)
          a(0) = t(2) + t(1)
      else
          a(1) = a(0)
          a(2) = t(0)
          a(0) = t(1)
      end if
      end
c
      subroutine dfst(n, a, t, ip, w)
      integer n, ip(0 : *), j, k, l, m, mh, nw, nc
      real*8 a(0 : n - 1), t(0 : n / 2 - 1), w(0 : *), xr, xi, yr, yi
      nw = ip(0)
      if (n .gt. 8 * nw) then
          nw = n / 8
          call makewt(nw, ip, w)
      end if
      nc = ip(1)
      if (n .gt. 2 * nc) then
          nc = n / 2
          call makect(nc, ip, w(nw))
      end if
      if (n .gt. 2) then
          m = n / 2
          mh = m / 2
          do j = 1, mh - 1
              k = m - j
              xr = a(j) + a(n - j)
              xi = a(j) - a(n - j)
              yr = a(k) + a(n - k)
              yi = a(k) - a(n - k)
              a(j) = xr
              a(k) = yr
              t(j) = xi + yi
              t(k) = xi - yi
          end do
          t(0) = a(mh) - a(n - mh)
          a(mh) = a(mh) + a(n - mh)
          a(0) = a(m)
          call dstsub(m, a, nc, w(nw))
          if (m .gt. 4) then
              call bitrv2(m, ip(2), a)
              call cftfsub(m, a, w)
              call rftfsub(m, a, nc, w(nw))
          else if (m .eq. 4) then
              call cftfsub(m, a, w)
          end if
          a(n - 1) = a(1) - a(0)
          a(1) = a(0) + a(1)
          do j = m - 2, 2, -2
              a(2 * j + 1) = a(j) - a(j + 1)
              a(2 * j - 1) = -a(j) - a(j + 1)
          end do
          l = 2
          m = mh
          do while (m .ge. 2)
              call dstsub(m, t, nc, w(nw))
              if (m .gt. 4) then
                  call bitrv2(m, ip(2), t)
                  call cftfsub(m, t, w)
                  call rftfsub(m, t, nc, w(nw))
              else if (m .eq. 4) then
                  call cftfsub(m, t, w)
              end if
              a(n - l) = t(1) - t(0)
              a(l) = t(0) + t(1)
              k = 0
              do j = 2, m - 2, 2
                  k = k + 4 * l
                  a(k - l) = -t(j) - t(j + 1)
                  a(k + l) = t(j) - t(j + 1)
              end do
              l = 2 * l
              mh = m / 2
              do j = 1, mh - 1
                  k = m - j
                  t(j) = t(m + k) + t(m + j)
                  t(k) = t(m + k) - t(m + j)
              end do
              t(0) = t(m + mh)
              m = mh
          end do
          a(l) = t(0)
      end if
      a(0) = 0
      end
c
c -------- initializing routines --------
c
      subroutine makewt(nw, ip, w)
      integer nw, ip(0 : *), j, nwh
      real*8 w(0 : nw - 1), delta, x, y
      ip(0) = nw
      ip(1) = 1
      if (nw .gt. 2) then
          nwh = nw / 2
          delta = atan(1.0d0) / nwh
          w(0) = 1
          w(1) = 0
          w(nwh) = cos(delta * nwh)
          w(nwh + 1) = w(nwh)
          if (nwh .gt. 2) then
              do j = 2, nwh - 2, 2
                  x = cos(delta * j)
                  y = sin(delta * j)
                  w(j) = x
                  w(j + 1) = y
                  w(nw - j) = y
                  w(nw - j + 1) = x
              end do
              call bitrv2(nw, ip(2), w)
          end if
      end if
      end
c
      subroutine makect(nc, ip, c)
      integer nc, ip(0 : *), j, nch
      real*8 c(0 : nc - 1), delta
      ip(1) = nc
      if (nc .gt. 1) then
          nch = nc / 2
          delta = atan(1.0d0) / nch
          c(0) = cos(delta * nch)
          c(nch) = 0.5d0 * c(0)
          do j = 1, nch - 1
              c(j) = 0.5d0 * cos(delta * j)
              c(nc - j) = 0.5d0 * sin(delta * j)
          end do
      end if
      end
c
c -------- child routines --------
c
      subroutine bitrv2(n, ip, a)
      integer n, ip(0 : *), j, j1, k, k1, l, m, m2
      real*8 a(0 : n - 1), xr, xi, yr, yi
      ip(0) = 0
      l = n
      m = 1
      do while (8 * m .lt. l)
          l = l / 2
          do j = 0, m - 1
              ip(m + j) = ip(j) + l
          end do
          m = m * 2
      end do
      m2 = 2 * m
      if (8 * m .eq. l) then
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 - m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              j1 = 2 * k + m2 + ip(k)
              k1 = j1 + m2
              xr = a(j1)
              xi = a(j1 + 1)
              yr = a(k1)
              yi = a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
          end do
      else
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + m2
                  xr = a(j1)
                  xi = a(j1 + 1)
                  yr = a(k1)
                  yi = a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
          end do
      end if
      end
c
      subroutine bitrv2conj(n, ip, a)
      integer n, ip(0 : *), j, j1, k, k1, l, m, m2
      real*8 a(0 : n - 1), xr, xi, yr, yi
      ip(0) = 0
      l = n
      m = 1
      do while (8 * m .lt. l)
          l = l / 2
          do j = 0, m - 1
              ip(m + j) = ip(j) + l
          end do
          m = m * 2
      end do
      m2 = 2 * m
      if (8 * m .eq. l) then
          do k = 0, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 - m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + 2 * m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 2 * k + ip(k)
              a(k1 + 1) = -a(k1 + 1)
              j1 = k1 + m2
              k1 = j1 + m2
              xr = a(j1)
              xi = -a(j1 + 1)
              yr = a(k1)
              yi = -a(k1 + 1)
              a(j1) = yr
              a(j1 + 1) = yi
              a(k1) = xr
              a(k1 + 1) = xi
              k1 = k1 + m2
              a(k1 + 1) = -a(k1 + 1)
          end do
      else
          a(1) = -a(1)
          a(m2 + 1) = -a(m2 + 1)
          do k = 1, m - 1
              do j = 0, k - 1
                  j1 = 2 * j + ip(k)
                  k1 = 2 * k + ip(j)
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
                  j1 = j1 + m2
                  k1 = k1 + m2
                  xr = a(j1)
                  xi = -a(j1 + 1)
                  yr = a(k1)
                  yi = -a(k1 + 1)
                  a(j1) = yr
                  a(j1 + 1) = yi
                  a(k1) = xr
                  a(k1 + 1) = xi
              end do
              k1 = 2 * k + ip(k)
              a(k1 + 1) = -a(k1 + 1)
              a(k1 + m2 + 1) = -a(k1 + m2 + 1)
          end do
      end if
      end
c
      subroutine cftfsub(n, a, w)
      integer n, j, j1, j2, j3, l
      real*8 a(0 : n - 1), w(0 : *)
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      l = 2
      if (n .gt. 8) then
          call cft1st(n, a, w)
          l = 8
          do while (4 * l .lt. n)
              call cftmdl(n, l, a, w)
              l = 4 * l
          end do
      end if
      if (4 * l .eq. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i - x2i
              a(j1) = x1r - x3i
              a(j1 + 1) = x1i + x3r
              a(j3) = x1r + x3i
              a(j3 + 1) = x1i - x3r
          end do
      else
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = a(j + 1) - a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = a(j + 1) + a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end
c
      subroutine cftbsub(n, a, w)
      integer n, j, j1, j2, j3, l
      real*8 a(0 : n - 1), w(0 : *)
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      l = 2
      if (n .gt. 8) then
          call cft1st(n, a, w)
          l = 8
          do while (4 * l .lt. n)
              call cftmdl(n, l, a, w)
              l = 4 * l
          end do
      end if
      if (4 * l .eq. n) then
          do j = 0, l - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = -a(j + 1) - a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = -a(j + 1) + a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i - x2i
              a(j2) = x0r - x2r
              a(j2 + 1) = x0i + x2i
              a(j1) = x1r - x3i
              a(j1 + 1) = x1i - x3r
              a(j3) = x1r + x3i
              a(j3 + 1) = x1i + x3r
          end do
      else
          do j = 0, l - 2, 2
              j1 = j + l
              x0r = a(j) - a(j1)
              x0i = -a(j + 1) + a(j1 + 1)
              a(j) = a(j) + a(j1)
              a(j + 1) = -a(j + 1) - a(j1 + 1)
              a(j1) = x0r
              a(j1 + 1) = x0i
          end do
      end if
      end
c
      subroutine cft1st(n, a, w)
      integer n, j, k1, k2
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      x0r = a(0) + a(2)
      x0i = a(1) + a(3)
      x1r = a(0) - a(2)
      x1i = a(1) - a(3)
      x2r = a(4) + a(6)
      x2i = a(5) + a(7)
      x3r = a(4) - a(6)
      x3i = a(5) - a(7)
      a(0) = x0r + x2r
      a(1) = x0i + x2i
      a(4) = x0r - x2r
      a(5) = x0i - x2i
      a(2) = x1r - x3i
      a(3) = x1i + x3r
      a(6) = x1r + x3i
      a(7) = x1i - x3r
      wk1r = w(2)
      x0r = a(8) + a(10)
      x0i = a(9) + a(11)
      x1r = a(8) - a(10)
      x1i = a(9) - a(11)
      x2r = a(12) + a(14)
      x2i = a(13) + a(15)
      x3r = a(12) - a(14)
      x3i = a(13) - a(15)
      a(8) = x0r + x2r
      a(9) = x0i + x2i
      a(12) = x2i - x0i
      a(13) = x0r - x2r
      x0r = x1r - x3i
      x0i = x1i + x3r
      a(10) = wk1r * (x0r - x0i)
      a(11) = wk1r * (x0r + x0i)
      x0r = x3i + x1r
      x0i = x3r - x1i
      a(14) = wk1r * (x0i - x0r)
      a(15) = wk1r * (x0i + x0r)
      k1 = 0
      do j = 16, n - 16, 16
          k1 = k1 + 2
          k2 = 2 * k1
          wk2r = w(k1)
          wk2i = w(k1 + 1)
          wk1r = w(k2)
          wk1i = w(k2 + 1)
          wk3r = wk1r - 2 * wk2i * wk1i
          wk3i = 2 * wk2i * wk1r - wk1i
          x0r = a(j) + a(j + 2)
          x0i = a(j + 1) + a(j + 3)
          x1r = a(j) - a(j + 2)
          x1i = a(j + 1) - a(j + 3)
          x2r = a(j + 4) + a(j + 6)
          x2i = a(j + 5) + a(j + 7)
          x3r = a(j + 4) - a(j + 6)
          x3i = a(j + 5) - a(j + 7)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          x0r = x0r - x2r
          x0i = x0i - x2i
          a(j + 4) = wk2r * x0r - wk2i * x0i
          a(j + 5) = wk2r * x0i + wk2i * x0r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j + 2) = wk1r * x0r - wk1i * x0i
          a(j + 3) = wk1r * x0i + wk1i * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j + 6) = wk3r * x0r - wk3i * x0i
          a(j + 7) = wk3r * x0i + wk3i * x0r
          wk1r = w(k2 + 2)
          wk1i = w(k2 + 3)
          wk3r = wk1r - 2 * wk2r * wk1i
          wk3i = 2 * wk2r * wk1r - wk1i
          x0r = a(j + 8) + a(j + 10)
          x0i = a(j + 9) + a(j + 11)
          x1r = a(j + 8) - a(j + 10)
          x1i = a(j + 9) - a(j + 11)
          x2r = a(j + 12) + a(j + 14)
          x2i = a(j + 13) + a(j + 15)
          x3r = a(j + 12) - a(j + 14)
          x3i = a(j + 13) - a(j + 15)
          a(j + 8) = x0r + x2r
          a(j + 9) = x0i + x2i
          x0r = x0r - x2r
          x0i = x0i - x2i
          a(j + 12) = -wk2i * x0r - wk2r * x0i
          a(j + 13) = -wk2i * x0i + wk2r * x0r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j + 10) = wk1r * x0r - wk1i * x0i
          a(j + 11) = wk1r * x0i + wk1i * x0r
          x0r = x1r + x3i
          x0i = x1i - x3r
          a(j + 14) = wk3r * x0r - wk3i * x0i
          a(j + 15) = wk3r * x0i + wk3i * x0r
      end do
      end
c
      subroutine cftmdl(n, l, a, w)
      integer n, l, j, j1, j2, j3, k, k1, k2, m, m2
      real*8 a(0 : n - 1), w(0 : *)
      real*8 wk1r, wk1i, wk2r, wk2i, wk3r, wk3i
      real*8 x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i
      m = 4 * l
      do j = 0, l - 2, 2
          j1 = j + l
          j2 = j1 + l
          j3 = j2 + l
          x0r = a(j) + a(j1)
          x0i = a(j + 1) + a(j1 + 1)
          x1r = a(j) - a(j1)
          x1i = a(j + 1) - a(j1 + 1)
          x2r = a(j2) + a(j3)
          x2i = a(j2 + 1) + a(j3 + 1)
          x3r = a(j2) - a(j3)
          x3i = a(j2 + 1) - a(j3 + 1)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j2) = x0r - x2r
          a(j2 + 1) = x0i - x2i
          a(j1) = x1r - x3i
          a(j1 + 1) = x1i + x3r
          a(j3) = x1r + x3i
          a(j3 + 1) = x1i - x3r
      end do
      wk1r = w(2)
      do j = m, l + m - 2, 2
          j1 = j + l
          j2 = j1 + l
          j3 = j2 + l
          x0r = a(j) + a(j1)
          x0i = a(j + 1) + a(j1 + 1)
          x1r = a(j) - a(j1)
          x1i = a(j + 1) - a(j1 + 1)
          x2r = a(j2) + a(j3)
          x2i = a(j2 + 1) + a(j3 + 1)
          x3r = a(j2) - a(j3)
          x3i = a(j2 + 1) - a(j3 + 1)
          a(j) = x0r + x2r
          a(j + 1) = x0i + x2i
          a(j2) = x2i - x0i
          a(j2 + 1) = x0r - x2r
          x0r = x1r - x3i
          x0i = x1i + x3r
          a(j1) = wk1r * (x0r - x0i)
          a(j1 + 1) = wk1r * (x0r + x0i)
          x0r = x3i + x1r
          x0i = x3r - x1i
          a(j3) = wk1r * (x0i - x0r)
          a(j3 + 1) = wk1r * (x0i + x0r)
      end do
      k1 = 0
      m2 = 2 * m
      do k = m2, n - m2, m2
          k1 = k1 + 2
          k2 = 2 * k1
          wk2r = w(k1)
          wk2i = w(k1 + 1)
          wk1r = w(k2)
          wk1i = w(k2 + 1)
          wk3r = wk1r - 2 * wk2i * wk1i
          wk3i = 2 * wk2i * wk1r - wk1i
          do j = k, l + k - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              x0r = x0r - x2r
              x0i = x0i - x2i
              a(j2) = wk2r * x0r - wk2i * x0i
              a(j2 + 1) = wk2r * x0i + wk2i * x0r
              x0r = x1r - x3i
              x0i = x1i + x3r
              a(j1) = wk1r * x0r - wk1i * x0i
              a(j1 + 1) = wk1r * x0i + wk1i * x0r
              x0r = x1r + x3i
              x0i = x1i - x3r
              a(j3) = wk3r * x0r - wk3i * x0i
              a(j3 + 1) = wk3r * x0i + wk3i * x0r
          end do
          wk1r = w(k2 + 2)
          wk1i = w(k2 + 3)
          wk3r = wk1r - 2 * wk2r * wk1i
          wk3i = 2 * wk2r * wk1r - wk1i
          do j = k + m, l + (k + m) - 2, 2
              j1 = j + l
              j2 = j1 + l
              j3 = j2 + l
              x0r = a(j) + a(j1)
              x0i = a(j + 1) + a(j1 + 1)
              x1r = a(j) - a(j1)
              x1i = a(j + 1) - a(j1 + 1)
              x2r = a(j2) + a(j3)
              x2i = a(j2 + 1) + a(j3 + 1)
              x3r = a(j2) - a(j3)
              x3i = a(j2 + 1) - a(j3 + 1)
              a(j) = x0r + x2r
              a(j + 1) = x0i + x2i
              x0r = x0r - x2r
              x0i = x0i - x2i
              a(j2) = -wk2i * x0r - wk2r * x0i
              a(j2 + 1) = -wk2i * x0i + wk2r * x0r
              x0r = x1r - x3i
              x0i = x1i + x3r
              a(j1) = wk1r * x0r - wk1i * x0i
              a(j1 + 1) = wk1r * x0i + wk1i * x0r
              x0r = x1r + x3i
              x0i = x1i - x3r
              a(j3) = wk3r * x0r - wk3i * x0i
              a(j3 + 1) = wk3r * x0i + wk3i * x0r
          end do
      end do
      end
c
      subroutine rftfsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      m = n / 2
      ks = 2 * nc / m
      kk = 0
      do j = 2, m - 2, 2
          k = n - j
          kk = kk + ks
          wkr = 0.5d0 - c(nc - kk)
          wki = c(kk)
          xr = a(j) - a(k)
          xi = a(j + 1) + a(k + 1)
          yr = wkr * xr - wki * xi
          yi = wkr * xi + wki * xr
          a(j) = a(j) - yr
          a(j + 1) = a(j + 1) - yi
          a(k) = a(k) + yr
          a(k + 1) = a(k + 1) - yi
      end do
      end
c
      subroutine rftbsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr, xi, yr, yi
      a(1) = -a(1)
      m = n / 2
      ks = 2 * nc / m
      kk = 0
      do j = 2, m - 2, 2
          k = n - j
          kk = kk + ks
          wkr = 0.5d0 - c(nc - kk)
          wki = c(kk)
          xr = a(j) - a(k)
          xi = a(j + 1) + a(k + 1)
          yr = wkr * xr + wki * xi
          yi = wkr * xi - wki * xr
          a(j) = a(j) - yr
          a(j + 1) = yi - a(j + 1)
          a(k) = a(k) + yr
          a(k + 1) = yi - a(k + 1)
      end do
      a(m + 1) = -a(m + 1)
      end
c
      subroutine dctsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      m = n / 2
      ks = nc / n
      kk = 0
      do j = 1, m - 1
          k = n - j
          kk = kk + ks
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          xr = wki * a(j) - wkr * a(k)
          a(j) = wkr * a(j) + wki * a(k)
          a(k) = xr
      end do
      a(m) = c(0) * a(m)
      end
c
      subroutine dstsub(n, a, nc, c)
      integer n, nc, j, k, kk, ks, m
      real*8 a(0 : n - 1), c(0 : nc - 1), wkr, wki, xr
      m = n / 2
      ks = nc / n
      kk = 0
      do j = 1, m - 1
          k = n - j
          kk = kk + ks
          wkr = c(kk) - c(nc - kk)
          wki = c(kk) + c(nc - kk)
          xr = wki * a(k) - wkr * a(j)
          a(k) = wkr * a(k) + wki * a(j)
          a(j) = xr
      end do
      a(m) = c(0) * a(m)
      end
c
