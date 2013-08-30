      SUBROUTINE YTOA(Y,A)
      IMPLICIT NONE
      INCLUDE "constants.h"
      INCLUDE "n_power.for"
      INTEGER NS,NSC,I1,I2,I3,K,INDEX
      DOUBLE PRECISION Y,A,PREF

      PARAMETER(NS=2**NN/6)

      DIMENSION A(0:NS,0:NS,0:NS),Y(NT),K(3),PREF(2),INDEX(2)

      DO 10 I1=0,NS
         DO 20 I2=0,NS
            DO 30 I3=0,NS
               A(I1,I2,I3)=0D0
               IF(I2.GT.I3) THEN
                  K(1)=2*I1
                  K(2)=2*I2
                  K(3)=2*I3
                  CALL MAPTOFD(K,INDEX,PREF)
                  IF(INDEX(1).NE.0) A(I1,I2,I3)=-1D0*PREF(1)*Y(INDEX(1))
                  IF(INDEX(2).NE.0) A(I1,I2,I3)=A(I1,I2,I3)-
     C                 PREF(2)*Y(INDEX(2))
               ELSE
                  IF((I1.EQ.I2).AND.(I1.EQ.I3)) GOTO 30
                  IF((I1.EQ.NS).OR.(I3.EQ.NS)) GOTO 30
                  K(1)=2*I1+1
                  K(2)=2*I2+1
                  K(3)=2*I3+1
                  CALL MAPTOFD(K,INDEX,PREF)
                  IF(INDEX(1).NE.0) A(I1,I2,I3)=PREF(1)*Y(INDEX(1))
                  IF(INDEX(2).NE.0) A(I1,I2,I3)=A(I1,I2,I3)+
     C                 PREF(2)*Y(INDEX(2))
               END IF
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE

      RETURN
      END

C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE ATOY(A,Y)
      IMPLICIT NONE
      INCLUDE "constants.h"
      INCLUDE "n_power.for"
      INTEGER NS,NSC,K,I1,I2,I3,I
      DOUBLE PRECISION Y,A

      PARAMETER(NS=2**NN/6)

      DIMENSION A(0:NS,0:NS,0:NS),Y(NT),K(3)

c      COMMON /VOR1/A

      DO 10 I=1,NT
         CALL NK(I,K)
         IF(MOD(K(1),2).EQ.0) THEN
            I1=K(1)/2
            I2=K(2)/2
            I3=K(3)/2
            Y(I)=A(I1,I3,I2)
         ELSE
            I1=(K(1)-1)/2
            I2=(K(2)-1)/2
            I3=(K(3)-1)/2
            Y(I)=A(I1,I2,I3)
         END IF
 10   CONTINUE

      RETURN
      END

C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE MAPTOFD(K,INDEX,PREF)
      IMPLICIT NONE
      INCLUDE "constants.h"
      INTEGER K,SGN,INDEX,P,I,D,Q
      DOUBLE PRECISION PREF

      DIMENSION K(3),P(3),INDEX(2),PREF(2),Q(3)

      IF((K(2)*K(3)).EQ.0) THEN
         INDEX(1)=0
         INDEX(2)=0
         RETURN
      END IF
      SGN=SIGN(1,K(2)*K(3))
      DO 10 I=1,3
         P(I)=ABS(K(I))
 10   CONTINUE
      IF((P(2).EQ.P(3)).AND.(MOD(P(1),2).EQ.0)) THEN
         INDEX(1)=0
         INDEX(2)=0
         RETURN
      END IF
      IF(P(2).GT.P(3)) THEN
         D=P(2)
         P(2)=P(3)
         P(3)=D
         IF(MOD(P(1),2).EQ.0) SGN=-1*SGN
      END IF

      IF(P(3).GE.P(1)) THEN
         IF(P(3).GT.P(2)) THEN
            CALL KN(P,INDEX(1))
            INDEX(2)=0
            PREF(1)=DFLOAT(SGN)
         ELSE
            CALL RROT(P,Q)
            CALL KN(Q,INDEX(1))
            INDEX(2)=0
            PREF(1)=DFLOAT(-2*SGN*P(2))/DFLOAT(P(1))
         END IF
      ELSE
         CALL LROT(P,Q)
         CALL KN(Q,INDEX(1))
         PREF(1)=DFLOAT(-SGN*P(2))/DFLOAT(P(1))
         D=P(1)
         P(1)=P(3)
         P(3)=D
         IF(MOD(P(1),2).EQ.0) SGN=-1*SGN
         CALL KN(P,INDEX(2))
         PREF(2)=DFLOAT(-SGN*P(1))/DFLOAT(P(3))
      END IF

      RETURN
      END

C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE KN(K,INDEX)
      IMPLICIT NONE
      INCLUDE "constants.h"
      INTEGER K,INDEX
      DIMENSION K(3)

      IF(MOD(K(1),2).EQ.0) THEN
         INDEX=(2*(K(3)/2-1)**3+3*(K(3)/2-1)**2-5*(K(3)/2-1))/6
     C        +(K(2)/2-1)*(K(3)/2+1)+(K(1)/2+1)
            ELSE
         INDEX=(((K(3)-3)/2)**3+3*((K(3)-3)/2)**2+2*((K(3)-3)/2))/3
     C              +((K(2)-1)/2)*((K(3)+1)/2)+(K(1)+1)/2 +NE
      END IF

      RETURN
      END

C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE NK(INDEX,K)
      IMPLICIT NONE
      INCLUDE "constants.h"
      INTEGER INDEX,K,K3L,FLAG,I,II,III,D
      DIMENSION K(3),K3L(NT)

      COMMON /AUX/ K3L

      DATA FLAG/0/
      IF(FLAG.NE.0) GO TO 1
      FLAG=1
      III=0
      DO 10 I=2,N/2
         DO 20 II=1,(I+1)*(I-1)
            III=III+1
            K3L(III)=I
 20      CONTINUE
 10   CONTINUE
      DO 30 I=1,(N-1)/2
         DO 40 II=1,(I+1)*I
            III=III+1
            K3L(III)=I
 40      CONTINUE
 30   CONTINUE
      
 1    IF(INDEX.LE.NE) THEN

         K(3)=2*K3L(INDEX)
         D=INDEX-(2*(K3L(INDEX)-1)**3+3*(K3L(INDEX)-1)**2
     C        -5*(K3L(INDEX)-1))/6
         K(2)=2*((D-1)/(K3L(INDEX)+1) +1)
         K(1)=2*MOD(D-1,K3L(INDEX)+1)

         ELSE

         K(3)=2*K3L(INDEX)+1
         D=INDEX-NE-((K3L(INDEX)-1)**3+3*(K3L(INDEX)-1)**2
     C        +2*(K3L(INDEX)-1))/3
         K(2)=2*((D-1)/(K3L(INDEX)+1)) +1
         K(1)=2*MOD(D-1,K3L(INDEX)+1)+1

      END IF

      RETURN
      END

C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE RROT(A,B)
      IMPLICIT NONE
      INTEGER A,B
      DIMENSION A(3),B(3)
      B(1)=A(3)
      B(2)=A(1)
      B(3)=A(2)
      RETURN
      END

C23456789012345678901234567890123456789012345678901234567890123456789012

      SUBROUTINE LROT(A,B)
      IMPLICIT NONE
      INTEGER A,B
      DIMENSION A(3),B(3)
      B(1)=A(2)
      B(2)=A(3)
      B(3)=A(1)
      RETURN
      END
