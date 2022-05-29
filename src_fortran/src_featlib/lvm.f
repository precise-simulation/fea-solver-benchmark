************************************************************************
*                                                                      *
* LVMn                                                                 *
*                                                                      *
* Purpose  Vector multiply and add                                     *
*          n=1 double precision version                                *
*          n=2 single precision version                                *
*          n=3 mixed precision version                                 *
*                                                                      *
* Subroutines/functions called   None                                  *
*                                                                      *
* Version from  02/20/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8                                                         *
* DX1      R*8    Input vectors                                        *
* DX2      R*8                                                         *
* VX       R*4                                                         *
* VX1      R*4    Input vectors (single precision)                     *
* VX2      R*4                                                         *
* NX       I*4    Number of elements                                   *
* A1       R*8    Coefficients of linear combinations, see below       *
* A2       R*8                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    DX(IX) := A1*DX1(IX)*DX2(IX)+A2*DX(IX)        (n=1)  *
* VX       R*8    VX(IX) := A1*VX1(IX)*VX2(IX)+A2*VX(IX)        (n=2)  *
* DX       R*8    DX(IX) := A1*DX1(IX)*DBLE(VX2(IX))+A2*DX(IX)  (n=3)  *
*                                                                      *
************************************************************************
C
      SUBROUTINE LVM1(DX1,DX2,DX,NX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),DX1(*),DX2(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LVM1  ','02/20/91')
C
      IF (A2.EQ.0D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 10 IX=1,NX
10      DX(IX)=DX1(IX)*DX2(IX)
       ELSE IF (A1.EQ.-1D0) THEN
        DO 11 IX=1,NX
11      DX(IX)=-DX1(IX)*DX2(IX)
       ELSE
        DO 12 IX=1,NX
12      DX(IX)=A1*DX1(IX)*DX2(IX)
       ENDIF
C
      ELSE IF (A2.EQ.1D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 20 IX=1,NX
20      DX(IX)=DX(IX)+DX1(IX)*DX2(IX)
       ELSE IF (A1.EQ.-1D0) THEN
        DO 21 IX=1,NX
21      DX(IX)=DX(IX)-DX1(IX)*DX2(IX)
       ELSE
        DO 22 IX=1,NX
22      DX(IX)=DX(IX)+A1*DX1(IX)*DX2(IX)
       ENDIF
C
      ELSE IF (A2.EQ.-1D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 30 IX=1,NX
30      DX(IX)=-DX(IX)+DX1(IX)*DX2(IX)
       ELSE IF (A1.EQ.-1D0) THEN
        DO 31 IX=1,NX
31      DX(IX)=-DX(IX)-DX1(IX)*DX2(IX)
       ELSE
        DO 32 IX=1,NX
32      DX(IX)=-DX(IX)+A1*DX1(IX)*DX2(IX)
       ENDIF
C
      ELSE
C
       DO 40 IX=1,NX
40     DX(IX)=A2*DX(IX)+A1*DX1(IX)*DX2(IX)
C
      ENDIF
C
      END
C
C
C
      SUBROUTINE LVM2(VX1,VX2,VX,NX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      REAL A10,A20
      DIMENSION VX(*),VX1(*),VX2(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LVM2  ','02/20/91')
C
      IF (A2.EQ.0D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 10 IX=1,NX
10      VX(IX)=VX1(IX)*VX2(IX)
       ELSE IF (A1.EQ.-1D0) THEN
        DO 11 IX=1,NX
11      VX(IX)=-VX1(IX)*VX2(IX)
       ELSE
        A10=REAL(A1)
        DO 12 IX=1,NX
12      VX(IX)=A10*VX1(IX)*VX2(IX)
       ENDIF
C
      ELSE IF (A2.EQ.1D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 20 IX=1,NX
20      VX(IX)=VX(IX)+VX1(IX)*VX2(IX)
       ELSE IF (A1.EQ.-1D0) THEN
        DO 21 IX=1,NX
21      VX(IX)=VX(IX)-VX1(IX)*VX2(IX)
       ELSE
        A10=REAL(A1)
        DO 22 IX=1,NX
22      VX(IX)=VX(IX)+A10*VX1(IX)*VX2(IX)
       ENDIF
C
      ELSE IF (A2.EQ.-1D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 30 IX=1,NX
30      VX(IX)=-VX(IX)+VX1(IX)*VX2(IX)
       ELSE IF (A1.EQ.-1D0) THEN
        DO 31 IX=1,NX
31      VX(IX)=-VX(IX)-VX1(IX)*VX2(IX)
       ELSE
        A10=REAL(A1)
        DO 32 IX=1,NX
32      VX(IX)=-VX(IX)+A10*VX1(IX)*VX2(IX)
       ENDIF
C
      ELSE
C
       A10=REAL(A1)
       A20=REAL(A2)
       DO 40 IX=1,NX
40     VX(IX)=A20*VX(IX)+A10*VX1(IX)*VX2(IX)
C
      ENDIF
C
      END
C
C
C
      SUBROUTINE LVM3(DX1,VX2,DX,NX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),DX1(*),VX2(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LVM3  ','02/20/91')
C
      IF (A2.EQ.0D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 10 IX=1,NX
10      DX(IX)=DX1(IX)*DBLE(VX2(IX))
       ELSE IF (A1.EQ.-1D0) THEN
        DO 11 IX=1,NX
11      DX(IX)=-DX1(IX)*DBLE(VX2(IX))
       ELSE
        DO 12 IX=1,NX
12      DX(IX)=A1*DX1(IX)*DBLE(VX2(IX))
       ENDIF
C
      ELSE IF (A2.EQ.1D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 20 IX=1,NX
20      DX(IX)=DX(IX)+DX1(IX)*DBLE(VX2(IX))
       ELSE IF (A1.EQ.-1D0) THEN
        DO 21 IX=1,NX
21      DX(IX)=DX(IX)-DX1(IX)*DBLE(VX2(IX))
       ELSE
        DO 22 IX=1,NX
22      DX(IX)=DX(IX)+A1*DX1(IX)*DBLE(VX2(IX))
       ENDIF
C
      ELSE IF (A2.EQ.-1D0) THEN
C
       IF (A1.EQ.1D0) THEN
        DO 30 IX=1,NX
30      DX(IX)=-DX(IX)+DX1(IX)*DBLE(VX2(IX))
       ELSE IF (A1.EQ.-1D0) THEN
        DO 31 IX=1,NX
31      DX(IX)=-DX(IX)-DX1(IX)*DBLE(VX2(IX))
       ELSE
        DO 32 IX=1,NX
32      DX(IX)=-DX(IX)+A1*DX1(IX)*DBLE(VX2(IX))
       ENDIF
C
      ELSE
C
       DO 40 IX=1,NX
40     DX(IX)=A2*DX(IX)+A1*DX1(IX)*DBLE(VX2(IX))
C
      ENDIF
C
      END
