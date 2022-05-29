************************************************************************
*                                                                      *
* XS2C                                                                 *
*                                                                      *
* Purpose  Call of S2C                                                 *
*                                                                      *
* Subroutines/functions called  S2C                                    *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* PARX     FUNC                                                        *
* PARY     FUNC   Description of domain                                *
* TMAX     FUNC                                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* S2C                                                                  *
*                                                                      *
* Purpose  Check of triangulation                                      *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DCORVG   R*8    Cartesian coordinates of interior vertices           *
*                 Parameter of vretices on the boundary                *
* KVERT    I*4    Numbers of vertices in each element, counterclockwise*
* KADJ     I*4    Number of neighbouring element, 0 for boundary edges *
* KNPR     I*4    0    for interior vertices                           *
*                 IBCT for nodal points on boundary component IBCT     *
* BADJ     L*4    .FALSE.: KADJ  not available - skip checks for KADJ  *
* PARX     FUNC                                                        *
* PARY     FUNC   Description of domain                                *
* TMAX     FUNC                                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    <0  Invalid triangulation                            *
*                 -163  Two points not belonging to the same boundary  *
*                       component join a boundary edge                 *
*                 -154  KVERT contains 2 vertices with the same number *
*                 -159  The area of an element is 0                    *
*                 -158  The vertices of one triangle are numbered      *
*                       in the clockwise sense                         *
*                 -162  There are neighboured triangles without a      *
*                       joining edge                                   *
*                 -150  NVE invalid                                    *
*                 -151  NEL, NVT, or NBCT <=0                          *
*                 -152  KNPR(IVT) > NBCT                               *
*                 -153  Parameter of boundary vertex IVT too large     *
*                 -155  KVERT contains invalid zero entries            *
*                 -156  KVERT contains entries larger than NVT         *
*                 -157  KVERT(NNVE,.) not zero for NVE = 3             *
*                 -161  KADJ contains entries larger than NEL          *
*                 -160  KADJ(NNVE,.)  not zero for NVE = 3             *
*                                                                      *
************************************************************************
C
      SUBROUTINE XS2C(PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL PARX,PARY,TMAX
      SAVE /TRIAA/,/ERRCTL/
C
      IF (ICHECK.GE.997) CALL OTRC('XS2C  ','01/02/89')
C
      BADJ=LADJ.GT.0
      CALL S2C(DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LADJ)),
     *         KWORK(L(LNPR)),BADJ,PARX,PARY,TMAX)
C
      END
C
C
C
      SUBROUTINE S2C(DCORVG,KVERT,KADJ,KNPR,BADJ,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4,NNVE0=5)
      DIMENSION DCORVG(2,*),KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*)
      DIMENSION KN(NNVE0),KN1(NNVE0),X(NNVE),Y(NNVE)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='S2C'
      IF (ICHECK.GE.997) CALL OTRC('S2C   ','01/02/89')
C
      NVE0=3
      IF (KVERT(NNVE,1).NE.0) NVE0=NNVE
      IF (NVE.NE.NVE0) CALL WERR(-150,'S2C   ')
      IF (NEL.LE.0.OR.NVT.LE.0.OR.NBCT.LE.0) CALL WERR(-151,'S2C   ')
      IF (IER.NE.0) GOTO 99999
C
      DO 10 IVT=1,NVT
      IBCT=KNPR(IVT)
      IF (IBCT.NE.0) THEN
       IF (IBCT.GT.NBCT.OR.IBCT.LT.0) THEN
        WRITE (CPARAM,'(2I15)') IVT,IBCT
        CALL WERR(-152,'S2C   ')
       ELSE
        IF (DCORVG(1,IVT).GT.TMAX(IBCT).OR.DCORVG(1,IVT).LT.0D0) THEN
         WRITE (CPARAM,'(I15,D25.16)') IVT,DCORVG(1,IVT)
         CALL WERR(-153,'S2C   ')
        ENDIF
       ENDIF
      ENDIF
10    CONTINUE
      IF (IER.LT.0) GOTO 99999
C
      DO 100 IEL=1,NEL
      WRITE (CPARAM,'(I15)') IEL
C
      DO 110 IVE=1,NVE
110   KN(IVE)=KVERT(IVE,IEL)
      KN(NVE+1)=KN(1)
      DO 111 IVE=1,NVE
      IF (KN(IVE).EQ.KN(IVE+1)) CALL WERR(-154,'S2C   ')
      IF (KN(IVE).EQ.0)         CALL WERR(-155,'S2C   ')
      IF (KN(IVE).GT.NVT)       CALL WERR(-156,'S2C   ')
111   CONTINUE
      IF (NVE.EQ.3.AND.KVERT(4,IEL).NE.0) CALL WERR(-157,'S2C   ')
      IF (IER.NE.0) GOTO 99999
C
      DO 120 IVE=1,NVE
      X(IVE)=DCORVG(1,KN(IVE))
      IBCT=KNPR(KN(IVE))
      IF (IBCT.EQ.0) THEN
       Y(IVE)=DCORVG(2,KN(IVE))
      ELSE
       Y(IVE)=PARY(X(IVE),IBCT)
       X(IVE)=PARX(X(IVE),IBCT)
      ENDIF
120   CONTINUE
C
      IF (NVE.EQ.3) THEN
       AREA=(X(2)-X(1))*(Y(3)-Y(1))-(X(3)-X(1))*(Y(2)-Y(1))
       IF (AREA.LT.0D0) CALL WERR(-158,'S2C   ')
       IF (AREA.GE.0D0.AND.AREA.LT.1D-70) CALL WERR(-159,'2SC   ')
      ELSE
       DO 121 IVE=1,NVE-1
       IVE1=IVE+1
       IVE2=MOD(IVE1,NVE)+1
       AREA =(X(IVE1)-X(IVE))*(Y(IVE2)-Y(IVE))
     *      -(X(IVE2)-X(IVE))*(Y(IVE1)-Y(IVE))
       IF (AREA.LT.0D0) THEN
         CALL WERR(-158,'S2C   ')
         WRITE (*,*) 'IVT = ',KN(1),KN(2),KN(3),KN(4)
         WRITE (*,*) 'X   = ',X(1),X(2),X(3),X(4)
         WRITE (*,*) 'Y   = ',Y(1),Y(2),Y(3),Y(4)
       ENDIF
121    IF (AREA.GE.0D0.AND.AREA.LT.1D-70) CALL WERR(-159,'S2C   ')
      ENDIF
      IF (IER.NE.0) GOTO 99999
C
      IF (.NOT.BADJ) GOTO 100
      IF (NVE.EQ.3.AND.KADJ(NNVE,IEL).NE.0) THEN
       CALL WERR(-160,'S2C   ')
       GOTO 99999
      ENDIF
C
      DO 130 IME=1,NVE
      JME=KADJ(IME,IEL)
      IF (JME.GT.NEL) THEN
       CALL WERR(-161,'S2C   ')
       GOTO 99999
      ENDIF
      IF (JME.NE.0) THEN
       K1=KN(IME)
       K=KN(IME+1)
       DO 131 JJVE=1,NVE
131    KN1(JJVE)=KVERT(JJVE,JME)
       KN1(NVE+1)=KN1(1)
       DO 132 JJVE=1,NVE
132    IF (KN1(JJVE).EQ.K) GOTO 133
135    CALL WERR(-162,'S2C   ')
       GOTO 99999
133    IF (KN1(JJVE+1).EQ.K1) GOTO 130
       GOTO 135
      ENDIF
130   CONTINUE
C
      DO 140 IME=1,NVE
      IF (KADJ(IME,IEL).EQ.0) THEN
       IF (KNPR(KN(IME)).NE.KNPR(KN(IME+1)).OR.KNPR(KN(IME)).EQ.0) THEN
        CALL WERR(-163,'S2C   ')
        GOTO 99999
       ENDIF
      ENDIF
140   CONTINUE
C
100   CONTINUE
C
99999 END
