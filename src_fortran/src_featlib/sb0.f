************************************************************************
*                                                                      *
* XSB0                                                                 *
*                                                                      *
* Purpose  Adjust dimensions of DCORVG, KVERT, KADJ and KNPR           *
*          Call of SB0                                                 *
*                                                                      *
* Subroutines/functions called  SB0, ZLEN, ZDISP, ZNEW, ZCPY           *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IDISP    I*4    =1 Release free space on the arrays                  *
*                    after determination of the new subdivision        *
* For the description of the remaining parameters see SB0              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* SB0                                                                  *
*                                                                      *
* Purpose  Regular refinement of a given quadrilateral subdivision     *
*          of a two-dimensional domain                                 *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  06/26/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DCORVG   R*8    Cartesian coordinates of interior vertices           *
*                 Parameter of vertices on the boundary                *
* KVERT    I*4    Numbers of vertices in each element, counterclockwise*
* KADJ     I*4    Number of neighbouring element, 0 for boundary edges *
*                                                                      *
*                 Convention for numbers of vertices and neighbours    *
*                                                                      *
*                             P4     N3      P3                        *
*                                * - - - - *                           *
*                                !         !                           *
*                             N4 !         ! N2                        *
*                                !         !                           *
*                                * - - - - *                           *
*                             P1     N1      P2                        *
*                                                                      *
* KNPR     I*4    0    for interior vertices                           *
*                 IBCT for nodal points on boundary component IBCT     *
* KMM      I*4    KMM(1,IBCT) and KMM(2,IBCT) contain the numbers      *
*                 of vertices with minimum and maximum parameter       *
*                                                                      *
* NNEL     I*4    Maximum dimension provided for KVERT                 *
* NNVT     I*4    Maximum dimension provided for DCORVG, KNPR          *
*                                                                      *
* NFINE    I*4    Number of regular subdivisions                       *
* S2DI     SUBR   S2DI(X1,Y1,X2,Y1,XT,YT) gives the coordinates        *
*                 of the new vertex on the line (X1,Y1)-(X2,Y2)        *
*                 S2DI=S2DI0 takes the midpoint                        *
* S2DB     SUBR   S2DB(X1,X2,XT,BMM,IBCT,PARX,PARY,TMAX(IBCT))         *
*                 gives the parameter of the new boundary vertex       *
*                 on the boundary segment joining the vertices with    *
*                 parameter X1 and X2                                  *
*                 BMM=.TRUE. means that X1,X2 are the minimum and      *
*                 maximum parameter on boundary component IBCT         *
* PARX     FUNC                                                        *
* PARY     FUNC   Double precision functions describing the domain     *
* TMAX     FUNC                                                        *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8    As above                                             *
* KVERT    I*4    As above                                             *
* KADJ     I*4    As above                                             *
* KNPR     I*4    As above                                             *
* KMM      I*4    As above                                             *
*                                                                      *
* NNEL     I*4    Dimension needed for KVERT, KADJ                     *
* NNVT     I*4    Dimension needed for DCORVG, KNPR                    *
*                                                                      *
* NEL      I*4    Number of elements of final triangulation            *
* NVT      I*4    Number if vertices of final triangulation            *
*                                                                      *
* IER      I*4     0  No error                                         *
*                  1  NNEL or NNVT too small                           *
*                     Requirements                                     *
*                     NNEL >=NEL(NFINE)                                *
*                     NNVT >=NVT(NFINE)                                *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSB0(NFINE,IDISP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNVE=4,NNARR=999)
      INCLUDE 'nnwork.inc'
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XSB0'
      IF (ICHECK.GE.997) CALL OTRC('XSB0  ','01/02/89')
C
      CALL ZLEN(LVERT,ILEN)
      NNEL=ILEN/NNVE
      CALL ZLEN(LADJ,ILEN)
      NNEL1=ILEN/NNVE
      NNEL=MIN(NNEL,NNEL1)
      CALL ZLEN(LCORVG,ILEN)
      NNVT=ILEN/2
      NNEL1=NNEL
      NNVT1=NNVT
C
      CALL SB0(DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LNPR)),NFINE,NNEL,NNVT)
C
      IF (IER) 99999,1,2
C
1     IF (IDISP.EQ.1) THEN
       CALL ZDISP(NNVE*NNEL,LVERT,'KVERT ')
       CALL ZDISP(NNVE*NNEL,LADJ,'KADJ  ')
       CALL ZDISP(2*NNVT,LCORVG,'DCORVG')
       CALL ZDISP(NNVT,LNPR,'KNPR  ')
      ENDIF
      GOTO 99999
C
2     IF (NNEL.GT.NNEL1) THEN
       CALL ZNEW(NNVE*NNEL,3,LV1,'KVERT ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LVERT,'KVOLD ',LV1,'KVERT ')
       CALL ZDISP(0,LVERT,'KVOLD ')
       LVERT=LV1
C
       CALL ZNEW(NNVE*NNEL,3,LA1,'KADJ  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LADJ,'KAOLD ',LA1,'KADJ  ')
       CALL ZDISP(0,LADJ,'KAOLD ')
       LADJ=LA1
      ENDIF
C
      IF (NNVT.GT.NNVT1) THEN
       CALL ZNEW(2*NNVT,1,LCV1,'DCORVG')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LCORVG,'DCVOLD',LCV1,'DCORVG')
       CALL ZDISP(0,LCORVG,'DCVOLD')
       LCORVG=LCV1
      ENDIF
C
      IF (NNVT.GT.NNVT1) THEN
       CALL ZNEW(NNVT,3,LNPR1,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LNPR,'KNPOLD',LNPR1,'KNPR  ')

       CALL ZDISP(0,LNPR,'KNPOLD')
       LNPR=LNPR1
      ENDIF
C
      CALL SB0(DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *         KWORK(L(LADJ)),KWORK(L(LNPR)),NFINE,NNEL,NNVT)
C
99999 END
C
C
C
      SUBROUTINE SB0(DCORVG,KVERT,KADJ,KNPR,NFINE,NNEL,NNVT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*),KEL(NNVE),KEL1(NNVE),X(NNVE),Y(NNVE)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KNPR(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
      DATA IVT/0/
C
      NNEL1=NNEL
      NNVT1=NNVT
      NNEL=NEL
      NNVT=NVT
      NNMT=NNEL+NNVT+NBCT-2
C
      IF (NFINE.GT.0) THEN
       DO 1 IFINE=1,NFINE
       NNVT=NNVT+NNMT+NNEL
       NNEL=4*NNEL
1      NNMT=NNVT+NNEL+NBCT-2
      ENDIF
C
      IF (NNEL.GT.NNEL1.OR.NNVT.GT.NNVT1) THEN
       IER=1
       GOTO 99999
      ENDIF
C
      IER=0
C
      DO 100 IFINE=1,NFINE
      NVT1=NVT
      NEL1=NEL
      DO 110 IEL=1,NEL1
C
C *** Calculate coordinates of parameter of new vertices
C
      DO 120 IVE1=1,NVE
C
      IF (KADJ(IVE1,IEL).LT.IEL .AND.
     *    KADJ(IVE1,IEL).NE.0) GOTO 120   ! Already done
C
      NVT1=NVT1+1
C
      IVT1=KVERT(IVE1,IEL)
C
      IVE2=MOD(IVE1,NVE)+1
      IVT2=KVERT(IVE2,IEL)
C
      DCORVG(1,NVT1)=0.5D0*(DCORVG(1,IVT1)+DCORVG(1,IVT2))
      DCORVG(2,NVT1)=0.5D0*(DCORVG(2,IVT1)+DCORVG(2,IVT2))
C
      IF (KADJ(IVE1,IEL).NE.0) THEN
       KNPR(NVT1)=0
      ELSE
       KNPR(NVT1)=1
      ENDIF
120   CONTINUE
C
C *** Actualization of KVERT
C
      KEL(1)=IEL
      DO 130 I=1,3
130   KEL(I+1)=NEL+I
C
      DO 131 I=1,4
131   KVERT(1,KEL(I))=KVERT(I,IEL)
C
      DO 140 IVE1=1,NVE
      IVE2=MOD(IVE1,NVE)+1
      IEL1=KADJ(IVE1,IEL)
C
      IF (IEL1.GE.IEL.OR.IEL1.EQ.0) THEN
       NVT=NVT+1
       KVERT(2,KEL(IVE1))=NVT
       KVERT(NVE,KEL(IVE2))=NVT
C
      ELSE
C
C *** Find out the number of the vertex (already exists)
C
       KEL1(1)=IEL1
       DO 141 JEL=1,3
141    KEL1(JEL+1)=KADJ(2,KEL1(JEL))
C
       DO 142 I=1,4
       IF (KADJ(NVE,KEL1(I)).EQ.IEL) THEN
        IVT=KVERT(NVE,KEL1(I))
        GOTO 145
       ENDIF
142    CONTINUE
C
145    KVERT(2,KEL(IVE1))=IVT
       KVERT(NVE,KEL(IVE2))=IVT
C
      ENDIF
C
140   CONTINUE
C
C *** Actualization of KADJ
C
      DO 150 IVE1=1,NVE
      IVE2=MOD(IVE1,NVE)+1
      JEL=KADJ(IVE1,IEL)
C
      IF (JEL.EQ.0.OR.JEL.GE.IEL) THEN
       KADJ(1,KEL(IVE1))=JEL
       KADJ(NVE,KEL(IVE2))=JEL
C
      ELSE
C
       KEL1(1)=JEL
       DO 151 I=1,3
151    KEL1(I+1)=KADJ(2,KEL1(I))
       DO 155 I=1,4
       IF (KADJ(NVE,KEL1(I)).EQ.IEL) THEN
        KADJ(1,KEL(IVE1))=KEL1(I)
        KADJ(NVE,KEL1(I))=KEL(IVE1)
        I1=I-1
        IF (I1.EQ.0) I1=4
        KADJ(1,KEL1(I1))=KEL(IVE2)
        KADJ(NVE,KEL(IVE2))=KEL1(I1)
        GOTO 150
       ENDIF
155    CONTINUE
      ENDIF
C
150   CONTINUE
C
      DO 159 I=1,4
      I1=MOD(I,NVE)+1
      KADJ(2,KEL(I))=KEL(I1)
159   KADJ(3,KEL(I1))=KEL(I)
      NEL=NEL+3
C
110   CONTINUE
C
C *** Calculation of coordinates and numbers of midpoints
C
      DO 160 IEL=1,NEL1
C
      KEL(1)=IEL
      DO 161 I=1,3
161   KEL(I+1)=KADJ(2,KEL(I))
C
      DO 162 IVE1=1,NVE
      IVT=KVERT(2,KEL(IVE1))
      IBCT=KNPR(IVT)
      X(IVE1)=DCORVG(1,IVT)
      Y(IVE1)=DCORVG(2,IVT)
162   CONTINUE
C
      NVT=NVT+1
      X21=X(2)-X(1)
      X31=X(3)-X(1)
      X24=X(2)-X(4)
      Y21=Y(2)-Y(1)
      Y31=Y(3)-Y(1)
      Y24=Y(2)-Y(4)
      ALPHA=(X21*Y24-Y21*X24)/(X31*Y24-Y31*X24)
      DCORVG(1,NVT)=X(1)+ALPHA*X31
      DCORVG(2,NVT)=Y(1)+ALPHA*Y31
      KNPR(NVT)=0
C
      DO 163 I=1,4
163   KVERT(3,KEL(I))=NVT
160   CONTINUE
C
100   CONTINUE
C
99999 END
