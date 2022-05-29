************************************************************************
*                                                                      *
* XS2M                                                                 *
*                                                                      *
* Purpose  Adjust dimensions of DCORMG, KMID, KADJ, KMEL, and KNPR     *
*          Call of S2M, S2MEL                                          *
*                                                                      *
* Subroutines/functions called  S2M, S2MEL, ZLEN, ZNEW, ZDISP, ZCPY    *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IDISP    I*4    =1 Release free space on the arrays                  *
* IADJ     I*4    =0 Array KADJ is used to store KMID                  *
*                    Old information is destroyed                      *
*                                                                      *
* For the description of the remaining parameters see S2M              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* S2M                                                                  *
*                                                                      *
* Purpose  Determination of numbers and coordinates of midpoints       *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DCORVG   R*8    Cartesian coordinates of interior vertices           *
*                 Parameter of vertices on the boundary                *
* KVERT    I*4    Numbers of vertices in each element, counterclockwise*
* KADJ     I*4    Number of neighbouring element, 0 for boundary edges *
* KNPR     I*4    0    for interior vertices                           *
*                 IBCT for nodal points on boundary component IBCT     *
* KMM      I*4    KMM(1,IBCT) and KMM(2,IBCT) contain the numbers      *
*                 of vertices with minimum and maximum parameter       *
* IMID     I*4    =1 Determine KMID, KNPR for midpoints and            *
*                    DCORVG(2,.) for midpoints on the boundary         *
*                 =2 Calculate KNPR and DCORMG for midpoints           *
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
* DCORMG   R*8    Coordinates of midpoints of edges if IMID=2          *
* KMID     I*4    Numbers of midpoints (or edges, resp.)               *
*                 Range of numbers of edges is NVT+1,...,NVT+NMT       *
* KNPR     I*4    As above for 1<=INPR<=NVT                            *
*                 For NVT+1,...,NVT+NMT it contains                    *
*                 0   for interior midpoints (or edges)                *
*                 IVT for boundary midpoints, where IVT is the         *
*                     preceeding boundary vertex                       *
*                                                                      *
************************************************************************
C
      SUBROUTINE XS2M(IMID,IADJ,IDISP,S2DI,S2DB,PARX,PARY,TMAX)
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
      EXTERNAL S2DI,S2DB,PARX,PARY,TMAX
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XS2M'
      IF (ICHECK.GE.997) CALL OTRC('XS2M  ','04/12/91')
C
      IF (IMID.EQ.0) THEN
       NMT=0
       GOTO 99999
      ENDIF
C
      NMT=NVE*NEL
        DO 1 IEL=1,NEL
        IPTR=(IEL-1)*NNVE
        DO 1 IVE=1,NVE
        I=IPTR+IVE-1
        IF (KWORK(L(LADJ)+I).EQ.0) NMT=NMT+1
1     CONTINUE
      NMT=NMT/2
C
CC      NMT=NVT+NEL+NBCT-2
      NVMT=NVT+NMT
C
      CALL ZLEN(LNPR,ILEN)
      IF (NVMT.GT.ILEN) THEN
       CALL ZNEW(NVMT,3,LNPR1,'KNPR  ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LNPR,'KNPOLD',LNPR1,'KNPR  ')
       CALL ZDISP(0,LNPR,'KNPOLD')
       LNPR=LNPR1
      ENDIF
C
      CALL ZLEN(LCORMG,ILEN)
      IF ((IMID.GE.2).AND.(2*NMT.GT.ILEN)) THEN
       CALL ZNEW(2*NMT,1,LCM1,'DCORMG')
       IF (IER.NE.0) GOTO 99999
       IF (ILEN.GT.0) CALL ZCPY(LCORMG,'DCMOLD',LCM1,'DCORMG')
       IF (LCORMG.NE.0) CALL ZDISP(0,LCORMG,'DCMOLD')
       LCORMG=LCM1
      ENDIF
C
      CALL ZLEN(LMID,ILEN)
      IF ((IADJ.GT.0).AND.(NNVE*NEL.GT.ILEN)) THEN
       CALL ZNEW(NNVE*NEL,3,LM1,'KMID  ')
       IF (IER.NE.0) GOTO 99999
       IF (ILEN.GT.0) CALL ZDISP(0,LMID,'KMOLD ')
       LMID=LM1
      ENDIF
C
      IF (IADJ.EQ.0) LMID=LADJ
C
      CALL ZLEN(LMEL,ILEN)
      IF ((2*NMT.GT.ILEN)) THEN
       CALL ZNEW(2*NMT,3,LM1,'KMEL  ')
       IF (IER.NE.0) GOTO 99999
       IF (ILEN.GT.0) CALL ZDISP(0,LMEL,'KMEOLD')
       LMEL=LM1
      ENDIF

C
      CALL S2M(DWORK(L(LCORVG)),DWORK(L(LCORMG)),KWORK(L(LVERT)),
     *         KWORK(L(LMID)),KWORK(L(LADJ)),KWORK(L(LNPR)),
     *         KWORK(L(LMM)),IMID,S2DI,S2DB,PARX,PARY,TMAX)
C
      IF (IADJ.EQ.0) LADJ=0
C
      CALL S2MEL(KWORK(L(LMEL)),KWORK(L(LMID)),IMID,NEL,NMT,NVT)
C
      IF (IDISP.EQ.1) THEN
       CALL ZDISP(NNVE*NEL,LMID,'KMID  ')
       CALL ZDISP(2*NMT,LMEL,'KMEL  ')
       IF (IMID.LE.1.AND.LCORMG.NE.0) CALL ZDISP(2*NMT,LCORMG,'DCORMG')
       CALL ZDISP(NVMT,LNPR,'KNPR  ')
       IF (LADJ.GT.0) CALL ZDISP(NNVE*NEL,LADJ,'KADJ  ')
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE S2M(DCORVG,DCORMG,KVERT,KMID,KADJ,KNPR,KMM,
     *               IMID,S2DI,S2DB,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      DIMENSION DCORVG(2,*),DCORMG(2,*)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*),KADJ(NNVE,*),KNPR(*),KMM(2,*)
      EXTERNAL S2DI,S2DB,PARX,PARY,TMAX
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /TRIAD/,/ERRCTL/
        DATA JVE2/0/
C
      IF (ICHECK.GE.997) CALL OTRC('S2M   ','01/02/89')
C
      NVMT=NVT
C
      DO 100 IEL=1,NEL
      DO 110 IVE1=1,NVE
      JEL=KADJ(IVE1,IEL)
C
      IF (JEL.EQ.0) THEN
C
       NVMT=NVMT+1
       KMID(IVE1,IEL)=NVMT
       IVT1=KVERT(IVE1,IEL)
       IVE2=MOD(IVE1,NVE)+1
       IVT2=KVERT(IVE2,IEL)
       IF (DCORVG(1,IVT1).GT.DCORVG(1,IVT2)) THEN
        JAUX=IVT2
        IVT2=IVT1
        IVT1=JAUX
       ENDIF
C
       IBCT=KNPR(IVT1)
       TTMAX=TMAX(IBCT)
       BMM=KMM(1,IBCT).EQ.IVT1.AND.KMM(2,IBCT).EQ.IVT2
       IF (BMM) THEN
        IVTAUX=IVT2
       ELSE
        IVTAUX=IVT1
       ENDIF
       KNPR(NVMT)=IVTAUX
       CALL S2DB(DCORVG(1,IVT1),DCORVG(1,IVT2),DCORVG(2,IVTAUX),
     *           BMM,IBCT,PARX,PARY,TTMAX)
C
       IF (IMID.GE.2) THEN
        DCORMG(1,NVMT-NVT)=PARX(DCORVG(2,IVTAUX),IBCT)
        DCORMG(2,NVMT-NVT)=PARY(DCORVG(2,IVTAUX),IBCT)
       ENDIF
C
      ELSE IF (JEL.GE.IEL) THEN
C
       NVMT=NVMT+1
       KMID(IVE1,IEL)=NVMT
       KNPR(NVMT)=0
C
       IF (IMID.GE.2) THEN
        IVT1=KVERT(IVE1,IEL)
        IVE2=MOD(IVE1,NVE)+1
        IVT2=KVERT(IVE2,IEL)
        IBCT1=KNPR(IVT1)
        IBCT2=KNPR(IVT2)
C
        IF (IBCT1.EQ.0) THEN
         X1=DCORVG(1,IVT1)
         Y1=DCORVG(2,IVT1)
        ELSE
         X1=PARX(DCORVG(1,IVT1),IBCT1)
         Y1=PARY(DCORVG(1,IVT1),IBCT1)
        ENDIF
        IF (IBCT2.EQ.0) THEN
         X2=DCORVG(1,IVT2)
         Y2=DCORVG(2,IVT2)
        ELSE
         X2=PARX(DCORVG(1,IVT2),IBCT2)
         Y2=PARY(DCORVG(1,IVT2),IBCT2)
        ENDIF
C
        CALL S2DI(X1,Y1,X2,Y2,DCORMG(1,NVMT-NVT),DCORMG(2,NVMT-NVT))
       ENDIF
C
      ELSE
C
       DO 120 JVE1=1,NVE
       JVE2=JVE1-1
       IF (JVE2.EQ.0) JVE2=NVE
       IVE2=MOD(IVE1,NVE)+1
       IF (KVERT(JVE1,JEL).EQ.KVERT(IVE1,IEL).AND.
     *     KVERT(JVE2,JEL).EQ.KVERT(IVE2,IEL)) GOTO 121
120    CONTINUE
121    KMID(IVE1,IEL)=KMID(JVE2,JEL)
       KNPR(KMID(IVE1,IEL))=KNPR(KMID(JVE2,JEL))
C
      ENDIF
C
110   CONTINUE
100   CONTINUE
C
      END


      SUBROUTINE S2MEL(KMEL,KMID,IMID,NEL,NMT,NVT)
C
C       The purpose of this routine is to determine the numbers of the
C       elements adjacent to each midpoint.  These numbers are stored on
C       KMEL(2,NMT)

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)

      DIMENSION KMEL(2,*),KMID(4,*)

      if(IMID.gt.0) then
        DO  IELNR=1,NEL
          N1=KMID(1,IELNR)-NVT
          N2=KMID(2,IELNR)-NVT
          N3=KMID(3,IELNR)-NVT
          N4=KMID(4,IELNR)-NVT
          IF (KMEL(1,N1).EQ.0) THEN
            KMEL(1,N1)=IELNR
          ELSE
            KMEL(2,N1)=IELNR
          ENDIF
          IF (KMEL(1,N2).EQ.0) THEN
            KMEL(1,N2)=IELNR
          ELSE
            KMEL(2,N2)=IELNR
          ENDIF
          IF (KMEL(1,N3).EQ.0) THEN
            KMEL(1,N3)=IELNR
          ELSE
            KMEL(2,N3)=IELNR
          ENDIF
          IF (KMEL(1,N4).EQ.0) THEN
            KMEL(1,N4)=IELNR
          ELSE
            KMEL(2,N4)=IELNR
          ENDIF
        ENDDO
      endif
      END
