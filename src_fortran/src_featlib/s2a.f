************************************************************************
*                                                                      *
* XS2A                                                                 *
*                                                                      *
* Purpose  Allocate workspace and call S2A                             *
*                                                                      *
* Subroutines/functions called  S2A, ZNEW, ZDISP                       *
*                                                                      *
* Version from  05/20/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
*                                                                      *
************************************************************************
*                                                                      *
* S2A                                                                  *
*                                                                      *
* Purpose  Calculate KADJ from KVERT information                       *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  05/20/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* KVERT   I*4                                                          *
* KAUX1   I*4     Auxiliary vectors                                    *
* KAUX2   I*4                                                          *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* KADJ    I*4                                                          *
*                                                                      *
************************************************************************
C
      SUBROUTINE XS2A
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=999,NNVE=4)
      INCLUDE 'nnwork.inc'
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /TRIAD/,/TRIAA/,/OUTPUT/,/CHAR/,/ERRCTL/
C
      SUB='XS2A'
      IF (ICHECK.GE.995) CALL OTRC('XS2A  ','05/28/89')
      IER=0
C
      CALL ZNEW(NNVE*NEL,3,LADJ,'KADJ  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(3*NNVE*NEL,-3,LAUX1,'KAUX1 ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NVT,3,LAUX2,'KAUX2 ')
      IF (IER.NE.0) GOTO 99999
C
      CALL S2A(KWORK(L(LVERT)),KWORK(L(LADJ)),
     *         KWORK(L(LAUX1)),KWORK(L(LAUX2)))
C
      CALL ZDISP(0,LAUX2,'KAUX2 ')
      CALL ZDISP(0,LAUX1,'KAUX1 ')
C
99999 END
C
C
C
      SUBROUTINE S2A(KVERT,KADJ,KAUX1,KAUX2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNVE=4)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*)
      DIMENSION KAUX1(3,*),KAUX2(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE /TRIAD/,/OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='S2A'
      IF (ICHECK.GE.995) CALL OTRC('S2A   ','05/28/89')
      IER=0
C
C *** Initialize auxiliary arrays
C
C *** Determine number of occurences of each vertex
C
      DO 1 IEL=1,NEL
      DO 1 IVE=1,NVE
      IVT=KVERT(IVE,IEL)
1     KAUX2(IVT)=KAUX2(IVT)+1
      DO 2 IVT=2,NVT
2     KAUX2(IVT)=KAUX2(IVT-1)+KAUX2(IVT)
C
      DO 3 IEL=1,NEL
      DO 3 IVE=1,NVE
      IVT1=KVERT(IVE,IEL)
      IVE1=MOD(IVE,NVE)+1
      IPTR=KAUX2(IVT1)
      KAUX1(1,IPTR)=IVT1
      KAUX1(2,IPTR)=KVERT(IVE1,IEL)
      KAUX1(3,IPTR)=IEL
      KAUX2(IVT1)=IPTR-1
3     CONTINUE
      DO 4 IVT=1,NVT
4     KAUX2(IVT)=KAUX2(IVT)+1
C *** Search for interior edges and determine entries of KADJ
C
      DO 30 IEDGE=1,NVE*NEL
      IVT1=KAUX1(1,IEDGE)
      IVT2=KAUX1(2,IEDGE)
      IEL =KAUX1(3,IEDGE)
      DO 31 I=KAUX2(IVT2),NVE*NEL
      IF (KAUX1(1,I).NE.IVT2) GOTO 30
      IF (KAUX1(2,I).EQ.IVT1) THEN
C
C *** Interior edge found
C
       DO 34 IVE=1,NVE
       IF (KVERT(IVE,IEL).EQ.IVT1) GOTO 35
34     CONTINUE
35     KADJ(IVE,IEL)=KAUX1(3,I)
       GOTO 30
      ENDIF
C
31    CONTINUE
C
30    CONTINUE
C
99999 END
