************************************************************************
*                                                                      *
* XS2V                                                                 *
*                                                                      *
* Purpose  Allocation of array KVEL                                    *
*          Call of S2V                                                 *
*                                                                      *
* Subroutines/functions called  S2V, ZNEW, ZDISP, ZLEN                 *
*                                                                      *
* Version from  09/28/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* S2V                                                                  *
*                                                                      *
* Purpose  Determine numbers of elements meeting at each vertex        *
*                                                                      *
* Subroutines/functions called   None                                  *
*                                                                      *
* Version from  09/28/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KVERT    I*4    Numbers of vertices in each element, counterclockwise*
* KADJ     I*4    Number of neighbouring element, 0 for boundary edges *
* ICHK     I*4    See below                                            *
* NVE      I*4                                                         *
* NEL      I*4    Parameters of subdivision                            *
* NVT      I*4                                                         *
* NVEL     I*4    Maximum number of elements meeting at a vertex       *
*                 Set to 1 for ICHK=1                                  *
*                 Not changed for ICHK=0                               *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KVEL     I*8    ICHK=1:                                              *
*                  KVEL(IVT,1) = #Elements meeting at vertex IVT       *
*                 ICHK=0:                                              *
*                  KVEL(IVT,.) = Numbers of elements meeting at IVT    *
*                                sorted in counterclockwise sense      *
* NVEL     I*4    Maximum number of elements meeting at a vertex       *
*                                                                      *
************************************************************************
C
      SUBROUTINE XS2V
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XS2V'
      IF (ICHECK.GE.997) CALL OTRC('XS2V  ','09/28/90')
C
      CALL ZLEN(LVEL,ILEN)
      IF (ILEN.LT.NVT) THEN
       IF (ILEN.GT.0) CALL ZDISP(0,LVEL,'KVEL  ')
       CALL ZNEW(NVT,3,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZCLEAR(LVEL,'KVEL  ')
      ENDIF
C
      NVEL=1
      CALL S2V(KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LVEL)),
     *         NVE,NEL,NVT,NVEL,1)
C
      CALL ZLEN(LVEL,ILEN)
      IF (ILEN.LT.NVEL*NVT) THEN
       CALL ZDISP(0,LVEL,'KVEL  ')
       CALL ZNEW (NVEL*NVT,3,LVEL,'KVEL  ')
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZCLEAR(LVEL,'KVEL  ')
      ENDIF
C
      CALL S2V(KWORK(L(LVERT)),KWORK(L(LADJ)),KWORK(L(LVEL)),
     *         NVE,NEL,NVT,NVEL,0)
C
99999 END
C
C
C
      SUBROUTINE S2V(KVERT,KADJ,KVEL,NVE,NEL,NVT,NVEL,ICHK)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
C
      PARAMETER (NNVE=4)
      DIMENSION KVERT(NNVE,*),KADJ(NNVE,*),KVEL(NVEL,*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.997) CALL OTRC('S2V   ','09/28/90')
C
      IF (ICHK.GT.0) THEN
C
C *** Determine number of elements meeting at vertex IVT
       DO 11 IEL=1,NEL
       DO 11 IVE=1,NVE
       IVT=KVERT(IVE,IEL)
11     KVEL(1,IVT)=KVEL(1,IVT)+1
       NVEL1=0
       DO 12 IVT=1,NVT
12     NVEL1=MAX(NVEL1,KVEL(1,IVT))
       NVEL=NVEL1
C
       GOTO 99999
      ENDIF
C
C *** Find one element containing the vertex IVT
      DO 20 IEL=1,NEL
      DO 21 IVE=1,NVE
      IVT=KVERT(IVE,IEL)
      IF (KVEL(1,IVT).EQ.0) KVEL(1,IVT)=IEL
21    CONTINUE
20    CONTINUE
C
C *** Find the remaining elements containing vertex IVT
      DO 30 IVT=1,NVT
      IEL0=KVEL(1,IVT)
      IEL=IEL0
C
C *** Search in counterclockwise sense first
      DO 33 IVEL=1,NVEL
C
      DO 31 IVE=1,NVE
      IVT1=KVERT(IVE,IEL)
      IF (IVT.EQ.IVT1) GOTO 32
31    CONTINUE
C
32    IADJ=MOD(IVE+NVE-2,NVE)+1
      IEL1=KADJ(IADJ,IEL)
C *** IEL1.EQ.0 means the vertex IVT is on the boundary
      IF (IEL1.EQ.0) GOTO 35
C *** Whole ring completed
      IF (IEL1.EQ.IEL0) GOTO 30
      KVEL(IVEL+1,IVT)=IEL1
      IEL=IEL1
33    CONTINUE
C
C *** Search in clockwise sense (for boundary vertices only)
35    IEL=IEL0
      DO 36 IVEL=NVEL,1,-1
      DO 37 IVE=1,NVE
      IVT1=KVERT(IVE,IEL)
      IF (IVT.EQ.IVT1) GOTO 38
37    CONTINUE
C
38    IEL1=KADJ(IVE,IEL)
      IF (IEL1.EQ.0) GOTO 39
      KVEL(IVEL,IVT)=IEL1
      IEL=IEL1
36    CONTINUE
C
C *** Cyclic permutation of the elements in KVEL until the last element
C *** found in clockwise sense is in the first position
39    DO 40 I=1,NVEL
      IF (KVEL(1,IVT).EQ.IEL) GOTO 30
      IAUX=KVEL(NVEL,IVT)
      DO 41 IVEL=NVEL,2,-1
41    KVEL(IVEL,IVT)=KVEL(IVEL-1,IVT)
      KVEL(1,IVT)=IAUX
40    CONTINUE
C
30    CONTINUE
C
99999 END
