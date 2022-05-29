************************************************************************
*                                                                      *
* ZNEW                                                                 *
*                                                                      *
* Purpose  Allocation of arrays on workspace DWORK                     *
*                                                                      *
* Subroutines/functions called   LCL1, LCL2, LCL3                      *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ILONG    I*4    Length of array to be allocated                      *
*                 For ILONG=0 the free part of DWORK                   *
*                 is completely reserved                               *
* ITYPE    I*4    Data type of array                                   *
*                  1  REAL*8     (Double precision)                    *
*                  2  REAL*4     (Single precision)                    *
*                  3  INTEGER*4                                        *
* ARR      C*6    Name of array - for messages only                    *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LNR      I*4    Number of array                                      *
*                 The address of the array is DWORK(L(LNR)),           *
*                 VWORK(L(LNR)), or KWORK(L(LNR))                      *
* ILONG    I*4    Set if ZNEW is called using ILONG=0                  *
*                 Returns number of free entries of type ITYPE         *
* IER      I*4    Error indicator                                      *
*                 -100  Not enough space, IWORK=NWORK                  *
*                       ZDISP must be called after ZNEW (ILONG=0)      *
*                 -101  Wrong value of ITYPE                           *
*                 -102  Tried to allocate more than  NNARR=999  arrays *
*                 -103  Not enough space, increase NWORK at least by N *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZNEW(ILONG,ITYPEA,LNR,ARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      DIMENSION VWORK(1),KWORK(1)
C
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TABLE/  KTYPE(NNARR),KLEN(NNARR),KLEN8(NNARR),IFLAG
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/,/TABLE/
C
      IF (ICHECK.GE.998) CALL OTRC('ZNEW  ','01/02/89')
      IER=0
C
      IF (IWORK.GE.NWORK) THEN
C ***  Error *** No space on DWORK ***
C ***  Call ZDISP first
       WRITE (CPARAM,'(A6)') ARR
       CALL WERR(-100,'ZNEW  ')
       GOTO 99999
      ENDIF
C
      ITYPE=ITYPEA
C
      IF ((ABS(ITYPE).LT.1).OR.(ABS(ITYPE).GT.3)) THEN
C ***  Error *** Wrong value of ITYPE ***
       WRITE (CPARAM,'(I15,A6)') ITYPE,ARR
       CALL WERR(-101,'ZNEW  ')
       GOTO 99999
      ENDIF
C
C *** Only allocation of new array if ITYPE < 0 - No initialization ***
C
      BZINIT=ITYPE.GT.0
      ITYPE=ABS(ITYPE)
C
C *** Block 1 - Allocation
C
C *** Look for a free entry in L ***
      DO 110 LNR=1,NNARR
      IF (L(LNR).EQ.0) GOTO 120
110   CONTINUE
C *** Error *** already NNARR arrays allocated ***
      WRITE (CPARAM,'(A6)') ARR
      CALL WERR(-102,'ZNEW  ')
      GOTO 99999
C
120   IF (ILONG.GT.0) THEN
C
C *** Calculation of space needed for allocation ***
      JLONGD=ILONG
      IF (ITYPE.GT.1) JLONGD=(JLONGD+1)/2
      JREQ=IWORK+JLONGD
C
      IF (JREQ.GT.NWORK) THEN
C ***  Error *** NWORK must be increased by JREQ-NWORK ***
       WRITE (CPARAM,'(I15,A6)') JREQ-NWORK,ARR
       CALL WERR(-103,'ZNEW  ')
       GOTO 99999
      ENDIF
C
C *** Set new value for KTYPE , KLEN , L , IWORK, IWMAX ***
      KTYPE(LNR)=ITYPE
      KLEN(LNR)=ILONG
      KLEN8(LNR)=JLONGD
      L(LNR)=IWORK+1
      IF (ITYPE.GT.1) L(LNR)=2*L(LNR)-1
      IWORK=JREQ
      IWMAX=MAX(IWMAX,IWORK)
C
C *** No initialization of new array ***
      IF (.NOT.BZINIT) GOTO 99998
C
C *** Block 2 - Initialization of new array by 0
C               according to data type
C
      GOTO (210,220,230),ITYPE
C *** REAL*8 data ***
210   CALL LCL1(DWORK(L(LNR)),KLEN(LNR))
      GOTO 99998
C *** REAL*4 data ***
220   CALL LCL2(VWORK(L(LNR)),KLEN(LNR))
      GOTO 99998
C *** INTEGER*4 data ***
230   CALL LCL3(KWORK(L(LNR)),KLEN(LNR))
      GOTO 99998
C
      ELSE
C
C *** ILONG = 0 ***
C
C *** Calculate number of free elements on DWORK ***
      ILONG=NWORK-IWORK
      KLEN8(LNR)=ILONG
      KTYPE(LNR)=ITYPE
      L(LNR)=IWORK+1
      IF (ITYPE.GT.1) THEN
       ILONG=ILONG*2
       L(LNR)=2*L(LNR)-1
      ENDIF
      KLEN(LNR)=ILONG
C *** reserve remaining part of DWORK - no initialization ***
      IFLAG=LNR
      IWORK=NWORK
C
      ENDIF
C
99998 WRITE (CPARAM,'(A6,2I15)') ARR,LNR,ILONG
      IF (IFLAG.NE.0) THEN
       CALL OMSG(2,'ZNEW  ')
      ELSE
       CALL OMSG(1,'ZNEW  ')
      ENDIF
C
99999 END
