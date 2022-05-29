************************************************************************
*                                                                      *
* ZCPY                                                                 *
*                                                                      *
* Purpose  Duplicate array on DWORK                                    *
*          Call  LCP1, LCP2, or LCP3  depending on data type           *
*                                                                      *
* Subroutines/functions called    ZTYPE, ZLEN, ZNEW, LCP1, LCP2, LCP3  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LNR1     I*4    Number of source array                               *
* ARR1     C*6    Name of source array                                 *
* ARR2     C*6    Name of target array                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LNR2     I*4    Number of target array                               *
* IER      I*4    Error indicator                                      *
*                 -104  Wrong input parameter LNR1 or LNR2             *
*                 -107  Wrong type of target array ARR2                *
*                 -108  Target array ARR2 too short                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZCPY(LNR1,ARR1,LNR2,ARR2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR1*6,ARR2*6
      DIMENSION VWORK(1),KWORK(1)
C
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE     (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      IF (ICHECK.GE.998) CALL OTRC('ZCPY  ','01/02/89')
      IER=0
C
      CALL ZTYPE(LNR1,ITYPE)
      IF (ITYPE.LE.0) THEN
       WRITE (CPARAM,'(I15,A6)') LNR1,ARR1
       CALL WERR(-104,'ZCPY  ')
       GOTO 99999
      ENDIF
      CALL ZLEN(LNR1,ILONG)
C
      IF (LNR2.EQ.0) THEN
       CALL ZNEW(ILONG,ITYPE,LNR2,ARR2)
       IF (IER.NE.0) GOTO 99999
      ELSE
       CALL ZLEN(LNR2,ILEN2)
       CALL ZTYPE(LNR2,ITYPE2)
       IF (ILEN2.LE.0) THEN
        WRITE (CPARAM,'(I15,A6)') LNR2,ARR2
        CALL WERR(-104,'ZCPY  ')
        GOTO 99999
       ELSE IF (ITYPE.NE.ITYPE2) THEN
        WRITE (CPARAM,'(A6,I15,A6,I15)') ARR1,LNR1,ARR2,LNR2
        CALL WERR(-107,'ZCPY  ')
        GOTO 99999
       ELSE IF (ILEN2.LT.ILONG) THEN
        WRITE (CPARAM,'(A6,I15,A6,I15)') ARR1,LNR1,ARR2,LNR2
        CALL WERR(-108,'ZCPY  ')
        GOTO 99999
       ENDIF
      ENDIF
C
      GOTO (10,20,30) ITYPE
C
10    CALL LCP1(DWORK(L(LNR1)),DWORK(L(LNR2)),ILONG)
      GOTO 50
C
20    CALL LCP2(VWORK(L(LNR1)),VWORK(L(LNR2)),ILONG)
      GOTO 50
C
30    CALL LCP3(KWORK(L(LNR1)),KWORK(L(LNR2)),ILONG)
C
50    WRITE (CPARAM,'(A6,I15,A6,I15)') ARR1,LNR1,ARR2,LNR2
      CALL OMSG(5,'ZCPY  ')
C
99999 END
