************************************************************************
*                                                                      *
* ZCLEAR                                                               *
*                                                                      *
* Purpose  Clear arrays                                                *
*          Call  LCL1, LCL2, or LCL3  depending on data type           *
*                                                                      *
* Subroutines/functions called   ZTYPE, ZLEN, LCL1, LCL2, LCL3         *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LNR      I*4    Number of array                                      *
* ARR      C*6    Name of array                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 -104  Wrong input parameter LNR                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZCLEAR(LNR,ARR)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      CHARACTER SUB*6,FMT*15,CPARAM*120,ARR*6
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ZCLEAR','01/02/89')
      IER=0
C
      CALL ZTYPE(LNR,ITYPE)
      CALL ZLEN(LNR,ILEN)
      GOTO (99999,10,20,30), ITYPE+1
C *** Wrong value of LNR ***
      WRITE (CPARAM,'(I15,A6)') LNR,ARR
      CALL WERR(-104,'ZCLEAR')
      GOTO 99999
C
10    CALL LCL1(DWORK(L(LNR)),ILEN)
      GOTO 99999
C
20    CALL LCL2(VWORK(L(LNR)),ILEN)
      GOTO 99999
C
30    CALL LCL3(KWORK(L(LNR)),ILEN)
C
99999 END
