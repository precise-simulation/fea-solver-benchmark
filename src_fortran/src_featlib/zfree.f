************************************************************************
*                                                                      *
* ZFREE                                                                *
*                                                                      *
* Purpose  Calculate free space on workspace DWORK                     *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  08/14/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ITYPE    I*4    Data type                                            *
*                  1  REAL*8     (Double precision)                    *
*                  2  REAL*4     (Single precision)                    *
*                  3  INTEGER*4                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IFREE    I*4    Number of free entries of type ITYPE                 *
* IER      I*4    Error indicator                                      *
*                 -106  Wrong value of ITYPE                           *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZFREE(ITYPE,IFREE)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /CHAR/,/ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ZFREE ','08/14/90')
      IER=0
C
      IF (ITYPE.LT.1.OR.ITYPE.GT.3) THEN
C ***  Error *** Wrong value of ITYPE ***
       WRITE (CPARAM,'(I15)') ITYPE
       CALL WERR(-106,'ZFREE ')
       GOTO 99999
      ENDIF
C
C *** Calculate number of free elements on DWORK ***
      IFREE=NWORK-IWORK
      IF (ITYPE.GT.1) IFREE=IFREE*2
C
99999 END
