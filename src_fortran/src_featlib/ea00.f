************************************************************************
*                                                                      *
* EA00                                                                 *
*                                                                      *
* Purpose  Returns number of element                                   *
*                                                                      *
* Subroutines/functions called   ELEnn                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ELE3     FUNC   Element used for NVE=3                               *
* ELE4     FUNC   Element used for NVE=4                               *
* NVE      I*8                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IELTYP   I*4    Number of element                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE EA00 (ELE3,ELE4,NVE,IELTYP)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /ERRCTL/ IER,ICHECK
C
      IF (ICHECK.GE.997) CALL OTRC('EA00  ','01/02/89')
C
      IELTYP=-1
      IF (NVE.EQ.3) THEN
       CALL ELE3(0D0,0D0,0D0,IELTYP)
      ELSE IF (NVE.EQ.4) THEN
       CALL ELE4(0D0,0D0,IELTYP)
      ENDIF
C
      END
