************************************************************************
*                                                                      *
* ZLEN                                                                 *
*                                                                      *
* Purpose  Determine number of elements of an array allocated on DWORK *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LNR      I*4    Number of array                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LENGTH   I*4    >0  Number of elements available                     *
*                  0  Array not yet allocated                          *
*                 -1  LNR invalid                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZLEN(LNR,LENGTH)
C
      PARAMETER (NNARR=999)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TABLE/  KTYPE(NNARR),KLEN(NNARR),KLEN8(NNARR),IFLAG
      SAVE /TABLE/
C
      IF (ICHECK.EQ.999) CALL OTRC('ZLEN  ','01/02/89')
C
      IF (LNR.GT.0.AND.LNR.LE.NNARR) THEN
       LENGTH=KLEN(LNR)
      ELSE IF (LNR.EQ.0) THEN
       LENGTH=0
      ELSE
       LENGTH=-1
      ENDIF
C
      END
