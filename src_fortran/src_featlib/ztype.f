************************************************************************
*                                                                      *
* ZTYPE                                                                *
*                                                                      *
* Purpose  Determine datatype of an array allocated on DWORK           *
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
* LTYPE    I*4    >0  Datatype of array                                *
*                  0  Array not yet allocated                          *
*                 -1  LNR invalid                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZTYPE(LNR,LTYPE)
C
      PARAMETER (NNARR=999)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TABLE/  KTYPE(NNARR),KLEN(NNARR),KLEN8(NNARR),IFLAG
      SAVE /TABLE/
C
      IF (ICHECK.EQ.999) CALL OTRC('ZTYPE ','01/02/89')
C
      IF (LNR.GT.0.AND.LNR.LE.NNARR) THEN
       LTYPE=KTYPE(LNR)
      ELSE IF (LNR.EQ.0) THEN
       LTYPE=0
      ELSE
       LTYPE=-1
      ENDIF
C
      END
