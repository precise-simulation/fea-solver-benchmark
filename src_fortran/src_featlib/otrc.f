************************************************************************
*                                                                      *
* OTRC                                                                 *
*                                                                      *
* Purpose  Output to file MTRC                                         *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* SUB      C*6    Name of subroutine                                   *
* VER      C*8    Date of version                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE OTRC(SUB,VER)
C
      CHARACTER SUB*6,VER*8
      LOGICAL BTRC
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      SAVE /OUTPUT/
      SAVE BTRC
      DATA BTRC/.FALSE./
C
      IF (SUB(1:5).EQ.'ZINIT') THEN
       BTRC=VER(1:1).EQ.'T'
      ELSE IF (BTRC) THEN
       WRITE (MTRC,10000) SUB,VER
10000  FORMAT(' Name:  ',A,'   Version from:  ',A)
      ENDIF
C
      END
