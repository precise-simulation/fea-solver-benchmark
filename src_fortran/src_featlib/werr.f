************************************************************************
*                                                                      *
* WERR                                                                 *
*                                                                      *
* Purpose  Set error parameter IER and call OERR                       *
*                                                                      *
* Subroutines/functions called   OERR                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IER0     I*4    Number of error                                      *
* SUB0     C*6    Name of calling routine                              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Number of error                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE WERR(IER0,SUB0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER*6 SUB0
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/ERRCTL/
C
      IER=IER0
      CALL OERR(IER,SUB0)
C
      END
