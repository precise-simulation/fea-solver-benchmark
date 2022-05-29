************************************************************************
*                                                                      *
* S2DI0                                                                *
*                                                                      *
* Purpose   Calculates the coordinates of the midpoint of a line given *
*           by the coordinates of the endpoints                        *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* X1,Y1    R*8    Coordinates of the endpoints of the line             *
* X2,Y2    R*8                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* XT,YT    R*8    Coordinates of the midpoint of the line              *
*                                                                      *
************************************************************************
C
      SUBROUTINE S2DI0(X1,Y1,X2,Y2,XT,YT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('S2DI0 ','01/02/89')
C
      XT=.5D0*(X1+X2)
      YT=.5D0*(Y1+Y2)
      END
