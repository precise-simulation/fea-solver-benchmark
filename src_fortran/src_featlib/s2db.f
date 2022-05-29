************************************************************************
*                                                                      *
* S2DB0                                                                *
*                                                                      *
* Purpose   Calculates the parameter value of the midpoint on the      *
*           boundary segment joining two vertices with given parameter *
*           values                                                     *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* X1       R*8    Parameter values                                     *
* X2       R*8                                                         *
* BMM      LOG    =.TRUE. if X1, X2 are minimum and maximum parameter  *
*                         values on boundary component IBCT            *
* IBCT            Boundary component                                   *
* PARX     R*8    Parametrization of the domain                        *
* PARY     R*8                                                         *
* TTMAX    R*8    Offset needed if BMM=.TRUE.                          *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* XT       R*8    Parameter value of midpoint                          *
*                                                                      *
************************************************************************
C
      SUBROUTINE S2DB0(X1,X2,XT,BMM,IBCT,PARX,PARY,TTMAX)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('S2DB0 ','01/02/89')
C
      IF (BMM) THEN
       XT=.5D0*(X1+X2+TTMAX)
       IF (XT.GT.TTMAX) XT=XT-TTMAX
      ELSE
       XT=.5D0*(X1+X2)
      ENDIF
C
      END
