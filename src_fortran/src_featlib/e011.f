************************************************************************
*                                                                      *
* E011                                                                 *
*                                                                      *
* Purpose  Calculation of values and derivatives of element  11        *
*          Bilinear element                                            *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* XI1      R*8    Evaluation point in cartesian coordinates            *
* XI2      R*8    with respect to the reference square                 *
* IPAR     I*4    Switch for desired action                            *
*                  0  Calculation of the desired values of DBAS        *
*                 -1  Set number of element                            *
*                 -2  Calculate values on the reference element        *
*                     for all cubature points and save them            *
*                 -3  same as 0, but use the values saved before       *
*                                                                      *
* BDER     LOG    Derivative J is calculated if BDER(J).EQ..TRUE.      *
*                 Multiindices are enumbered from 1 to 6               *
* DJAC     R*8    Jacobian                                             *
* DETJ     R*8    Determinant of the jacobian                          *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DBAS     R*8    Array of function values and derivatives for         *
*                 all basis functions (block /ELEM/)                   *
* IPAR     I*4    Set to  11  if called with IPAR = -1                 *
* IER      I*4    Error indicator                                      *
*                 -131  Desired derivative not available               *
*                 -132  Triangle has vanishing area                    *
*                 -130  Clockwise ordering of the corner points        *
*                                                                      *
************************************************************************
C
      SUBROUTINE E011(XI1,XI2,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4,NDOFL=4)
      PARAMETER (Q4=0.25D0)
      DIMENSION DHELP(NDOFL,2)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      SAVE /OUTPUT/,/ERRCTL/,/ELEM/,/CHAR/,/CUB/
C
      SUB='E011'
      IF (ICHECK.GE.998) CALL OTRC('E011  ','01/02/89')
C
      IER=0
      GOTO (10,1,2,3) -IPAR+1
C
1     IPAR=11
2     GOTO 99999
C
3     CONTINUE
10    IF (ICHECK.EQ.0) GOTO 1100
      IF (DETJ.LT.0D0) CALL WERR(-130,'E011  ')
      IF (DETJ.GT.0D0.AND.DETJ.LT.1D-70) CALL WERR(-132,'E011  ')
C *** No second order derivatives available
C *** Used for second order problems only
      IF (BDER(4).OR.BDER(5).OR.BDER(6)) CALL WERR(-131,'E011  ')
      IF (IER.NE.0) GOTO 99999
C
1100  IF (.NOT.BDER(1)) GOTO 1101
C *** Function values
      DBAS(1,1)= Q4*(1D0-XI1)*(1D0-XI2)
      DBAS(2,1)= Q4*(1D0+XI1)*(1D0-XI2)
      DBAS(3,1)= Q4*(1D0+XI1)*(1D0+XI2)
      DBAS(4,1)= Q4*(1D0-XI1)*(1D0+XI2)
1101  IF (.NOT.(BDER(2).OR.BDER(3))) GOTO 99999
C
C *** First order derivatives
      XJ=Q4/DETJ
      DHELP(1,1)=-(1D0-XI2)
      DHELP(2,1)= (1D0-XI2)
      DHELP(3,1)= (1D0+XI2)
      DHELP(4,1)=-(1D0+XI2)
      DHELP(1,2)=-(1D0-XI1)
      DHELP(2,2)=-(1D0+XI1)
      DHELP(3,2)= (1D0+XI1)
      DHELP(4,2)= (1D0-XI1)
      IF (.NOT.BDER(2)) GOTO 1102
      DBAS(1,2)= XJ*(DJAC(2,2)*DHELP(1,1)-DJAC(2,1)*DHELP(1,2))
      DBAS(2,2)= XJ*(DJAC(2,2)*DHELP(2,1)-DJAC(2,1)*DHELP(2,2))
      DBAS(3,2)= XJ*(DJAC(2,2)*DHELP(3,1)-DJAC(2,1)*DHELP(3,2))
      DBAS(4,2)= XJ*(DJAC(2,2)*DHELP(4,1)-DJAC(2,1)*DHELP(4,2))
1102  IF (.NOT.BDER(3)) GOTO 99999
      DBAS(1,3)=-XJ*(DJAC(1,2)*DHELP(1,1)-DJAC(1,1)*DHELP(1,2))
      DBAS(2,3)=-XJ*(DJAC(1,2)*DHELP(2,1)-DJAC(1,1)*DHELP(2,2))
      DBAS(3,3)=-XJ*(DJAC(1,2)*DHELP(3,1)-DJAC(1,1)*DHELP(3,2))
      DBAS(4,3)=-XJ*(DJAC(1,2)*DHELP(4,1)-DJAC(1,1)*DHELP(4,2))
C
99999 END
