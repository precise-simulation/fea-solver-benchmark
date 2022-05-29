************************************************************************
*                                                                      *
* XLSC1                                                                *
*                                                                      *
* Purpose  Scaling of a double precision vector                        *
*                                                                      *
* Subroutines/functions called  LSC1                                   *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX       I*4    Number of vector                                     *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* LSC1                                                                 *
*                                                                      *
* Purpose  Scaling of a double precision vector                        *
*                                                                      *
*                                                                      *
* Subroutines/functions called   LCL1                                  *
* BLAS                           DSCAL                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Vector                                               *
* NX       I*4    Length                                               *
* A        R*8    Scaling factor                                       *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Scaled vector                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLSC1(LX,NX,A)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /ERRCTL/,/CHAR/
C
      IF (ICHECK.GE.998) CALL OTRC('XLSC1 ','01/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE)
       CALL ZLEN(LX,ILEN)
       IF (NX.GT.ILEN) THEN
        CALL WERR(-121,'XLSC1 ')
        GOTO 99999
       ELSE IF (ITYPE.NE.1) THEN
        CALL WERR(-170,'XLSC1 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LSC1(DWORK(L(LX)),NX,A)
C
99999 END
C
C
C
      SUBROUTINE LSC1(DX,NX,A)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LSC1  ','01/05/89')
C
      IF (A.NE.1D0) THEN
       IF (A.EQ.0D0) THEN
        CALL LCL1(DX,NX)
       ELSE IF (A.EQ.-1D0) THEN
        DO 10 IX=1,NX
10      DX(IX)=-DX(IX)
       ELSE
        CALL DSCAL(NX,A,DX,1)
       ENDIF
      ENDIF
C
      END
