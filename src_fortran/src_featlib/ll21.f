************************************************************************
*                                                                      *
* XLL21                                                                *
*                                                                      *
* Purpose  l2-norm of a double precision vector                        *
*                                                                      *
* Subroutines/functions called  LL21                                   *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX       I*4    Number of vector                                     *
* NX       I*4    length                                               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
*                                                                      *
* LL21                                                                 *
*                                                                      *
* Purpose  l2-norm of a double precision vector                        *
*                                                                      *
*                                                                      *
* Subroutines/functions called   LSP1                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Vector                                               *
* NX       I*4    Length                                               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* XNORM    R*8    l2-norm                                              *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLL21(LX,NX,XNORM)
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
      IF (ICHECK.GE.998) CALL OTRC('XLL21 ','01/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE)
       CALL ZLEN(LX,ILEN)
       IF (NX.GT.ILEN) THEN
        CALL WERR(-121,'XLL21 ')
        GOTO 99999
       ELSE IF (ITYPE.NE.1) THEN
        CALL WERR(-170,'XLL21 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LL21(DWORK(L(LX)),NX,XNORM)
C
99999 END
C
C
C
      SUBROUTINE LL21(DX,NX,XNORM)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('LL21  ','01/02/89')
C
      CALL LSP1(DX,DX,NX,XNORM)
      XNORM=SQRT(XNORM)
      END
