************************************************************************
*                                                                      *
* XLLI1                                                                *
*                                                                      *
* Purpose  Maximum norm of a double precision vector                   *
*                                                                      *
* Subroutines/functions called  LLI1                                   *
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
* LLI1                                                                 *
*                                                                      *
* Purpose  Maximum norm of a double precision vector                   *
*                                                                      *
*                                                                      *
* Subroutines/functions called   none                                  *
* BLAS                           IDAMAX                                *
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
* XNORM    R*8    Maximum norm                                         *
* IND      I*4    Index of maximum component                           *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLLI1(LX,NX,XNORM,IND)
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
      IF (ICHECK.GE.998) CALL OTRC('XLLI1 ','01/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE)
       CALL ZLEN(LX,ILEN)
       IF (NX.GT.ILEN) THEN
        CALL WERR(-121,'XLLI1 ')
        GOTO 99999
       ELSE IF (ITYPE.NE.1) THEN
        CALL WERR(-170,'XLLI1 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LLI1(DWORK(L(LX)),NX,XNORM,IND)
C
99999 END
C
C
C
      SUBROUTINE LLI1(DX,NX,XNORM,IND)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('LLI1  ','01/02/89')
C
      IND=IDAMAX(NX,DX,1)
      XNORM=ABS(DX(IND))
      END
