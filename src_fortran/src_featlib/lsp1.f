************************************************************************
*                                                                      *
* XLSP1                                                                *
*                                                                      *
* Purpose  Scalar product of two double precision vectors              *
*          Call  LSP1                                                  *
*                                                                      *
* Subroutines/functions called    ZTYPE, ZLEN, LSP1                    *
*                                                                      *
* Version from  06/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX,LY    I*4    Numbers of vectors                                   *
* NX       I*4    Number of elements to be copied                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                                                                      *
************************************************************************
*                                                                      *
* LSP1                                                                 *
*                                                                      *
* Purpose  Scalar product                                              *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called   None                                  *
* BLAS                           DDOT                                  *
*                                                                      *
* Version from  06/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX,DY    R*8    Input vectors                                        *
* NX       I*4    Number of elements                                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* SP       R*8    Scalar product                                       *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLSP1(LX,LY,NX,SP)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      CHARACTER SUB*6,FMT*15,CPARAM*120
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('XLSP1 ','06/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NX.GT.ILEN1.OR.NX.GT.ILEN2) THEN
        CALL WERR(-121,'XLSP1 ')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XLSP1 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LSP1(DWORK(L(LX)),DWORK(L(LY)),NX,SP)
99999 END
C
C
C
      SUBROUTINE LSP1(DX,DY,NXY,SP)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('LSP1  ','06/02/89')
C
      SP=DDOT(NXY,DX,1,DY,1)
      END
