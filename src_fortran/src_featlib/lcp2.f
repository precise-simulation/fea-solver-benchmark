************************************************************************
*                                                                      *
* XLCP2                                                                *
*                                                                      *
* Purpose  Copy a single precision vector                              *
*          Call  LCP2                                                  *
*                                                                      *
* Subroutines/functions called    ZTYPE, ZLEN, ZNEW, LCP2              *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX       I*4    Number of source array                               *
* LY       I*4    Number of target array                               *
* NX       I*4    Number of elements to be copied                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                                                                      *
************************************************************************
*                                                                      *
* LCP2                                                                 *
*                                                                      *
* Purpose  Copy a vector                                               *
*          single precision version                                    *
*                                                                      *
* Subroutines/functions called   None                                  *
* BLAS                           SCOPY                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* VX       R*4    Input vector                                         *                       *
* NX       I*4    Number of elements to be copied                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* VY       R*4    Copy of vector                                       *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLCP2(LX,LY,NX)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1))
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('XLCP2 ','01/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NX.GT.ILEN1.OR.NX.GT.ILEN2) THEN
        CALL WERR(-121,'XLCP2 ')
        GOTO 99999
       ELSE IF (ITYPE1.NE.2.OR.ITYPE2.NE.2) THEN
        CALL WERR(-170,'XLCP2 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LCP2(VWORK(L(LX)),VWORK(L(LY)),NX)
99999 END
C
C
C
      SUBROUTINE LCP2(VX,VY,NX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION VX(*),VY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LCP2  ','01/02/89')
C
      CALL SCOPY(NX,VX,1,VY,1)
      END
