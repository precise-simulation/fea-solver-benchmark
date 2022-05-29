************************************************************************
*                                                                      *
* XLCP3                                                                *
*                                                                      *
* Purpose  Copy an integer vector                                      *
*          Call  LCP3                                                  *
*                                                                      *
* Subroutines/functions called    ZTYPE, ZLEN, ZNEW, LCP3              *
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
* LCP3                                                                 *
*                                                                      *
* Purpose  Copy a vector                                               *
*          Integer version                                             *
*                                                                      *
* Subroutines/functions called   None                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KX       I*4    Input vector                                         *         * NX       I*4    Number of elements to be copied                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLCP3(LX,LY,NX)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/CHAR/,/ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('XLCP3 ','01/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NX.GT.ILEN1.OR.NX.GT.ILEN2) THEN
        CALL WERR(-121,'XLCP3 ')
        GOTO 99999
       ELSE IF (ITYPE1.NE.3.OR.ITYPE2.NE.3) THEN
        CALL WERR(-170,'XLCP3 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LCP3(KWORK(L(LX)),KWORK(L(LY)),NX)
99999 END
C
C
C
      SUBROUTINE LCP3(KX,KY,NX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KX(*),KY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LCP3  ','01/02/89')
C
      DO 10 IX=1,NX
10    KY(IX)=KX(IX)
C
      END
