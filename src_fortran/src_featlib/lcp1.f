************************************************************************
*                                                                      *
* XLCP1                                                                *
*                                                                      *
* Purpose  Copy a double precision vector                              *
*          Call  LCP1                                                  *
*                                                                      *
* Subroutines/functions called    ZTYPE, ZLEN, ZNEW, LCP1              *
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
* LCP1(A)                                                              *
*                                                                      *
* Purpose  Copy a vector                                               *
*          double precision version                                    *
*                                                                      *
* Subroutines/functions called   None                                  *
* BLAS                           DCOPY                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*X    Input vector                                         *
* NX       I*4    Number of elements to be copied                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DY       R*8    Copy of vector                                       *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLCP1(LX,LY,NX)
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
      IF (ICHECK.GE.998) CALL OTRC('XLCP1 ','01/02/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       IF (NX.GT.ILEN1.OR.NX.GT.ILEN2) THEN
        CALL WERR(-121,'XLCP1 ')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XLCP1 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LCP1(DWORK(L(LX)),DWORK(L(LY)),NX)
99999 END
C
C
C
      SUBROUTINE LCP1(DX,DY,NX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LCP1  ','01/02/89')
C
      CALL DCOPY(NX,DX,1,DY,1)
      END
