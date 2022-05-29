************************************************************************
*                                                                      *
* XLCL3                                                                *
*                                                                      *
* Purpose  Clear arrays                                                *
*          Call  LCL3                                                  *
*                                                                      *
* Subroutines/functions called   ZTYPE, ZLEN, LCL3                     *
*                                                                      *
* Version from  04/19/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX       I*4    Number of array                                      *
* NX       I*4    Number of elements to be cleared                     *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 -104  Wrong input parameter LX                       *
*                                                                      *
************************************************************************
*                                                                      *
* LCL3                                                                 *
*                                                                      *
* Purpose  Clearing an integer vector                                  *
*                                                                      *
* Subroutines/functions called   None                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NX       I*4    Number of elements to be cleared                     *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KX       I*4    Cleared vector                                       *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLCL3(LX,NX)
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
      IF (ICHECK.GE.998) CALL OTRC('XLCL3 ','04/19/89')
      IER=0
C
      IF (ICHECK.GE.3) THEN
       CALL ZTYPE(LX,ITYPE)
       CALL ZLEN(LX,ILEN)
       IF (NX.GT.ILEN) THEN
        CALL WERR(-121,'XLCL3 ')
        GOTO 99999
       ELSE IF (ITYPE.NE.3) THEN
        CALL WERR(-170,'XLCL3 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LCL3(KWORK(L(LX)),NX)
99999 END
C
C
C
      SUBROUTINE LCL3(KX,NX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LCL3  ','01/02/89')
C
      DO 10 IX=1,NX
10    KX(IX)=0
C
      END
