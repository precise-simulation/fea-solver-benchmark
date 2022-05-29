************************************************************************
*                                                                      *
* XLLC1                                                                *
*                                                                      *
* Purpose  Linear combination of two double precision vectors          *
*          Call  LLC1                                                  *
*                                                                      *
* Subroutines/functions called   ZLEN, ZTYPE, LLC1                     *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX,LY    I*4    Number of input vectors                              *
* NX       I*4    Length of vectors                                    *
* A1,A2    R*8    DY := A1*DX + A2*DY                                  *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 -121 At least one of the vectors too small           *
*                      or invalid number                               *
*                 -170 Data type not single or double precision        *
*                                                                      *
************************************************************************
*                                                                      *
* LLC1                                                                 *
*                                                                      *
* Purpose  Linear combination of two double precision vectors          *
*                                                                      *
* Subroutines/functions called   None                                  *
* BLAS                           DAXPY,DCOPY,DSCAL                     *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX,DY    R*8    Input vectors                                        *
* NX       I*4    Length of vectors                                    *
* A1,A2    R*8    DY := A1*DX + A2*DY                                  *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DY       R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE XLLC1(LX,LY,NX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='XLLC1'
      IF (ICHECK.GE.998) CALL OTRC('XLLC1 ','01/02/89')
C
      IF (ICHECK.GE.3) THEN
       CALL ZLEN(LX,ILEN1)
       CALL ZLEN(LY,ILEN2)
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LY,ITYPE2)
       IF (NX.GT.ILEN1.OR.NX.GT.ILEN2) THEN
        CALL WERR(-121,'XLLC1 ')
        GOTO 99999
       ELSE IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR (-170,'XLLC1 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      CALL LLC1(DWORK(L(LX)),DWORK(L(LY)),NX,A1,A2)
C
99999 END
C
C
C
      SUBROUTINE LLC1(DX,DY,NX,A1,A2)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DX(*),DY(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('LLC1  ','01/02/89')
C
      IF (A2.EQ.0D0) THEN
       CALL DCOPY(NX,DX,1,DY,1)
       IF (A1.NE.1D0) CALL DSCAL(NX,A1,DY,1)
      ELSE IF (A2.EQ.1D0) THEN
       CALL DAXPY(NX,A1,DX,1,DY,1)
      ELSE
       A3=A1/A2
       CALL DAXPY(NX,A3,DX,1,DY,1)
       CALL DSCAL(NX,A2,DY,1)
      ENDIF
C
      END
