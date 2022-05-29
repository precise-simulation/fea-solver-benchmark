************************************************************************
*                                                                      *
* IA11n                                                                *
*                                                                      *
* Purpose  Preconditioning by scaling                                  *
*          Matrix stored in technique  n  (see Reference Manual)       *
*          Double/Double precision version                             *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8    Matrix stored in technique  n                        *
* KLD      I*4    Pointer vectors corresponding to the                 *
* KOP      I*4    storage technique                                    *
* DX       R*8    Vector                                               *
* NEQ      I*4    Number of equations (length of DX)                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Resulting vector                                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE IA113(DA,DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
      DIMENSION DA(*),DX(*)
C
      IF (ICHECK.GE.998) CALL OTRC('IA113 ','11/19/90')
C
      DO 1 IEQ=1,NEQ
1     DX(IEQ)=DX(IEQ)/DA(IEQ)
C
      END
C
C
C
      SUBROUTINE IA117(DA,KLD,DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
      DIMENSION DA(*),KLD(*),DX(*)
C
      IF (ICHECK.GE.998) CALL OTRC('IA117 ','01/02/89')
C
      DO 1 IEQ=1,NEQ
1     DX(IEQ)=DX(IEQ)/DA(KLD(IEQ))
C
      END
C
C
C
      SUBROUTINE IA11A(DA,KLD,KOP,DX,NEQ)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
      DIMENSION DA(*),KLD(*),DX(*),KOP(*)
C
      IF (ICHECK.GE.998) CALL OTRC('IA11A ','01/02/89')
C
      DO 1 IEQ=1,NEQ
1     DX(IEQ)=DX(IEQ)/DA(KLD(KOP(IEQ)))
C
      END
