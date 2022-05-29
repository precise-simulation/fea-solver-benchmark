************************************************************************
*                                                                      *
* ID11n                                                                *
*                                                                      *
* Purpose  SSOR Preconditioning                                        *
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
* KCOL     I*4    storage technique                                    *
* KOP      I*4                                                         *
* DX       R*8    Vector                                               *
* NEQ      I*4    Number of equations (length of DX)                   *
* OMEGA    R*8    Relaxation parameter                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Resulting vector                                     *
*                                                                      *3
************************************************************************
C
      SUBROUTINE ID117(DA,KCOL,KLD,DX,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID117 ','01/02/89')
C
      DO 1 IEQ=1,NEQ
      AUX=0D0
      ILD=KLD(IEQ)
      DO 4 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).GE.IEQ) GOTO 1
4     AUX=AUX+DA(ICOL)*DX(KCOL(ICOL))
1     DX(IEQ)=(DX(IEQ)-AUX*OMEGA)/DA(ILD)
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      ILD=KLD(IEQ)
      DO 12 ICOL=ILD+1,KLD(IEQ+1)-1
      IF (KCOL(ICOL).LE.IEQ) GOTO 12
      AUX=AUX+DA(ICOL)*DX(KCOL(ICOL))
12    CONTINUE
10    DX(IEQ)=DX(IEQ)-AUX*OMEGA/DA(ILD)
C
      END
C
C
C
      SUBROUTINE ID118(DA,KCOL,KLD,DX,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID118 ','01/02/89')
C
      DO 1 IEQ=1,NEQ
      ILD=KLD(IEQ)
      DX(IEQ)=DX(IEQ)*OMEGA/DA(ILD)
      DO 4 ICOL=ILD+1,KLD(IEQ+1)-1
4     DX(KCOL(ICOL))=DX(KCOL(ICOL))-DA(ICOL)*DX(IEQ)
1     CONTINUE
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      ILD=KLD(IEQ)
      DO 12 ICOL=ILD+1,KLD(IEQ+1)-1
12    AUX=AUX+DA(ICOL)*DX(KCOL(ICOL))
10    DX(IEQ)=DX(IEQ)-AUX*OMEGA/DA(ILD)
C
      END
C
C
C
      SUBROUTINE ID11A(DA,KCOL,KLD,KOP,DX,NEQ,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),KOP(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('ID11A ','01/02/89')
C
      DO 1 IEQ=1,NEQ
      AUX=0D0
      IOP=KOP(IEQ)
      ILD=KLD(IOP)
      DO 4 ICOL=ILD+1,KLD(IOP+1)-1
      IF (KCOL(ICOL).GT.0) GOTO 1
4     AUX=AUX+DA(ICOL)*DX(KCOL(ICOL)+IEQ)
1     DX(IEQ)=(DX(IEQ)-AUX*OMEGA)/DA(ILD)
C
      DO 10 IEQ=NEQ-1,1,-1
      AUX=0D0
      IOP=KOP(IEQ)
      ILD=KLD(IOP)
      DO 12 ICOL=KLD(IOP+1)-1,ILD+1,-1
      IF (KCOL(ICOL).LE.0) GOTO 10
      AUX=AUX+DA(ICOL)*DX(KCOL(ICOL)+IEQ)
12    CONTINUE
10    DX(IEQ)=DX(IEQ)-AUX*OMEGA/DA(ILD)
C
      END
