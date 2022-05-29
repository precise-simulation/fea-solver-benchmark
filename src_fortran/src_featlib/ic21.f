************************************************************************
*                                                                      *
* IC21n                                                                *
*                                                                      *
* Purpose  Smoothing using SOR-method                                  *
*          Double/double precision version                             *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  10/27/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix DA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Maximum number of smoothing steps                    *
* OMEGA    R*8    Relaxation parameter                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IC217(DA,KCOL,KLD,DX,DB,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IC217 ','10/27/89')
C
      DO 1 ITE=1,NIT
         DO 2 IEQ=1,NEQ
      AUX=DB(IEQ)
      DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
3     AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
      AUX=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
2     DX(IEQ)=AUX
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IC21A(DA,KCOL,KLD,KOP,DX,DB,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),KOP(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.997) CALL OTRC('IC21A ','10/27/89')
C
      DO 1 ITE=1,NIT
      DO 2 IEQ=1,NEQ
      AUX=DB(IEQ)
      JOP=KOP(IEQ)
      DO 3 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
3     AUX=AUX-DA(ICOL)*DX(KCOL(ICOL)+IEQ)
      AUX=OMEGA*(AUX/DA(KLD(JOP))-DX(IEQ))+DX(IEQ)
2     DX(IEQ)=AUX
1     CONTINUE
C
      END
