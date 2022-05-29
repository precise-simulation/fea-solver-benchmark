************************************************************************
*                                                                      *
* IA21n                                                                *
*                                                                      *
* Purpose  Smoothing using (damped) Jacobi-method                      *
*          Double/double precision version                             *
*                                                                      *
* Subroutines/functions called  LCP1, LLC1, LVM1                       *
*                                                                      *
* Version from  10/26/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix DA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* DD       R*8    Help vector                                          *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Number of smoothing steps                            *
* OMEGA    R*8    Damping factor                                       *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Smoothed vector                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE IA213(DA,KDIA,KDIAS,NDIA,DX,DB,DD,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KDIA(*),KDIAS(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA213 ','11/19/90')
C
      DO 1 IEQ=1,NEQ
1     DA(IEQ)=1D0/DA(IEQ)
C
      DO 2 ITE=1,NIT
      CALL LCP1(DB,DD,NEQ)
      DO 3 IDIA=2,NDIA
      J1=KDIA(IDIA)
      IF (J1.GT.0) THEN
       I1=1
       NEQ1=NEQ-J1
      ELSE
       I1=1-J1
       NEQ1=NEQ+J1
      ENDIF
      J0=KDIAS(IDIA)-I1
      CALL LVM1(DX(I1+J1),DA(I1+J0),DD(I1),NEQ1,-1D0,1D0)
3     CONTINUE
      CALL LVM1(DD,DA,DD,NEQ,1D0,0D0)
      CALL LLC1(DD,DX,NEQ,1D0-OMEGA,OMEGA)
2     CONTINUE
C
      DO 10 IEQ=1,NEQ
10    DA(IEQ)=1D0/DA(IEQ)
C
      END
C
C
C
      SUBROUTINE IA217(DA,KCOL,KLD,DX,DB,DD,NEQ,NIT,OMEGA)
      IMPLICIT NONE
      INTEGER NEQ,NIT,ITE,KCOL(*),KLD(NEQ+1),IER,ICHECK,IEQ,ICOL
      INTEGER M,MP1
      DOUBLE PRECISION DA(*),DX(NEQ),DB(NEQ),DD(NEQ),OMEGA
      COMMON /ERRCTL/ IER,ICHECK
C
      IF ( ICHECK.GE.998 ) CALL OTRC('IA217 ','10/26/89')
C
      DO 1 ITE=1,NIT
C
        CALL LCP1(DB,DD,NEQ)
C
!$OMP   PARALLEL DO
!$OMP*  PRIVATE(IEQ,ICOL)
!$OMP*  SHARED(DD,DA,DX,KLD,KCOL,NEQ)
        DO 2 IEQ=1,NEQ
          DO 3 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
            DD(IEQ) = DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL))
3         CONTINUE
2       CONTINUE
!$OMP   END PARALLEL DO
C
        IF ( OMEGA.EQ.1D0 ) THEN
          M = MOD(NEQ,7)
          IF ( M.EQ.0 ) GOTO 40
          DO 30 IEQ=1,M
            DX(IEQ) = DD(IEQ)/DA(KLD(IEQ))
30        CONTINUE
40        MP1 = M+1
!$OMP     PARALLEL DO
!$OMP*    PRIVATE(IEQ)
!$OMP*    SHARED(DD,DA,DX,KLD,MP1,NEQ)
          DO 50 IEQ=MP1,NEQ,7
            DX(IEQ  ) = DD(IEQ  )/DA(KLD(IEQ  ))
            DX(IEQ+1) = DD(IEQ+1)/DA(KLD(IEQ+1))
            DX(IEQ+2) = DD(IEQ+2)/DA(KLD(IEQ+2))
            DX(IEQ+3) = DD(IEQ+3)/DA(KLD(IEQ+3))
            DX(IEQ+4) = DD(IEQ+4)/DA(KLD(IEQ+4))
            DX(IEQ+5) = DD(IEQ+5)/DA(KLD(IEQ+5))
            DX(IEQ+6) = DD(IEQ+6)/DA(KLD(IEQ+6))
50        CONTINUE
!$OMP     END PARALLEL DO
c$$$          DO 4 IEQ=1,NEQ
c$$$            DX(IEQ) = DD(IEQ)/DA(KLD(IEQ))
c$$$4         CONTINUE
        ELSE
          DO 5 IEQ=1,NEQ
            DX(IEQ) = (1D0-OMEGA)*DX(IEQ)+OMEGA*DD(IEQ)/DA(KLD(IEQ))
5         CONTINUE
        ENDIF
C
1     CONTINUE
C
      END
C
C
C
      SUBROUTINE IA21A(DA,KCOL,KLD,KOP,DX,DB,DD,NEQ,NIT,OMEGA)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),KOP(*),DX(*),DB(*),DD(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IA21A ','10/26/89')
C
      DO 1 ITE=1,NIT
      CALL LCP1(DB,DD,NEQ)
      DO 2 IEQ=1,NEQ
      JOP=KOP(IEQ)
      DO 3 ICOL=KLD(JOP)+1,KLD(JOP+1)-1
3     DD(IEQ)=DD(IEQ)-DA(ICOL)*DX(KCOL(ICOL)+IEQ)
2     CONTINUE
      DO 4 IEQ=1,NEQ
4     DX(IEQ)=(1D0-OMEGA)*DX(IEQ)+OMEGA*DD(IEQ)/DA(KLD(KOP(IEQ)))
1     CONTINUE
C
      END
