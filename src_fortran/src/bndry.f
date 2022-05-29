************************************************************************
      SUBROUTINE BDRSET(DU,NVT,KNPR,DCORVG,UE)
************************************************************************
*    Purpose:  updates the solution vector DU and the right hand 
*              side DB for all Dirichlet boundary nodes
*-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NVT,KNPR(NVT),IVT
      DOUBLE PRECISION DU(NVT),DCORVG(2,NVT),X,Y,DUEX,UE
      EXTERNAL UE
C
      DO IVT=1,NVT
C
        IF (KNPR(IVT).NE.0) THEN
C
          X       = DCORVG(1,NVT)
          Y       = DCORVG(2,NVT)
          DUEX    = UE(X,Y)
C
          DU(IVT) = DUEX
C
        ENDIF
C
      ENDDO
C
      END
C
************************************************************************
      SUBROUTINE BDRSET0(DU,NVT,KNPR)
************************************************************************
*-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER NVT,KNPR(NVT),IVT
      DOUBLE PRECISION DU(NVT)
C
      DO IVT=1,NVT
C
        IF (KNPR(IVT).NE.0) THEN
C
          DU(IVT) = 0D0
C
        ENDIF
C
      ENDDO
C
      END
C
************************************************************************
      SUBROUTINE BDRYA(DA,KCOL,KLD,NVT,KNPR)
************************************************************************
*    Purpose:  updates the matrix entries for all Dirichlet boundary nodes
*-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER KCOL(*),KLD(*),NVT,KNPR(NVT),IVT,ICOL
      DOUBLE PRECISION DA(*)
C
      DO IVT=1,NVT
C
        IF (KNPR(IVT).NE.0) THEN
C
C         Set diagonal element to 1
          DA(KLD(IVT)) = 1D0
C
C         Zero other entries in row
          DO 2 ICOL=KLD(IVT)+1,KLD(IVT+1)-1
2           DA(ICOL) = 0D0
C
        ENDIF
C
      ENDDO
C
      END
