************************************************************************
*                                                                      *
* XUSQGRD - Create a grid on a square with dimensions LX, LY           *
*                                                                      *
************************************************************************
C
      SUBROUTINE XUSQGRD(NX,NY,LX,LY,LCORVG,NVT,NMT,LVERT,LNPR,NEL,IPAR)
      IMPLICIT NONE
C
      INTEGER NX,NY,LCORVG,NVT,NMT,LVERT,LNPR,NEL,IPAR
      DOUBLE PRECISION LX,LY
C
      INTEGER NNARR,NWORK,IWORK,IWMAX,L
      INTEGER          KWORK(1)
      REAL*4           VWORK(1)
      INCLUDE 'nnwork.inc'
      DOUBLE PRECISION DWORK
      PARAMETER (NNARR=999)
      COMMON    NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      NVT = (NX+1)*(NY+1)
      NMT = NX*(NY+1)+NY*(NX+1)
      NEL = NX*NY
      CALL ZNEW(2*NVT,1,LCORVG,'DCORVG')
      CALL ZNEW(4*NEL,3,LVERT ,'KVERT ')
      CALL ZNEW(NVT,  3,LNPR  ,'KNPR  ')
C
      CALL USQGRD(NX,NY,LX,LY,DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *            KWORK(L(LNPR)),IPAR)
C
      END
C
************************************************************************
C
      SUBROUTINE USQGRD(NX,NY,LX,LY,DCORVG,KVERT,KNPR,IPAR)
      IMPLICIT NONE
C
      INTEGER NX,NY,KVERT(4,*),IROW,ICOL,IVT,IEL,KNPR(*),IPAR
      DOUBLE PRECISION LX,LY,DCORVG(2,*),DX,DY,X,Y
C
      DX = LX/DBLE(NX)
      DY = LY/DBLE(NY)
      Y  = 0D0
C
C     Compute vertex coordinates.
      DO IROW=1,NY+1
C
        X = 0D0
C
        DO ICOL=1,NX+1
C
          IVT = ICOL+(NX+1)*(IROW-1)
C
          DCORVG(1,IVT) = X
          DCORVG(2,IVT) = Y
          X             = X+DX
C
C         Assign cell connectivities.
          IF (IROW.LT.NY+1.AND.ICOL.LT.NX+1) THEN
            IEL          = ICOL+NX*(IROW-1)
            KVERT(1,IEL) = IVT
          ENDIF
          IF (IROW.LT.NY+1.AND.ICOL.GT.1) THEN
            IEL          = ICOL-1+NX*(IROW-1)
            KVERT(2,IEL) = IVT
          ENDIF
          IF (IROW.GT.1.AND.ICOL.GT.1) THEN
            IEL          = ICOL-1+NX*(IROW-2)
            KVERT(3,IEL) = IVT
          ENDIF
          IF (IROW.GT.1.AND.ICOL.LT.NX+1) THEN
            IEL          = ICOL+NX*(IROW-2)
            KVERT(4,IEL) = IVT
          ENDIF
C
        ENDDO
C
        Y = Y+DY
C
      ENDDO
C
C     Bottom boundary
      DO IVT=1,NX+1
        KNPR(IVT) = 1
        IF (IPAR.NE.0) THEN
c          DCORVG(1,IVT) = DCORVG(1,IVT)/LX
c          DCORVG(2,IVT) = 0D0
        ENDIF
      ENDDO
c      KMM(1,1) = 2
c      KMM(2,1) = NX
C
C     Right boundary
      DO IVT=NX+1,(NX+1)*(NY+1),NX+1
        KNPR(IVT) = 2
        IF (IPAR.NE.0) THEN
c          DCORVG(1,IVT) = DCORVG(2,IVT)/LY
c          DCORVG(2,IVT) = 0D0
        ENDIF
      ENDDO
c      KMM(1,2) = NX+1
c      KMM(2,2) = (NX+1)*(NY)
C
C     Top boundary
      DO IVT=(NX+1)*(NY+1)-NX,(NX+1)*(NY+1)
        KNPR(IVT) = 3
        IF (IPAR.NE.0) THEN
c          DCORVG(1,IVT) = DCORVG(1,IVT)/LX
c          DCORVG(2,IVT) = 0D0
        ENDIF
      ENDDO
c      KMM(1,3) = (NX+1)*(NY+1)-NX+1
c      KMM(2,3) = (NX+1)*(NY+1)
C
C     Left boundary
      DO IVT=1,(NX+1)*(NY+1)-NX,NX+1
        KNPR(IVT) = 4
        IF (IPAR.NE.0) THEN
c          DCORVG(1,IVT) = DCORVG(2,IVT)/LY
c          DCORVG(2,IVT) = 0D0
        ENDIF
      ENDDO
c      KMM(1,4) = 1
c      KMM(2,4) = (NX+1)*(NY+1)-NX
C
      END
C
C
C
************************************************************************
      DOUBLE PRECISION FUNCTION PARX(T,IBCT)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /GEOM/   DLX,DLY
C
      IF (IBCT.EQ.1.OR.IBCT.EQ.3) THEN
        PARX = T*DLX
      ELSE
        PARX = 0D0
      ENDIF
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION PARY(T,IBCT)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      COMMON /GEOM/   DLX,DLY
C
      IF (IBCT.EQ.2.OR.IBCT.EQ.4) THEN
        PARY = T*DLY
      ELSE
        PARY = 0D0
      ENDIF
C
      END
C
************************************************************************
      DOUBLE PRECISION FUNCTION TMAX(IBCT)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
C
      TMAX = 1D0
C
      END
