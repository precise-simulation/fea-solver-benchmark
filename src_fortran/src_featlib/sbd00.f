************************************************************************
*                                                                      *
* SBD01, SBD02, SBD03                                                  *
*                                                                      *
* Purpose  Auxiliary routines to exchange information on boundary      *
*          vertices between DCORVG and DVBDP (DMBDP)                   *
*                                                                      *
*                                                                      *
* SBD01: Determine KVBD (without sorting), copy parameter values       *
*            from DCORVG to DVBDP (and DMBDP)                          *
* SBD02: Inverse of SBD01                                              *
* SBD03: Same as SBD01, but KVBD is already available and the order    *
*            is preserved                                              *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  08/15/90                                               *
*                                                                      *
* I/O     TYPE                                                         *
* -----   ----                                                         *
* DCORVG  R*8    See description of triangulation                      *
* KNPR    I*4                                                          *
* KVBD    I*4                                                          *
* DVBDP   R*8                                                          *
* DMBDP   R*8                                                          *
*                                                                      *
************************************************************************
*                                                                      *
* SBD04                                                                *
*                                                                      *
* Purpose  Replace parameter values for boundary vertices by           *
*          cartesian coordinates                                       *
*                                                                      *
* Subroutines/functions called  none                                   *
*                                                                      *
* Version from  08/15/90                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----   ----                                                         *
* DCORVG  R*8    See description of triangulation                      *
* KNPR    I*4                                                          *
* PARX    R*8                                                          *
* PARY    R*8                                                          *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* DCORVG  R*8    See description of triangulation                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE SBD01(DCORVG,KNPR,KVBD,DVBDP,DMBDP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DCORVG(2,*),KNPR(*),KVBD(*),DVBDP(*),DMBDP(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBD01'
      IF (ICHECK.GE.997) CALL OTRC('SBD01 ','08/15/90')
C
      IVBD=0
      DO 10 IVT=1,NVT
      IF (KNPR(IVT).GT.0) THEN
       IVBD=IVBD+1
       KVBD(IVBD)=IVT
       DVBDP(IVBD)=DCORVG(1,IVT)
       IF (NMT.GT.0) DMBDP(IVBD)=DCORVG(2,IVT)
      ENDIF
10    CONTINUE
C
      END



      SUBROUTINE SBD02(KVBD,DVBDP,DMBDP,DCORVG)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DCORVG(2,*),KVBD(*),DVBDP(*),DMBDP(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBD02'
      IF (ICHECK.GE.997) CALL OTRC('SBD02 ','08/15/90')
C
      DO 10 IVBD=1,NVBD
      IVT=KVBD(IVBD)
      DCORVG(1,IVT)=DVBDP(IVBD)
      IF (NMT.GT.0) THEN
       DCORVG(2,IVT)=DMBDP(IVBD)
      ELSE
       DCORVG(2,IVT)=0D0
      ENDIF
10    CONTINUE
C
      END




      SUBROUTINE SBD03(DCORVG,KVBD,DVBDP,DMBDP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DCORVG(2,*),KVBD(*),DVBDP(*),DMBDP(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBD03'
      IF (ICHECK.GE.997) CALL OTRC('SBD03 ','08/15/90')
C
      DO 10 IVBD=1,NVBD
      IVT=KVBD(IVBD)
      DVBDP(IVBD)=DCORVG(1,IVT)
      IF (NMT.GT.0) DMBDP(IVBD)=DCORVG(2,IVT)
10    CONTINUE
C
      END




      SUBROUTINE SBD04(DCORVG,KNPR,PARX,PARY)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DCORVG(2,*),KNPR(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='SBD04'
      IF (ICHECK.GE.997) CALL OTRC('SBD04 ','08/15/90')
C
      DO 10 IVT=1,NVT
      INPR=KNPR(IVT)
      IF (INPR.GT.0) THEN
       DCORVG(2,IVT)=PARY(DCORVG(1,IVT),INPR)
       DCORVG(1,IVT)=PARX(DCORVG(1,IVT),INPR)
      ENDIF
10    CONTINUE
C
      END
