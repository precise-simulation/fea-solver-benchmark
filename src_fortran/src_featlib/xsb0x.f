************************************************************************
*                                                                      *
* XSB0X                                                                *
*                                                                      *
* Purpose  Call of the following subdivision routines                  *
*                                                                      *
* Subroutines/functions called  XSBC, XSB0, XS2M, XS2V, XSVEB, ZDISP   *
*                                                                      *
* Version from  10/27/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ICHECK   I*4    =0 Skip check of subdivision                         *
* NFINE    I*4    Number of regular refinements (SB0)                  *
* IMID     I*4    >=1 Determine numbers of midpoints and KMEL          *
*                 >=2 Determine coordinates of midpoints               *
* IADJ     I*4    >0 Determine KADJ                                    *
*                 =0 Release array KADJ on return                      *
*                    after determination of the new subdivision        *
* IVEL     I*4    =1 Determine numbers of elements meeting at each vert*
* IDISP    I*4    =1 Release free space on all arrays after completion *
* IBDP     I*4    =-1 Adjust only NVBD and KBCT, (CALL XSVEB(-1,TMAX)) *
*                    after completion                                  *
*                 >=0 Store numbers of boundary vertices on KVBD       *
*                    (CALL XSVEB(0,TMAX))                              *
*                 >=1 Calculate sorted arrays KVBD and KEBD            *
*                    (CALL XSVEB(1,TMAX))                              *
*                 >=2 Store parameter values for boundary vertices     *
*                     (and midpoints) on DVBDP and DMBDP               *
*                                                                      *
* For the description of the remaining parameters see SB0              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8   Cartesian coordinates of vertices                     *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSB0X(NFINE,IMID,IADJ,IVEL,IDISP,IBDP,
     *                 S2DI,S2DB,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /ERRCTL/ IER,ICHECK
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL S2DI,S2DB,PARX,PARY,TMAX
      SAVE /ERRCTL/,/TRIAA/,/TRIAD/
C
      IF (ICHECK.GE.997) CALL OTRC('XSB0X ','10/27/89')
C
      IF (LVBD.EQ.0.OR.LVBDP.EQ.0.OR.NVBD.EQ.0) THEN
       CALL WERR(-164,'XOWA  ')
       GOTO 99999
      ENDIF
C
      CALL SBD02(KWORK(L(LVBD)),DWORK(L(LVBDP)),DWORK(L(LMBDP)),
     *           DWORK(L(LCORVG)))
C
      IF (NMT.GT.0) THEN
       CALL ZDISP(0,LMBDP,'DMBDP ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
      CALL ZDISP(0,LVBDP,'DVBDP ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LVBD,'KVBD  ')
      IF (IER.NE.0) GOTO 99999
      IF (LEBD.NE.0) CALL ZDISP(0,LEBD,'KEBD  ')
      IF (IER.NE.0) GOTO 99999
      NVBD=0
C
      IF (ICHECK.GT.0) THEN
       CALL XS2C(PARX,PARY,TMAX)
       IF (IER.NE.0) GOTO 99999
      ENDIF

C
      CALL XSB0(NFINE,IDISP,S2DI,S2DB,PARX,PARY,TMAX)
      IF (IER.NE.0) GOTO 99999

C
      IF (IVEL.GT.0) THEN
       CALL XS2V
       IF (IER.NE.0) GOTO 99999
      ENDIF

C
      IF (IMID.GT.0) THEN
       CALL XS2M(IMID,IADJ,IDISP,S2DI,S2DB,PARX,PARY,TMAX)
       IF (IER.NE.0) GOTO 99999
      ENDIF

C
      IF (IADJ.EQ.0.AND.LADJ.GT.0) CALL ZDISP(0,LADJ,'KADJ  ')
C
      IBDP1=IBDP
      IF (IBDP1.EQ.2) IBDP1=1
      CALL XSVEB(IBDP1,TMAX)

C
      IF (IBDP.EQ.2) THEN
       CALL ZNEW(NVBD,-1,LVBDP,'DVBDP ')
       IF (IER.NE.0) GOTO 99999
       IF (NMT.GT.0) THEN
        CALL ZNEW(NVBD,-1,LMBDP,'DMBDP ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       CALL SBD03(DWORK(L(LCORVG)),KWORK(L(LVBD)),
     *            DWORK(L(LVBDP)),DWORK(L(LMBDP)))

      ENDIF
C
      CALL SBD04(DWORK(L(LCORVG)),KWORK(L(LNPR)),PARX,PARY)
C
99999 END
