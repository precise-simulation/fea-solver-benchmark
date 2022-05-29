************************************************************************
*                                                                      *
* XM010L                                                               *
*                                                                      *
* Purpose  Allocate Workspace vector on DWORK                          *
*          Call M010L                                                  *
*                                                                      *
* Subroutines/functions called   ZNEW, ZDISP, ZCPY, ZLEN, M010L        *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* LX       I*4    Number of starting vector, length >= KNEQ(NLMAX)     *
* LB       I*4    Number of right hand side, length >= KNEQ(NLMAX)     *
* YDAX     SUBR                                                        *
* YPROL    SUBR                                                        *
* YREST    SUBR   External subroutines                                 *
* YPRSM    SUBR                                                        *
* YPOSM    SUBR                                                        *
* YEX      SUBR                                                        *
* YBC      SUBR                                                        *
* IDISP    I*4    =1: On return, vectors assigned to numbers           *
*                     LX  and  LB are truncated to length  KNEQ(NLMAX) *
* Hint:   Use  IDISP=0  if  M010  is called several times              *
*         with the same parameters and  NNWORK  is large enough        *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE XM010L(LX,LB,KNEQ,NIT,ITE,EPS,
     *                  YDAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC,IDISP,IREL,
     *                  TTSOL)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNLEV=19,NNARR=999)
      INCLUDE 'nnwork.inc'
      DIMENSION VWORK(1),KWORK(1)
      DIMENSION KNEQ(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL YDAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XM010L'
      IF (ICHECK.GE.997) CALL OTRC('XM010L','04/12/91')
      IER=0
C
      CALL ZLEN(LX,ILENX)
      CALL ZLEN(LB,ILENB)
      IF (ICHECK.GT.0) THEN
       IF (ICYCLE.LT.0) THEN
        CALL WERR(-181,'XM010 ')
        GOTO 99999
       ENDIF
       IF (NLMAX.LT.NLMIN.OR.NLMIN.LE.0.OR.NLMAX.GT.NNLEV) THEN
        CALL WERR(-182,'XM010 ')
        GOTO 99999
       ENDIF
       IF (ILENX.LT.KNEQ(NLMAX).OR.ILENB.LT.KNEQ(NLMAX)) THEN
        CALL WERR(-121,'XM010 ')
        GOTO 99999
       ENDIF
       CALL ZTYPE(LX,ITYPE1)
       CALL ZTYPE(LB,ITYPE2)
       IF (ITYPE1.NE.1.OR.ITYPE2.NE.1) THEN
        CALL WERR(-170,'XM010 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      IREQ=0
      DO 1 ILEV=NLMIN,NLMAX
1     IREQ=IREQ+KNEQ(ILEV)
      IF (ILENX.LT.IREQ) THEN
       CALL ZNEW(IREQ,-1,LX1,'DX    ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LX,'DXOLD ',LX1,'DX    ')
       CALL ZDISP(0,LX,'DXOLD ')
       IF (IER.NE.0) GOTO 99999
       LX=LX1
      ENDIF
      IF (ILENB.LT.IREQ) THEN
       CALL ZNEW(IREQ,-1,LB1,'DB    ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LB,'DBOLD ',LB1,'DB    ')
       CALL ZDISP(0,LB,'DBOLD ')
       IF (IER.NE.0) GOTO 99999
       LB=LB1
      ENDIF
      CALL ZNEW(IREQ,-1,LD,'DD    ')
      IF (IER.NE.0) GOTO 99999
C
      CALL ZNEW(NLMAX,-3,LOFFX,'KOFFX ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NLMAX,-3,LOFFB,'KOFFB ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NLMAX,-3,LOFFD,'KOFFD ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NLMAX,-3,LIT,'KIT   ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NLMAX,-3,LIT0,'KIT0  ')
      IF (IER.NE.0) GOTO 99999
C

      LLX=L(LX)-1
      LLB=L(LB)-1
      LLD=L(LD)-1
      DO 2 ILEV=NLMAX,NLMIN,-1
      KWORK(L(LOFFX)+ILEV-1)=LLX
      KWORK(L(LOFFB)+ILEV-1)=LLB
      KWORK(L(LOFFD)+ILEV-1)=LLD
      LLX=LLX+KNEQ(ILEV)
      LLB=LLB+KNEQ(ILEV)
2     LLD=LLD+KNEQ(ILEV)
C
      CALL M010L(DWORK(1),DWORK(1),DWORK(1),
     *           KWORK(L(LOFFX)),KWORK(L(LOFFB)),KWORK(L(LOFFD)),
     *           KNEQ,NIT,ITE,EPS,
     *           YDAX,YPROL,YREST,YPRSM,YPOSM,YEX,YBC,
     *           KWORK(L(LIT0)),KWORK(L(LIT)),IREL,TTSOL)
C
      IER1=IER
      CALL ZDISP(0,LIT0,'KIT0  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LIT,'KIT   ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LOFFD,'KOFFD ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LOFFB,'KOFFB ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LOFFX,'KOFFX ')
      IF (IER.NE.0) GOTO 99999
      CALL ZDISP(0,LD,'DD    ')
      IF (IER.NE.0) GOTO 99999
      IF (IDISP.EQ.1) THEN
       CALL ZDISP(KNEQ(NLMAX),LB,'DB    ')
       CALL ZDISP(KNEQ(NLMAX),LX,'DX    ')
      ENDIF
      IER=IER1
C
99999 END
C
************************************************************************
*                                                                      *
* M010L                                                                 *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          multigrid iteration                                         *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LL21 , LLC1                            *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting address of vectors containing the           *
* DB       R*8    solution and the right hand side, DD is used as      *
* DD       R*8    auxiliary vector only                                *
* KOFFX           The actual starting address of DX on level ILEV      *
* KOFFB           is DX(1+KOFFX(ILEV)) (analogously for DB and DD)     *
* KOFFD           Total space required for all vectors is              *
*                 KNEQ(NLMIN)+...+KNEQ(NLMAX)                          *
*                 DX(1+KOFFX(NLMAX)) contains initial solution         *
*                 DB(1+KOFFB(NLMAX)) contains right hand side          *
* KNEQ     I*4    Number of equations for all levels                   *
* NLMAX    I*4    Iteration uses levels NLMIN to NLMAX,                *
* NLMIN    I*4    NLMAX  is the finest level                           *
* NIT      I*4    Maximum number of iterations                         *
*                 Iteration completed after reaching the finest level  *
* EPS      R*8    Desired precision                                    *
*                 Stop if !!DEF!! < EPS                                *
* KPRSM    I*4    Number of pre -smoothing steps for all levels        *
* KPOSM    I*4    Number of post-smoothing steps for all levels        *
* ICYCLE   I*4    <0: special cycle types (not yet implemented)        *
*                 =0: F-Cycle                                          *
*                 =1: V-Cycle                                          *
*                 =2: W-Cycle                                          *
*                 >2: Cycle of higher order                            *
* DAX      SUBR   CALL DAX(DX,DAX,NEQ,A1,A2)                           *
*                 Returns  DAX := A1*A*DX+A2*DAX                       *
* DPROL    SUBR   DPROL(DX1,DX2)                                       *
*                 Returns  DX2 := Prolongation(DX1)  to higher level   *
* DREST    SUBR   DREST(DD1,DD2)                                       *
*                 Returns  DD2 := Restriction(DD1)  to lower level     *
* DPRSM    SUBR   DPRSM(DX,DB,DD,NEQ,NPRSM)                            *
*                 Returns  DX  after  NPRSM:=KPRSM(ILEV)               *
*                 pre-smoothing steps, DD  is used as auxiliary vector *
* DPOSM    SUBR   Same as above, used for post-smoothing               *
* DEX      SUBR   DEX(DX,DB,DD,NEQ)                                    *
*                 Returns exact solution                               *
* DBC      SUBR   DBC(DX,NEQ)                                          *
*                 Copies boundary data onto components of  DX          *
* KIT0     I*4    auxiliary vectors of length  NLMAX                   *
* KIT      I*4                                                         *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*4    Solution vector on DX(KOFFX(NLMAX))                  *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                                                                      *
************************************************************************
C
      SUBROUTINE M010L(DX,DB,DD,KOFFX,KOFFB,KOFFD,KNEQ,NIT,ITE,EPS,
     *                DAX,DPROL,DREST,DPRSM,DPOSM,DEX,DBC,KIT0,KIT,IREL,
     *                TTSOL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=999,NNLEV=19)
      DIMENSION DX(*),DB(*),DD(*),KOFFX(*),KOFFB(*),KOFFD(*)
      DIMENSION KNEQ(*),KIT0(*),KIT(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      SAVE /ERRCTL/,/CHAR/,/OUTPUT/,/MGTRD/,/MGPAR/,/MGTIME/
C
c$$$      ITMIN = 3
C
      SUB='M010L '
      IF (ICHECK.GE.997) CALL OTRC('M010L ','08/25/90')
      IER=0
C
      BREL=IREL.EQ.1
      BMSG2=M.GE.2.OR.MT.GE.2
      MT0=MT
      MT=0
C
      BTIME=IMTIME.GT.0
      IF (BTIME) THEN
       IF (IMTIME.EQ.1) THEN
        TTMG=0D0
        TTS=0D0
        TTE=0D0
        TTD=0D0
        TTP=0D0
        TTR=0D0
       ENDIF
       CALL ZTIME(TTMG0)
      ENDIF
C
      NIT0=MAX(ITE,0)
      IF (ICYCLE.LT.0) THEN
       CALL WERR(-181,'M010  ')
       GOTO 99999
      ENDIF
      IF (NLMAX.LT.NLMIN.OR.NLMIN.LE.0.OR.NLMAX.GT.NNLEV) THEN
       CALL WERR(-182,'M010  ')
       GOTO 99999
      ENDIF
      ILEV=NLMAX
C
C *** special case - only one level
      IF (NLMIN.EQ.NLMAX) THEN
       CALL LCP1(DB(1+KOFFB(NLMAX)),DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
       IF (BTIME) CALL ZTIME(TTE0)
       CALL DEX(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
     *          DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       IF (BTIME) THEN
        CALL ZTIME(TTE1)
        TTE=TTE+TTE1-TTE0
       ENDIF
       GOTO 99999
      ENDIF
C
C *** level counts
      KIT0(NLMAX)=1
      DO 2 ILEV=NLMIN+1,NLMAX-1
      IF (ICYCLE.EQ.0) THEN
       KIT0(ILEV)=2
      ELSE
       KIT0(ILEV)=ICYCLE
      ENDIF
2     CONTINUE
C
      CALL ZTIME(TTSOL0)
      IF (KPRSM(NLMAX).GT.0)
     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
     *            DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
C
      IF (BTIME) THEN
       CALL ZTIME(TTS1)
       TTS=TTS+TTS1-TTSOL0
      ENDIF
C
C
      CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
      IF (BTIME) CALL ZTIME(TTD0)
      CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *         -1D0,1D0)
C
      IF (BTIME) THEN
       CALL ZTIME(TTD1)
       TTD=TTD+TTD1-TTD0
      ENDIF
C
      CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
      DEFOLD=DEF
C *** FD  is considered as initial defect
C *** (after at least one pre-smoothing step)
      FD=DEF
      IF (DEF.LE.EPS.AND..NOT.BREL) THEN
       ITE=0
       GOTO 1000
      ENDIF
C
C
C *** Start multigrid iteration
C
      DO 100 ITE=1,NIT
C
C *** initialize level counts for all levels
      DO 101 ILEV=NLMIN,NLMAX
101   KIT(ILEV)=KIT0(ILEV)
C
      ILEV=NLMAX
C
110   IF (ILEV.NE.NLMIN) THEN
C
C *** Pre-smoothing
C
       IF (ILEV.NE.NLMAX) THEN
C ***  defect on finest level already available
C
C
        IF (BTIME) CALL ZTIME(TTS0)
        IF (KPRSM(ILEV).GT.0)
     *   CALL DPRSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *              DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPRSM(ILEV))
C
        IF (BTIME) THEN
         CALL ZTIME(TTS1)
         TTS=TTS+TTS1-TTS0
        ENDIF
C
C
        CALL LCP1(DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV))
        IF (BTIME) CALL ZTIME(TTD0)
        CALL DAX(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV),
     *           -1D0,1D0)
C
        IF (BTIME) THEN
         CALL ZTIME(TTD1)
         TTD=TTD+TTD1-TTD0
        ENDIF
C
       ENDIF
C
       ILEV=ILEV-1
C ***  restriction of defect
C
       IF (BTIME) CALL ZTIME(TTR0)
       CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)))
C
       IF (BTIME) THEN
        CALL ZTIME(TTR1)
        TTR=TTR+TTR1-TTR0
       ENDIF
C
C ***  choose zero as initial vector on lower level
       CALL LCL1(DX(1+KOFFX(ILEV)),KNEQ(ILEV))
       CALL DBC(DB(1+KOFFB(ILEV)),KNEQ(ILEV))
       GOTO 110
      ENDIF
C
C *** Exact solution on lowest level
      CALL LCP1(DB(1+KOFFB(NLMIN)),DX(1+KOFFX(NLMIN)),KNEQ(NLMIN))
      IF (BTIME) CALL ZTIME(TTE0)
      CALL DEX(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),DD(1+KOFFD(NLMIN)),
     *         KNEQ(NLMIN))
C
      IF (BTIME) THEN
       CALL ZTIME(TTE1)
       TTE=TTE+TTE1-TTE0
      ENDIF
C
C
130   IF (ILEV.NE.NLMAX) THEN
       ILEV=ILEV+1
C ***  DPROL  returns  DD:=PROL(DX)
C
C
       IF (BTIME) CALL ZTIME(TTP0)
       CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)))
C
       IF (BTIME) THEN
        CALL ZTIME(TTP1)
        TTP=TTP+TTP1-TTP0
       ENDIF
C
       CALL DBC(DD(1+KOFFD(ILEV)),KNEQ(ILEV))
       CALL LLC1(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV),1D0,1D0)
C
C ***  Post-smoothing
       IF (BTIME) CALL ZTIME(TTS0)
       IF (KPOSM(ILEV).GT.0)
     * CALL DPOSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),
     *            KNEQ(ILEV),KPOSM(ILEV))
C
       IF (BTIME) THEN
        CALL ZTIME(TTS1)
        TTS=TTS+TTS1-TTS0
       ENDIF
C
       KIT(ILEV)=KIT(ILEV)-1
       IF (KIT(ILEV).EQ.0) THEN
        IF (ICYCLE.EQ.0) THEN
         KIT(ILEV)=1
        ELSE
         KIT(ILEV)=KIT0(ILEV)
        ENDIF
        GOTO 130
       ELSE
        GOTO 110
       ENDIF
      ENDIF
C
C
      IF (BTIME) CALL ZTIME(TTS0)
      IF (KPRSM(NLMAX).GT.0)
     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
     *           DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
C
      IF (BTIME) THEN
       CALL ZTIME(TTS1)
       TTS=TTS+TTS1-TTS0
      ENDIF
C
C
C
      CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
      IF (BTIME) CALL ZTIME(TTD0)
      CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *         -1D0,1D0)
C
      IF (BTIME) THEN
       CALL ZTIME(TTD1)
       TTD=TTD+TTD1-TTD0
      ENDIF
C
      CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
      IF (BMSG2) THEN
       MT=MT0
       WRITE (CPARAM,'(I15,D25.16)') ITE,DEF
       CALL OMSG(73,'M010  ')
       MT=0
      ENDIF
      DEFOLD=DEF
c$$$      IF (ITE.GT.ITMIN) THEN
      IF (BREL) THEN
       IF (DEF.LE.FD*EPS.AND.ITE.GE.NIT0) GOTO 1000
      ELSE
       IF (DEF.LE.EPS) GOTO 1000
      ENDIF
c$$$      ENDIF
C
100   CONTINUE
C
      MT=MT0
      WRITE (CPARAM,'(I15,2D25.16)') NIT,DEF,DEF/FD
      CALL OMSG(71,'M010  ')
      CALL OMSG(72,'M010  ')
      IER=1
      GOTO 99999
C
1000  CONTINUE
      IER=0
      MT=MT0
      IF (FD.GE.1D-70) FD=DEF/FD
      WRITE (CPARAM,'(I15,2D25.16)') ITE,DEF,FD
      CALL OMSG(72,'M010  ')
      WRITE (CPARAM,'(D25.16)') FD**(1D0/DBLE(ITE))
      CALL OMSG(76,'M010  ')
C
99999 MT=MT0
      IF (BTIME) THEN
       CALL ZTIME(TTMG1)
       TTMG=TTMG+TTMG1-TTMG0
      ENDIF
      CALL ZTIME(TTSOL1)
      TTSOL = TTSOL1-TTSOL0
C
      END
C
C
C
************************************************************************
*   Perform matrix-vector-operation
************************************************************************
      SUBROUTINE YAXU(DX,DAX,NEQ,A1,A2)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19)
      INCLUDE 'nnwork.inc'
      DIMENSION DX(*),DAX(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLNPR(NNLEV),KLVERT(NNLEV),KLADJ(NNLEV)
C-----------------------------------------------------------------------
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
      CALL LAX17(DWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *           KWORK(L(KLLDA(ILEV))),NEQ,DX,DAX,A1,A2)
C
      END
C
C
C
************************************************************************
*   Set Dirichlet boundary components of DX to zero
************************************************************************
      SUBROUTINE YDBCU(DX,NEQ)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19)
      INCLUDE 'nnwork.inc'
      DIMENSION DX(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLNPR(NNLEV),KLVERT(NNLEV),KLADJ(NNLEV)
C-----------------------------------------------------------------------
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
      CALL BDRSET0(DX,NEQ,KWORK(L(KLNPR(ILEV))))
C
      END
C
C
C
************************************************************************
C     Coarse grid solver
************************************************************************
      SUBROUTINE YEXU(DX,DB,DD,NEQ)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19)
      INCLUDE 'nnwork.inc'
      DIMENSION DX(*),DB(*),DD(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGPAR2/ RLXSMF,DMPSLC,RLXSLC,ISM,ISOLC,NSLC,IREST,
     *                KNX(NNLEV),KNY(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLNPR(NNLEV),KLVERT(NNLEV),KLADJ(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
C-----------------------------------------------------------------------
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
      IF     (ISOLC.EQ.1) THEN   ! Jacobi
C
        CALL IA017(DWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),DX,DB,DD,
     *             NEQ,NSLC,ITE,DMPSLC,RLXSLC)
C
      ELSEIF (ISOLC.EQ.2) THEN   ! SOR(GS if RLXSLC=1)
C
        CALL IC017(DWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),DX,DB,
     *             NEQ,NSLC,ITE,DMPSLC,RLXSLC)
C
      ELSEIF (ISOLC.EQ.3) THEN   ! 4-Color SOR(GS if RLXSLC=1)
C
        CALL IC017_4C(DWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),DX,DB,
     *             NEQ,NSLC,ITE,DMPSLC,RLXSLC,KNX(ILEV)+1)
C
      ENDIF
C
      END
C
C
C
************************************************************************
C     Smoother
************************************************************************
      SUBROUTINE YSMU(DX,DB,DD,NEQ,NITSM)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19)
      INCLUDE 'nnwork.inc'
      DIMENSION DX(*),DB(*),DD(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGPAR2/ RLXSMF,DMPSLC,RLXSLC,ISM,ISOLC,NSLC,IREST,
     *                KNX(NNLEV),KNY(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLNPR(NNLEV),KLVERT(NNLEV),KLADJ(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
C-----------------------------------------------------------------------
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
      IF (ISM.EQ.1) THEN       ! Jacobi
C
        CALL IA217(DWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),DX,DB,DD,NEQ,NITSM,RLXSMF)
C
      ELSEIF (ISM.EQ.2) THEN   ! SOR(GS if RLXSLC=1)
C
        CALL IC217(DWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NITSM,RLXSMF)
C
      ELSEIF (ISM.EQ.3) THEN   ! 4-Color SOR(GS if RLXSLC=1)
C
        CALL IC217_4C(DWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NITSM,RLXSMF,
     *             KNX(ILEV)+1)
C
      ELSEIF (ISM.EQ.4) THEN   ! SSOR
C
        CALL ID217(DWORK(L(KLA(ILEV))),KWORK(L(KLCOLA(ILEV))),
     *             KWORK(L(KLLDA(ILEV))),DX,DB,NEQ,NITSM,RLXSMF)
C
      ENDIF
C
      END
C
C
************************************************************************
*     PROLQ1( DX(1+KOFFX(ILEV-1)), DD(1+KOFFD(ILEV)) )
*     DX2 = Prolongation (interpolation) of DX1 from coarse to fine grid.
************************************************************************
      SUBROUTINE PROLQ1(DX1,DX2)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19)
      INCLUDE 'nnwork.inc'
      DOUBLE PRECISION DX1(*),DX2(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGPAR2/ RLXSMF,DMPSLC,RLXSLC,ISM,ISOLC,NSLC,IREST,
     *                KNX(NNLEV),KNY(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLNPR(NNLEV),KLVERT(NNLEV),KLADJ(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /SET/    ISORT
C-----------------------------------------------------------------------
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
      IF (ISORT.NE.1) THEN
C                  1 = Coarse, 2 = Fine
C            MP011(DUC,DUF,KVERTC,KVERTF,KADJC,KADJF,NVTC,NELC)
        CALL MP011(DX1,DX2,KWORK(L(KLVERT(ILEV-1))),
     *                     KWORK(L(KLVERT(ILEV))),
     *             KWORK(L(KLADJ(ILEV-1))),KWORK(L(KLADJ(ILEV))),
     *             KNVT(ILEV-1),KNEL(ILEV-1))
C
      ELSE
C
        CALL MP011L(DX1,DX2,KNX(ILEV))
C
      ENDIF
C
      END
C
      SUBROUTINE MP011L(DUC,DUF,NXF)
      IMPLICIT NONE
      DOUBLE PRECISION DUC(*),DUF(*),DC1,DC2,DC3,DC4
      INTEGER NXF,NXC,IEL,IC,JC,I,J,IEQC,IEQF
C
      NXC = NXF/2
      IEL = 0
      DO IC=1,NXC
        DO JC=1,NXC
          IEL  = IEL+1
          IEQC = IEL+INT((IEL-1)/NXC)
C
          DC1  = DUC(IEQC)
          DC2  = DUC(IEQC+1)
          DC3  = DUC(IEQC+NXC+2)
          DC4  = DUC(IEQC+NXC+1)
C
          I = 2*(IC-1)+1;
          J = 2*(JC-1)+1;
          IEQF = J+(I-1)*(NXF+1)
C
          DUF(IEQF)    = DC1
          DUF(IEQF+1)  = 0.5D0*(DC1+DC2)
          DUF(IEQF+2)  = DC2
C
          DUF(IEQF+NXF+1) = 0.5D0*(DC1+DC4)
          DUF(IEQF+NXF+2) = 0.25D0*(DC1+DC2+DC3+DC4)
          DUF(IEQF+NXF+3) = 0.5D0*(DC2+DC3)
C
          DUF(IEQF+2*NXF+2) = DC4
          DUF(IEQF+2*NXF+3) = 0.5D0*(DC3+DC4)
          DUF(IEQF+2*NXF+4) = DC3
        ENDDO
      ENDDO
C
      END
C
C
C
C                MP011(DUC,DUF,KVERTC,KVERTF,KADJC,KADJF,NVTC,NELC)
      SUBROUTINE MP011(DU1,DU2,KVERT1,KVERT2,KADJ1,KADJ2,NVT1,NEL1)
C
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (Q2=.5D0,Q4=.25D0)
      DIMENSION DU1(*),DU2(*)
      DIMENSION KVERT1(NNVE,*),KVERT2(NNVE,*)
      DIMENSION KADJ1(NNVE,*),KADJ2(NNVE,*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('MP011 ','01/25/90')
C
      CALL LCP1(DU1,DU2,NVT1)
C
      DO 10 IEL1=1,NEL1
C
      DUH1=DU1(KVERT1(1,IEL1))
      DUH2=DU1(KVERT1(2,IEL1))
      DUH3=DU1(KVERT1(3,IEL1))
      DUH4=DU1(KVERT1(4,IEL1))
C
      IELH1=IEL1
      IELH2=KADJ2(2,IELH1)
      IELH3=KADJ2(2,IELH2)
      IELH4=KADJ2(2,IELH3)
C
      IF (KADJ1(1,IEL1).LT.IEL1)
     *   DU2(KVERT2(2,IEL1))=Q2*(DUH1+DUH2)
C
      IF (KADJ1(2,IEL1).LT.IEL1)
     *   DU2(KVERT2(2,IELH2))=Q2*(DUH2+DUH3)
C
      IF (KADJ1(3,IEL1).LT.IEL1)
     *   DU2(KVERT2(2,IELH3))=Q2*(DUH3+DUH4)
C
      IF (KADJ1(4,IEL1).LT.IEL1)
     *   DU2(KVERT2(2,IELH4))=Q2*(DUH4+DUH1)
C
      DU2(KVERT2(3,IEL1))=Q4*(DUH1+DUH2+DUH3+DUH4)
C
10    CONTINUE
C
      END
C
C
C
************************************************************************
*     RESTQ1( DD(1+KOFFD(ILEV+1)), DB(1+KOFFB(ILEV)) )
*     DD1 = Restriction of DD2 from fine to coarse grid.
************************************************************************
      SUBROUTINE RESTQ1(DD2,DD1)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19)
      INCLUDE 'nnwork.inc'
      DOUBLE PRECISION DD1(*),DD2(*)
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGPAR2/ RLXSMF,DMPSLC,RLXSLC,ISM,ISOLC,NSLC,IREST,
     *                KNX(NNLEV),KNY(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLNPR(NNLEV),KLVERT(NNLEV),KLADJ(NNLEV)
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /SET/    ISORT
C-----------------------------------------------------------------------
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
C
      IF (ISORT.NE.1) THEN
C                  2 = Fine, 1 = Coarse
C            MR011(DXF,DXC,KVERTF,KADJF,NVTC,NELF)
        CALL MR011(DD2,DD1,KWORK(L(KLVERT(ILEV+1))),
     *             KWORK(L(KLADJ(ILEV+1))),
     *             KNVT(ILEV),KNEL(ILEV+1))
C
      ELSE
C
        CALL MR011L(DD2,DD1,KNX(ILEV+1),KNY(ILEV+1),IREST)
C
      ENDIF
C
      END
C
      SUBROUTINE MR011L(DXF,DXC,NXF,NYF,IREST)
      IMPLICIT NONE
      DOUBLE PRECISION Q2,Q4,DXF(*),DXC(*)
      PARAMETER (Q2=0.5D0,Q4=0.25D0)
      INTEGER NXF,NYF,IREST,I,J,IEQC,IEQF
C
C     Pointwise (straight) Injection
      IF (IREST.EQ.1) THEN
C
        IEQC = 0
        DO I=1,NXF+1,2
          DO J=1,NXF+1,2
            IEQC = IEQC+1
            IEQF = J+(I-1)*(NXF+1);
            DXC(IEQC) = 4D0*DXF(IEQF)
          ENDDO
        ENDDO
C
C     Half weighted restriction
      ELSEIF (IREST.EQ.2) THEN
C
        CALL LCL1(DXC,(NXF/2+1)*(NYF/2+1))
C
        IEQC = 0
        DO 1 I=1,NXF+1,2
          DO 2 J=1,NXF+1,2
C
            IEQC = IEQC+1
            IF (I.EQ.1.OR.I.EQ.NYF+1.OR.J.EQ.1.OR.J.EQ.NXF+1) THEN
              GOTO 2
            ENDIF
C
            IEQF = J+(I-1)*(NXF+1)
            DXC(IEQC) = DXC(IEQC)+Q2*DXF(IEQF-1)
            DXC(IEQC) = DXC(IEQC)+2D0*DXF(IEQF)
            DXC(IEQC) = DXC(IEQC)+Q2*DXF(IEQF+1)
C
            IEQF = IEQF-(NXF+1)
            DXC(IEQC) = DXC(IEQC)+Q2*DXF(IEQF)
C
            IEQF = IEQF+2*(NXF+1)
            DXC(IEQC) = DXC(IEQC)+Q2*DXF(IEQF)
2         CONTINUE
1       CONTINUE
C
C       Treat boundary nodes
        IEQC = 0
        DO  I=1,NXF+1,2
          DO J=1,NXF+1,2
            IEQC = IEQC+1
            IF (I.EQ.1.OR.I.EQ.NYF+1.OR.J.EQ.1.OR.J.EQ.NXF+1) THEN
              DXC(IEQC) = 0D0
            ENDIF
          ENDDO
        ENDDO
C
C     Full weighted restriction
      ELSEIF (IREST.EQ.3) THEN
C
        CALL LCL1(DXC,(NXF/2+1)*(NYF/2+1))
C
        IEQC = 0
        DO 5 I=1,NXF+1,2
          DO 6 J=1,NXF+1,2
C
            IEQC = IEQC+1
            IF (I.EQ.1.OR.I.EQ.NYF+1.OR.J.EQ.1.OR.J.EQ.NXF+1) THEN
              GOTO 6
            ENDIF
C
            IEQF = J+(I-1)*(NXF+1)
            DXC(IEQC) = DXC(IEQC)+Q2*DXF(IEQF-1)
            DXC(IEQC) = DXC(IEQC)+   DXF(IEQF)
            DXC(IEQC) = DXC(IEQC)+Q2*DXF(IEQF+1)
C
            IEQF = IEQF-(NXF+1)
            DXC(IEQC) = DXC(IEQC)+Q4*DXF(IEQF-1)
            DXC(IEQC) = DXC(IEQC)+Q2*DXF(IEQF)
            DXC(IEQC) = DXC(IEQC)+Q4*DXF(IEQF+1)
C
            IEQF = IEQF+2*(NXF+1)
            DXC(IEQC) = DXC(IEQC)+Q4*DXF(IEQF-1)
            DXC(IEQC) = DXC(IEQC)+Q2*DXF(IEQF)
            DXC(IEQC) = DXC(IEQC)+Q4*DXF(IEQF+1)
C
6         CONTINUE
5       CONTINUE
C
C       Treat boundary nodes
        IEQC = 0
        DO  I=1,NXF+1,2
          DO J=1,NXF+1,2
            IEQC = IEQC+1
            IF (I.EQ.1.OR.I.EQ.NYF+1.OR.J.EQ.1.OR.J.EQ.NXF+1) THEN
              DXC(IEQC) = 0D0
            ENDIF
          ENDDO
        ENDDO
C
      ENDIF
C
      END
C
C
C
C                MR011(DXF,DXC,KVERTF,KADJF,NVTC,NELF)
      SUBROUTINE MR011(DX2,DX1,KVERT2,KADJ2,NVT1,NEL2)
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      PARAMETER (Q4=.25D0)
      DIMENSION DX1(*),DX2(*)
      DIMENSION KVERT2(NNVE,*),KADJ2(NNVE,*)
C
      CALL LCP1(DX2,DX1,NVT1)
C
      DO 10 IEL2=1,NEL2
C
      I1=KVERT2(1,IEL2)
      I2=KVERT2(2,IEL2)
      I3=KVERT2(3,IEL2)
      I4=KVERT2(4,IEL2)
C
      DX1(I1)=DX1(I1)+Q4*(DX2(I2)+DX2(I3)+DX2(I4))
      IF (KADJ2(1,IEL2).EQ.0) DX1(I1)=DX1(I1)+Q4*DX2(I2)
      IF (KADJ2(4,IEL2).EQ.0) DX1(I1)=DX1(I1)+Q4*DX2(I4)
C
10    CONTINUE
C
      END
C
C
C
C     4-color SOR
      SUBROUTINE IC217_4C(DA,KCOL,KLD,DX,DB,NEQ,NIT,OMEGA,NX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DB(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('IC2174','10/27/89')
C
      NY = NEQ/NX  ! Number of sections/horizontal rows with NX nodes
                   ! equal to the number of rows of vertices in grid
C
      DO 1 ITE=1,NIT
C
C     Loop over odd grid rows
      DO 20 IROW=3,NY-1,2
C
C       Loop over odd grid columns
        DO 25 ICOL=3,NX-1,2
C
          IEQ = ICOL+NX*(IROW-1)
C
          AUX = DB(IEQ)
          DO ICOLA=KLD(IEQ)+1,KLD(IEQ+1)-1
            AUX = AUX - DA(ICOLA)*DX(KCOL(ICOLA))
          ENDDO
          AUX     = OMEGA*(AUX/DA(KLD(IEQ)) - DX(IEQ))+DX(IEQ)
          DX(IEQ) = AUX
25      CONTINUE
C
C       Loop over even grid columns
        DO 28 ICOL=2,NX-1,2
C
          IEQ = ICOL+NX*(IROW-1)
C
          AUX = DB(IEQ)
          DO ICOLA=KLD(IEQ)+1,KLD(IEQ+1)-1
            AUX = AUX - DA(ICOLA)*DX(KCOL(ICOLA))
          ENDDO
          AUX     = OMEGA*(AUX/DA(KLD(IEQ)) - DX(IEQ))+DX(IEQ)
          DX(IEQ) = AUX
28      CONTINUE
C
20    CONTINUE
C
C     Loop over even grid rows
      DO 30 IROW=2,NY-1,2
C
C       Loop over odd grid columns
        DO 35 ICOL=3,NX-1,2
C
          IEQ = ICOL+NX*(IROW-1)
C
          AUX = DB(IEQ)
          DO ICOLA=KLD(IEQ)+1,KLD(IEQ+1)-1
            AUX = AUX - DA(ICOLA)*DX(KCOL(ICOLA))
          ENDDO
          AUX     = OMEGA*(AUX/DA(KLD(IEQ)) - DX(IEQ))+DX(IEQ)
          DX(IEQ) = AUX
35      CONTINUE
C
C       Loop over even grid columns
        DO 38 ICOL=2,NX-1,2
C
          IEQ = ICOL+NX*(IROW-1)
C
          AUX = DB(IEQ)
          DO ICOLA=KLD(IEQ)+1,KLD(IEQ+1)-1
            AUX = AUX - DA(ICOLA)*DX(KCOL(ICOLA))
          ENDDO
          AUX     = OMEGA*(AUX/DA(KLD(IEQ)) - DX(IEQ))+DX(IEQ)
          DX(IEQ) = AUX
38      CONTINUE
C
30    CONTINUE
C
1     CONTINUE
C
      END
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* XMSB2                                                                *
*                                                                      *
* Purpose  Generate sequence of meshes for multigrid applications      *
*          by successive calls of XSB0X  (two-level ordering)          *
*                                                                      *
* Subroutines/functions called  XSB0X, WERR, ZCPY                      *
*                                                                      *
* Version from  05/11/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    desired number of l                                  *
* IMID     I*4                                                         *
* IADJ     I*4                                                         *
* IVEL     I*4                                                         *
* IDISP    I*4                                                         *
* IBDP     I*4     See parameter list of XSB0X                         *
* S2DI     SUBR                                                        *
* S2DB     SUBR                                                        *
* PARX     FNCT                                                        *
* PARY     FNCT                                                        *
* TMAX     FNCT                                                        *
* Coarse grid on /TRIAD and /TRIAA/                                    *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DCORVG   R*8   Cartesian coordinates of vertices                     *                                                                      *                                                                      *
************************************************************************
C
      SUBROUTINE XMSB2(IMID,IADJ,IVEL,IDISP,IBDP,
     *                 S2DI,S2DB,PARX,PARY,TMAX)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19)
      INCLUDE 'nnwork.inc'
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGTRA/  KLCVG(NNLEV),KLCMG(NNLEV),KLVERT(NNLEV),
     *                KLMID(NNLEV),KLADJ(NNLEV),KLVEL(NNLEV),
     *                KLMEL(NNLEV),KLNPR(NNLEV),KLMM(NNLEV),
     *                KLVBD(NNLEV),KLEBD(NNLEV),KLBCT(NNLEV),
     *                KLVBDP(NNLEV),KLMBDP(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C
      EXTERNAL S2DI,S2DB,PARX,PARY,TMAX
      SAVE
C
      SUB='XMSB2 '
      IF (ICHECK.GE.997) CALL OTRC('XMSB2 ','11/05/91')
C
      NLEV=NLMAX
C
      IF (NLEV.GT.NNLEV) THEN
       CALL WERR(-180,'XMSB2 ')
       GOTO 99999
      ENDIF
C
C
      IADJ0=1
      IBDP0=2
C
C
      DO 10 ILEV=1,NLEV
      IF (ILEV.EQ.NLEV) THEN
       IADJ0=IADJ
       IBDP0=IBDP
      ENDIF
      IF (ILEV.EQ.1) THEN
       CALL XSB0X(0,IMID,IADJ0,IVEL,IDISP,IBDP0,
     *            S2DI,S2DB,PARX,PARY,TMAX)
      ELSE
       CALL XSB0X(1,IMID,IADJ0,IVEL,IDISP,IBDP0,
     *            S2DI,S2DB,PARX,PARY,TMAX)
      ENDIF
C
************************************************************************
C
C      GREPS=0.20D0
C      IF (ILEV.EQ.NLEV) CALL GRDIST(DWORK(L(LCORVG)),KWORK(L(LNPR)),
C     *                              NVT,GREPS)
C
************************************************************************
C
      IF (ILEV.GE.2)
     * CALL CHCOOR(DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LMID)),
     *             KWORK(L(LADJ)),KWORK(L(LNPR)),NEL,NVT)
************************************************************************
C
      IF (ILEV.GE.NLMIN) THEN
C
C *** Save dimensions for all levels between NLMIN and NLMAX
C
      KNEL(ILEV) =NEL
      KNVT(ILEV) =NVT
      KNMT(ILEV) =NMT
      KNVEL(ILEV)=NVEL
      KNVBD(ILEV)=NVBD
C
C
C *** Save arrays for all levels between NLMIN and NLMAX
C
C
      IF (ILEV.LT.NLEV) THEN
C
C       CALL ZCPY(LCORVG,'DCORVG',KLCVG(ILEV),'DCVG0 ')
C       IF (IER.NE.0) GOTO 99999
       IF (LCORMG.NE.0) THEN
        CALL ZCPY(LCORMG,'DCORMG',KLCMG(ILEV),'DCMG0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LVERT,'KVERT ',KLVERT(ILEV),'KVERT0')
       IF (IER.NE.0) GOTO 99999
       IF (LMID.NE.0) THEN
        CALL ZCPY(LMID,'KMID  ',KLMID(ILEV),'KMID0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IADJ.EQ.1) THEN
        CALL ZCPY(LADJ,'KADJ  ',KLADJ(ILEV),'KADJ0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LVEL.NE.0) THEN
        CALL ZCPY(LVEL,'KVEL  ',KLVEL(ILEV),'KVEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (LMEL.NE.0) THEN
        CALL ZCPY(LMEL,'KMEL  ',KLMEL(ILEV),'KMEL0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LNPR,'KNPR  ',KLNPR(ILEV),'KNPR0 ')
       IF (IER.NE.0) GOTO 99999
       CALL ZCPY(LMM,'KMM   ',KLMM(ILEV),'KMM0  ')
       IF (IER.NE.0) GOTO 99999
       IF (IBDP.GE.0) THEN
        CALL ZCPY(LVBD,'KVBD  ',KLVBD(ILEV),'KVBD0 ')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IBDP.GE.1) THEN
       CALL ZCPY(LEBD,'KEBD  ',KLEBD(ILEV),'KEBD0 ')
       IF (IER.NE.0) GOTO 99999
       ENDIF
       CALL ZCPY(LBCT,'KLBCT ',KLBCT(ILEV),'KBCT0 ')
       IF (IER.NE.0) GOTO 99999
       IF (IBDP.GE.2) THEN
        CALL ZCPY(LVBDP,'DVBDP ',KLVBDP(ILEV),'DVBDP0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
       IF (IBDP.GE.2.AND.LMBDP.NE.0) THEN
        CALL ZCPY(LMBDP,'DMBDP ',KLMBDP(ILEV),'DMBDP0')
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
      ELSE
C
       KLCVG(ILEV) =LCORVG
       KLCMG(ILEV) =LCORMG
       KLVERT(ILEV)=LVERT
       KLMID(ILEV) =LMID
       KLADJ(ILEV) =LADJ
       KLVEL(ILEV) =LVEL
       KLMEL(ILEV) =LMEL
       KLNPR(ILEV) =LNPR
       KLMM(ILEV)  =LMM
       KLVBD(ILEV) =LVBD
       KLEBD(ILEV) =LEBD
       KLBCT(ILEV) =LBCT
       KLVBDP(ILEV)=LVBDP
       KLMBDP(ILEV)=LMBDP
C
      ENDIF
C
      DO 20 ILEVH=1,NLEV-1
      KLCVG(ILEVH)=LCORVG
20    CONTINUE
C
      ENDIF
C
10    CONTINUE
C
C
99999 END
C
C
C
C************************************************************************
      SUBROUTINE CHCOOR(DCORVG,KVERT,KMID,KADJ,KNPR,NEL,NVT)
      IMPLICIT REAL*8 (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION DCORVG(2,*),KVERT(4,*),KMID(4,*),KADJ(4,*),KNPR(*)
C
C
C
      DO 100 IEL=1,NEL/4
C
      IEL1=IEL
      IEL2=KADJ(2,IEL1)
      IEL3=KADJ(2,IEL2)
      IEL4=KADJ(2,IEL3)
C
      IADJ1=KADJ(1,IEL1)
      IADJ2=KADJ(1,IEL2)
      IADJ3=KADJ(1,IEL3)
      IADJ4=KADJ(1,IEL4)
C
      IVT1=KVERT(2,IEL1)
      IVT2=KVERT(2,IEL2)
      IVT3=KVERT(2,IEL3)
      IVT4=KVERT(2,IEL4)
C
      NADJ0=0
      IF (IADJ1.EQ.0) NADJ0=NADJ0+1
      IF (IADJ2.EQ.0) NADJ0=NADJ0+1
      IF (IADJ3.EQ.0) NADJ0=NADJ0+1
      IF (IADJ4.EQ.0) NADJ0=NADJ0+1
C
C
      IF (NADJ0.EQ.0) GOTO 100
C
C
      PX1=DCORVG(1,IVT1)
      PX2=DCORVG(1,IVT2)
      PX3=DCORVG(1,IVT3)
      PX4=DCORVG(1,IVT4)
C
      PY1=DCORVG(2,IVT1)
      PY2=DCORVG(2,IVT2)
      PY3=DCORVG(2,IVT3)
      PY4=DCORVG(2,IVT4)
C
      IVTM=KVERT(3,IEL)
      PXM=DCORVG(1,IVTM)
      PYM=DCORVG(2,IVTM)
C
C
      IF (NADJ0.EQ.1) THEN
       IF ((IADJ1.EQ.0).OR.(IADJ3.EQ.0)) THEN
        PX=0.5D0*(PX1+PX3)
        PY=0.5D0*(PY1+PY3)
       ELSE
        PX=0.5D0*(PX2+PX4)
        PY=0.5D0*(PY2+PY4)
       ENDIF
       DCORVG(1,IVTM)=PX
       DCORVG(2,IVTM)=PY
      ENDIF
C
C
      IF (NADJ0.GT.1) THEN
       PX=0.25D0*(PX1+PX2+PX3+PX4)
       PY=0.25D0*(PY1+PY2+PY3+PY4)
       DCORVG(1,IVTM)=PX
       DCORVG(2,IVTM)=PY
      ENDIF
C
100   CONTINUE
C
C
C
      END
