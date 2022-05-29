************************************************************************
      PROGRAM  FEATFEM
************************************************************************
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19,NNAB=21)
      INCLUDE 'nnwork.inc'
      CHARACTER SUB*6,FMT*15,CPARAM*1
C
C *** IO file names
C     CDATA  - Input parameter file   [hpcfem.dat]
C     CFILE  - Data output log file   [*.out]
      CHARACTER CDATA*60,CFILE*60,CMSG*45,CLINE*79
C
C-----------------------------------------------------------------------
C     C O M M O N S
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM(120)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      COMMON /GEOM/   DLX,DLY
C
C *** NBLOCA=1 means that the mass matrix consists of one block
      PARAMETER (NBLOCA=1)
C
C *** Matrix and load vector names (for messages only)
      CHARACTER ARRA*6,ARRF*6
C
C *** KABA  - structure of the mass matrix
      DIMENSION KABA(2,NNAB,NBLOCA)
C *** KABAN - number of terms in the integrand for the mass matrix
      DIMENSION KABAN(NBLOCA)
C *** KB    - structure of linear form
      DIMENSION KB(NNAB,1)
C
C *** Boolean for declaring constant coefficents
      DIMENSION BCONA(NBLOCA)
C *** Boolean not to convert to single precision
      DIMENSION BSNGLA(NBLOCA)
C-----------------------------------------------------------------------
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGPAR2/ RLXSMF,DMPSLC,RLXSLC,ISM,ISOLC,NSLC,IREST,
     *                KNX(NNLEV),KNY(NNLEV)
      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),KLLDA(NNLEV),
     *                KLNPR(NNLEV),KLVERT(NNLEV),KLADJ(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      COMMON /SET/    ISORT
      INTEGER KNEQ(NNLEV),KNA(NNLEV)
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C-----------------------------------------------------------------------
      EXTERNAL E011,COEF1,UE
      EXTERNAL YAXU,PROLQ1,RESTQ1,YSMU,YEXU,YDBCU
      EXTERNAL PARX,PARY,TMAX,S2DI0,S2DB0
C
      DATA ARRA/'DA    '/
      DATA ARRF/'DF    '/
      DATA BCONA/.TRUE./
      DATA BSNGL/.FALSE./
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CALL ZTIME(TTTSUM)
      SUB = 'MAIN  '
C
      CALL ZTIME(TTT00)
C
C     Initialization of BLANK COMMON and of I/O devices
      CALL ZINIT(NNWORK,'feat.msg','data/featfem.err',
     *                             'data/featfem.prt',
     *                             'data/featfem.sys',
     *                             'data/featfem.trc')
C
      CALL ZTIME(TTT1)
      TTLC = TTLC+TTT1-TTT0
C
C=======================================================================
C     Read input data file
C=======================================================================
C *** Read unit numbers and file names on /FILES/:
C-----------------------------------------------------------------------
C
      CLINE = '----------------------------------------'//
     *        '---------------------------------------'
C
      MSHOW = 0
C
      MDATA = 79
      CDATA = 'data/featfem.dat'
C
      CALL OF0(MDATA,CDATA,1)
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
      CALL FTOUT(MSHOW,MFILE,MTERM,'         INPUT DATA (FEATFEM)')
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
      READ(MDATA,*) CFILE
      MFILE = 63
      CMSG = 'MFILE                             : MFILE  = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,MFILE)
      CMSG = 'CFILE                             : CFILE  = '
      CALL FTOUTC(MSHOW,MFILE,MTERM,CMSG,CFILE)
C
      CALL OF0(MFILE,CFILE,1)
C
C-----------------------------------------------------------------------
C *** Read values for /OUTPUT/
C-----------------------------------------------------------------------
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*) M
      M = ABS(M)
C
      READ(MDATA,*) MT
      MT = ABS(MT)
C
      READ(MDATA,*) ICHECK
      ICHECK = ABS(ICHECK)
C
      READ(MDATA,*) MSHOW
      MSHOW = ABS(MSHOW)
C
      READ(MDATA,*) IOUT
C
      CMSG = 'message level for file output     : M      = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,M)
      CMSG = 'message level for terminal output : MT     = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,MT)
      CMSG = 'error tracing level               : ICHECK = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,ICHECK)
      CMSG = 'message output level              : MSHOW  = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,MSHOW)
      CMSG = 'write solution to file            : IOUT   = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,IOUT)
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
C-----------------------------------------------------------------------
C *** Read values for grid
C-----------------------------------------------------------------------
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*)  NX
      READ(MDATA,*)  NY
      IF (NY.LE.0) NY = NX
      READ(MDATA,*)  DLX
      READ(MDATA,*)  DLY
      READ(MDATA,*)  ISORT
C
      CMSG = 'number of horizontal cells        : NX     = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,NX)
      CMSG = 'number of verticalc cells         : NY     = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,NY)
      CMSG = 'horizontal domain length          : DLX    = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,DLX)
      CMSG = 'vertical domain length            : DLY    = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,DLY)
      CMSG = 'grid sorting/ordering             : ISORT  = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,ISORT)
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
C-----------------------------------------------------------------------
C *** Read values for assembly
C-----------------------------------------------------------------------
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*)  IASM
      READ(MDATA,*)  ICUB
C
      CMSG = 'assembly routine                  : IASM   = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,IASM)
      CMSG = 'cubature rule                     : ICUB   = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,ICUB)
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
C-----------------------------------------------------------------------
C *** Read values for linear solver
C-----------------------------------------------------------------------
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*)  ISOL
      READ(MDATA,*)  NIT
      READ(MDATA,*)  OMEGA
      READ(MDATA,*)  EPS
C
      CMSG = 'linear solver method              : ISOL   = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,ISOL)
      ISOLMG = 6
      CMSG = 'maximum number of iterations      : NIT    = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,NIT)
      CMSG = 'relaxation parameter              : OMEGA  = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,OMEGA)
      CMSG = 'tolerance/stopping criteria       : EPS    = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,EPS)
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
C-----------------------------------------------------------------------
C *** Read values for multigrid
C-----------------------------------------------------------------------
C
      READ(MDATA,*)
      READ(MDATA,*)
      READ(MDATA,*)
C
      READ(MDATA,*)  NLEV
      IF (ISOL.NE.ISOLMG) NLEV = 1
      READ(MDATA,*)  ICYCLE
      READ(MDATA,*)  IRELMG
      READ(MDATA,*)  ISM
      READ(MDATA,*)  NPRESM
      READ(MDATA,*)  NPOSSM
      READ(MDATA,*)  RLXSMF
      READ(MDATA,*)  NSMFAC
      READ(MDATA,*)  IREST
      READ(MDATA,*)  ISOLC
      READ(MDATA,*)  DMPSLC
      READ(MDATA,*)  NSLC
      READ(MDATA,*)  RLXSLC
      READ(MDATA,*)  IMTIME
C
      CMSG = 'number of grid levels             : NLEV   = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,NLEV)
      CMSG = 'mg-cycle: 0=F, 1=V, 2=W)          : ICYCLE = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,ICYCLE)
      CMSG = 'monitor relative/absolute defect  : IRELMG = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,IRELMG)
      CMSG = 'smoother: 1=JAC, 2=SOR, 3=SSOR    : ISM    = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,ISM)
      CMSG = 'number of presmoothing steps      : NPRESM = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,NPRESM)
      CMSG = 'number of postsmoothing steps     : NPOSSM = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,NPOSSM)
      CMSG = 'relaxation for the mg solver      : RLXSMF = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,RLXSMF)
      CMSG = 'factor for sm. on coarser levels  : NSMFAC = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,NSMFAC)
      CMSG = 'restricton operator               : IREST  = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,IREST)
      CMSG = 'coarse grid solver 1=JAC, 2=SOR   : ISOLC  = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,ISOLC)
      CMSG = 'damping of residuals for mg-it.   : DMPSLC = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,DMPSLC)
      CMSG = 'maximum # of c-grid solver it.    : NSLC   = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,NSLC)
      CMSG = 'relaxation for the c-grid solver  : RLXSLC = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,RLXSLC)
      CMSG = 'check mg timings                  : IMTIME = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,IMTIME)
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
C=======================================================================
C     Grid generation
C=======================================================================
C
      NVE = 4
C
      IF (ISORT.EQ.1) THEN
        DO I=NLEV,1,-1
          NLMIN = I
          IF (NX/(2**(NLEV-I)).EQ.2.OR.NY/(2**(NLEV-I)).EQ.2) GOTO 2
        ENDDO
2       CONTINUE
        NLMAX  = NLEV-(NLMIN-1)
        NLMIN  = 1
        NLEV   = NLMAX
      ELSE
        DO I=1,1000
          IF (NX/(2**I).LE.3.OR.NY/(2**I).LE.3) GOTO 3
        ENDDO
3       CONTINUE
        NLMAX = I+1
        NLMIN = MAX(1,NLMAX-8)
        NLEV  = NLMAX-NLMIN+1
        NX0   = NX/(2**(NLEV-1))
        NY0   = NY/(2**(NLEV-1))
c$$$        WRITE (*,*) NLMIN,NLMAX,NLEV,NX0,NY0
      ENDIF
C
      TTGRID = 0D0
      TTPTR  = 0D0
      TTASM  = 0D0
      TTBDR  = 0D0
C
C
      IF (ISORT.NE.1) THEN
        CALL XUSQGRD(NX0,NY0,DLX,DLY,LCORVG,NVT,NMT,LVERT,LNPR,NEL,1)
        CALL XS2A
      ENDIF
      NBCT = 4
      LMID = 0
C
C *** Loop over all multigrid levels
      DO 100 ILEV=NLMIN,NLMAX   !!!
C
      KPRSM(ILEV) = NPRESM*NSMFAC**(NLMAX-ILEV)
      KPOSM(ILEV) = NPOSSM*NSMFAC**(NLMAX-ILEV)
C
      CALL ZTIME(TTT0)
      IF (ISORT.EQ.1) THEN
        CALL XUSQGRD(NX/(2**(NLEV-ILEV)),NY/(2**(NLEV-ILEV)),
     *               DLX,DLY,LCORVG,NVT,NMT,LVERT,LNPR,NEL,0)
        LADJ = 0
      ELSEIF (ILEV.GT.NLMIN) THEN
        CALL ZNEW(4*4*NEL,3,LV1,'KVERT ')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY(LVERT,'KVOLD ',LV1,'KVERT ')
        LVERT=LV1
        CALL ZNEW(4*4*NEL,3,LA1,'KADJ  ')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY(LADJ,'KAOLD ',LA1,'KADJ  ')
        LADJ=LA1
        CALL ZNEW(2*4*4*NEL,1,LCV1,'DCORVG')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY(LCORVG,'DCVOLD',LCV1,'DCORVG')
        CALL ZDISP(0,LCORVG,'DCORVG')
        LCORVG=LCV1
        CALL ZNEW(4*4*NEL,3,LNPR1,'KNPR  ')
        IF (IER.NE.0) GOTO 99999
        CALL ZCPY(LNPR,'KNPOLD',LNPR1,'KNPR  ')
        LNPR=LNPR1
C
        IDISP = 1
        CALL XSB0(1,IDISP)
C
      ENDIF
C
      KNEL(ILEV)   = NEL
      KNVT(ILEV)   = NVT
      NMT = 0 !??? Bug or assembly crashes
      KLVERT(ILEV) = LVERT
      KLADJ(ILEV)  = LADJ
      KLNPR(ILEV)  = LNPR
      IF (ISORT.EQ.1) THEN
        KNX(ILEV)    = NX/(2**(NLEV-ILEV))
        KNY(ILEV)    = NY/(2**(NLEV-ILEV))
      ELSE
        KNX(ILEV)    = NX/(2**(NLMAX-ILEV))
        KNY(ILEV)    = NY/(2**(NLMAX-ILEV))
      ENDIF
C
      CALL ZTIME(TTT1)
      TTGRID = TTGRID+TTT1-TTT0
C
C     Store pointer to last grid element on DWORK and
C     the total number of stored elements on DWORK
      IWORKG = IWORK
      IWMAXG = IWMAX
C
C=======================================================================
C     Allocate space for arrays
C=======================================================================
C
      IF (ILEV.EQ.NLMAX) THEN
        CALL ZNEW(NVT,1,LU,'DU    ')
        IF (IER.NE.0) GOTO 99998
      ENDIF
C
C======================================================================
C     Matrix and adress pointer initializations
C======================================================================
C
C     ==================================================================
C *** Calculate the pointer vectors for the matrix
C *** Storage technique 7 - symmetry of the matrix is neglected
C *** which facilitates the implementation boundary conditions
C     ==================================================================
C     Outputs:
C     LCOLA  - Pointers to start of new columns (adress)
C     LLDA   - Pointers to start of new rows (adress)
C     NA     - Number of entries in matrix
C     NEQ    - Number of equations
C
C     Inputs:
C     E011   - Element subroutine
C     ISYMM  - Matrix symmetry switch (0 => no symmetry assumed)
C
C     Subcall to AP7 takes the adresses L(LVERT) and L(LMID)
C     in COMMON /TRIAA/
C     ------------------------------------------------------------------
C
      CALL ZTIME(TTT0)
      LCOLA = 0
      LLDA  = 0
      NA    = 0
      NEQ   = 0
      ISYMM = 0
C
      CALL XAP7(LCOLA,LLDA,NA,NEQ,E011,ISYMM)
      IF (IER.NE.0) GOTO 99998
      CALL ZTIME(TTT1)
      TTPTR = TTPTR+TTT1-TTT0
C
      KLCOLA(ILEV) = LCOLA
      KLLDA(ILEV)  = LLDA
C
C=======================================================================
C     Feat assembly of a matrix
C=======================================================================
C
      CALL ZTIME(TTT0)
      IF (IASM.EQ.1) THEN
C
C     ==================================================================
C *** Assembly of a diffusion matrix
C     ==================================================================
C     XAB07 Assemble bilinear matrix, quadrilaterals, area integral,
C           storage technique 7 (Initialization also done)
C
C     LA     - Number of array (LNR)
C     LCOLA  - Column pointer
C     LLDA   - Pointer to start of new rows
C     NA     - Number of non-zero elements per matrix
C     NEQ    - Number of equations (rows) per matrix
C     NBLOCA - 1 => One matrix block stored on LA
C     ICRLA  - 1 => Clear all matrices
C     E011   - Element subroutine
C     COEFFA - Coefficient subroutine
C     BCONA  - .TRUE. => the block has constant coefficients
C     KABA   - Bilinear form structure
C     KABAN  - Pairs of multi-indices occuring in the bilinear forms
C              specified separately for each matrix
C     ABS(ICUBA) - The number of quadrature formula
C     ISYMMA - 0 => No symmetry assumed
C     ARRA   - Name of array LA, used for messages only
C     BSNGLA - .FALSE. => don't change to single precision
C     ------------------------------------------------------------------
C
      LA          = 0
      ICLR        = 1
      KABA(1,1,1) = 2
      KABA(2,1,1) = 2
      KABA(1,2,1) = 3
      KABA(2,2,1) = 3
      KABAN(1)    = 2
      ISYMM       = 0
C
      CALL XAB07(LA,LCOLA,LLDA,NA,NEQ,NBLOCA,ICLR,E011,
     *           COEF1,BCONA,KABA,KABAN,ICUB,ISYMM,
     *           ARRA,BSNGL)
      IF (IER.NE.0) GOTO 99998
C
      KNA(ILEV) = NA
C
C=======================================================================
C     Custom assembly of a matrix
C=======================================================================
C
      ELSEIF (IASM.EQ.2) THEN
C
        CALL ZNEW(NA,1,LA,'DA    ')
        IF (IER.NE.0) GOTO 99998
        CALL ZNEW(16,1,LENTRY,'DENTRY')
        IF (IER.NE.0) GOTO 99998
C
        CALL ASMMAT(DWORK(L(LA)), KWORK(L(LCOLA)),KWORK(L(LLDA)),NA,NEQ,
     *              KWORK(L(LVERT)),KWORK(L(LMID)),DWORK(L(LCORVG)),
     *              E011,COEF1,BCONA,ICUB,DWORK(L(LENTRY)))
C
      ENDIF
C
      KLA(ILEV)  = LA
      KNEQ(ILEV) = NEQ
C
C     Release grid storage arrays
      IF (ILEV.LT.NLMAX.AND.ISORT.EQ.1) THEN
C
        CALL ZDISP(0,LVERT, 'KVERT ')
        IF (IER.NE.0) GOTO 99998
C
        CALL ZDISP(0,LCORVG,'DCORVG')
        IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
      CALL ZTIME(TTT1)
      TTASM = TTASM+TTT1-TTT0
C
C=======================================================================
C     Calculate right hand side (load) vector
C=======================================================================
C
      IF (ILEV.EQ.NLMAX) THEN
C
        CALL ZTIME(TTT0)
        NBLOCF  = 1
        KB(1,1) = 1
        KBN     = 1
        BCONF   = .TRUE.
        LB      = 0
C
        CALL XVB0(LB,NEQ,NBLOCF,ICLR,E011,
     *            COEF1,BCONF,KB,KBN,ICUB,ARRF,BSNGL)
        CALL ZTIME(TTT1)
        TTRHS = TTT1-TTT0
        IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
      IF (ILEV.EQ.NLMAX.AND.IOUT.EQ.0) THEN
C
        CALL ZDISP(0,LCORVG,'DCORVG')
        IF (IER.NE.0) GOTO 99998
C
      ENDIF
C
C=======================================================================
C     Set Dirichlet boundary conditions
C=======================================================================
C
      CALL ZTIME(TTT0)
      IF (ILEV.EQ.NLMAX) THEN
        CALL BDRSET(DWORK(L(LU)),NVT,KWORK(L(LNPR)),
     *              DWORK(L(LCORVG)),UE)
        CALL BDRSET(DWORK(L(LB)),NVT,KWORK(L(LNPR)),
     *              DWORK(L(LCORVG)),UE)
      ENDIF
      CALL BDRYA(DWORK(L(LA)),KWORK(L(LCOLA)),KWORK(L(LLDA)),
     *           NVT,KWORK(L(LNPR)))
      CALL ZTIME(TTT1)
      TTBDR = TTBDR+TTT1-TTT0
      IF (IER.NE.0) GOTO 99998
C
C
100   CONTINUE  !!!
C *** END Loop over all multigrid levels
C
C=======================================================================
C     Display grid statistics
C=======================================================================
C
      CMSG = 'time for grid generation          : TTGRID = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,TTGRID)
C
      DO 11 I=NLMIN,NLMAX
C
C       Display/write subdivision information
        IF (MSHOW.GE.2)
     *    WRITE (MTERM,'(A,I2,I12,I12,I12,I12)')
     *      ' ILEV,NVT,NEL,NEQ,NA: ',
     *      I,KNVT(I),KNEL(I),KNEQ(I),KNA(I)
        IF (MSHOW.GE.0)
     *    WRITE (MFILE,'(A,I2,I12,I12,I12,I12)')
     *      ' ILEV,NVT,NEL,NEQ,NA: ',
     *      I,KNVT(I),KNEL(I),KNEQ(I),KNA(I)
C
11    CONTINUE
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
C=======================================================================
C     Display timings
C=======================================================================
C
      CMSG = 'time for FEAT2D ptr allocation    : TTPTR  = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,TTPTR)
C
      CMSG = 'time for matrix assembly          : TTASM  = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,TTASM)
C
      CMSG = 'time for FEAT rhs assembly        : TTRHS  = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,TTRHS)
C
      CMSG = 'time for boundary conditions      : TTBDR  = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,TTBDR)
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
C=======================================================================
C     Sparse MV
C=======================================================================
C
      CALL ZNEW(NEQ,-1,LH,'DH    ')
      IF (IER.NE.0) GOTO 99998
      CALL ZTIME(TTT0)
      DO I=1,100
      CALL LAX17(DWORK(L(KLA(NLMAX))),KWORK(L(KLCOLA(NLMAX))),
     *           KWORK(L(KLLDA(NLMAX))),NEQ,DWORK(L(LU)),
     *           DWORK(L(LH)),1D0,1D0)
      ENDDO
      CALL ZTIME(TTT1)
      TTSPMV = (TTT1-TTT0)/100D0
      CALL ZDISP(0,LH,'DH    ')
      IF (IER.NE.0) GOTO 99998

      CMSG = 'time for sparse MV                : TTSPMV  = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,TTSPMV)
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
C=======================================================================
C     Compute solution
C=======================================================================
C
      CALL ZTIME(TTT0)
      LD = 0
C
C     Jacobi method.
      IF (ISOL.EQ.1) THEN
C
        CALL ZNEW(NEQ,-1,LD,'DD    ')
        IF (IER.NE.0) GOTO 99998
C
        CALL IA017(DWORK(L(LA)),KWORK(L(LCOLA)),KWORK(L(LLDA)),
     *             DWORK(L(LU)),DWORK(L(LB)),DWORK(L(LD)),
     *             NEQ,NIT,ITE,EPS,OMEGA)
C
C     Gauss-Seidel method.
      ELSEIF (ISOL.EQ.2) THEN
C
        CALL IB017(DWORK(L(LA)),KWORK(L(LCOLA)),KWORK(L(LLDA)),
     *             DWORK(L(LU)),DWORK(L(LB)),
     *             NEQ,NIT,ITE,EPS)
C
C     SOR method.
      ELSEIF (ISOL.EQ.3) THEN
C
        CALL IC017(DWORK(L(LA)),KWORK(L(LCOLA)),KWORK(L(LLDA)),
     *             DWORK(L(LU)),DWORK(L(LB)),
     *             NEQ,NIT,ITE,EPS,OMEGA)
C
C     4-color SOR
      ELSEIF (ISOL.EQ.4) THEN
C
        CALL IC017_4C(DWORK(L(LA)),KWORK(L(LCOLA)),KWORK(L(LLDA)),
     *             DWORK(L(LU)),DWORK(L(LB)),
     *             NEQ,NIT,ITE,EPS,OMEGA,NX+1)
C
C     Preconditioned Conjugate-Gradient (PCG) algorithm.
      ELSEIF (ISOL.EQ.5) THEN
C
        CALL ZNEW(4*NEQ,-1,LD,'DD    ')
        IF (IER.NE.0) GOTO 99998
C
        CALL IE017(DWORK(L(LA)),KWORK(L(LCOLA)),KWORK(L(LLDA)),
     *             DWORK(L(LU)),DWORK(L(LB)),NEQ,NIT,ITE,EPS,OMEGA,
     *             DWORK(L(LD)),DWORK(L(LD)+NVT),
     *             DWORK(L(LD)+2*NVT),DWORK(L(LD)+3*NVT))
C
C     Multigrid
      ELSEIF (ISOL.EQ.ISOLMG) THEN
C
        CALL ZTIME(TTT22)
        IDISP = 1
        ITE = 0   ! In: Min numbe of iterations
        CALL XM010L(LU,LB,KNEQ,NIT,ITE,EPS,
     *              YAXU,PROLQ1,RESTQ1,YSMU,
     *              YSMU,YEXU,YDBCU,IDISP,IRELMG,TTSOL)
        CALL ZTIME(TTT22)
C
      ENDIF
C
      IF (IER.GT.0) WRITE (*,*) CPARAM
      IF (ISOL.NE.ISOLMG) THEN
        CALL ZTIME(TTT1)
        TTSOL = TTT1-TTT0
      ENDIF
C
C     Calculate final defect norm
      NEQ = KNEQ(NLMAX)
      IF (LD.EQ.0) CALL ZNEW(NEQ,-1,LD,'DD    ')
      CALL LCP1(DWORK(L(LB)),DWORK(L(LD)),NEQ)
      CALL LAX17(DWORK(L(KLA(NLMAX))),KWORK(L(KLCOLA(NLMAX))),
     *           KWORK(L(KLLDA(NLMAX))),NEQ,DWORK(L(LU)),
     *           DWORK(L(LD)),-1D0,1D0)
      CALL LL21(DWORK(L(LD)),NEQ,DNM)
C
      CALL ZDISP(0,LD,'DD    ')
      IF (IER.NE.0) GOTO 99998
C
      CMSG = 'L2 norm of final defect           :        = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,DNM)
      CMSG = 'number of iterations              :        = '
      CALL FTOUTI(MSHOW,MFILE,MTERM,CMSG,ITE)
      CMSG = 'time for solution                 :        = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,TTSOL)
!      READ(CPARAM,'(I15,D25.16)') ITE,DPREC
!      CMSG = 'final precision reached           :        = '
!      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,DPREC)
C
C=======================================================================
C     Write mg timings
C=======================================================================
C
      IF (ISOL.EQ.ISOLMG.AND.IMTIME.EQ.1) THEN
C
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
      CMSG = 'mg time                           : TTMG   = '
      CALL FTOUTD(MSHOW,MFILE,MTERM,CMSG,TTMG)
      CMSG = 'mg time for smooting              : TTS    = '
      CALL FTOUTP(MSHOW,MFILE,MTERM,CMSG,TTS,1.D2*TTS/TTMG)
      CMSG = 'mg time for defect calculation    : TTD    = '
      CALL FTOUTP(MSHOW,MFILE,MTERM,CMSG,TTD,1.D2*TTD/TTMG)
      CMSG = 'mg time for prolongation          : TTP    = '
      CALL FTOUTP(MSHOW,MFILE,MTERM,CMSG,TTP,1.D2*TTP/TTMG)
      CMSG = 'mg time for restriction           : TTR    = '
      CALL FTOUTP(MSHOW,MFILE,MTERM,CMSG,TTR,1.D2*TTR/TTMG)
      CMSG = 'mg time for solver                : TTE    = '
      CALL FTOUTP(MSHOW,MFILE,MTERM,CMSG,TTE,1.D2*TTE/TTMG)
      TTO = TTMG-(TTS+TTD+TTP+TTR+TTE)
      CMSG = 'mg time for other                 : TTO    = '
      CALL FTOUTP(MSHOW,MFILE,MTERM,CMSG,TTO,1.D2*TTO/TTMG)
C
      ENDIF
C
C=======================================================================
C     Write solution to file
C=======================================================================
C
      IF (IOUT.GT.0.AND.NEL.LE.100000) THEN
        CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
        MUNIT = 80
        CFILE = 'output.gmv'
C
        CMSG = 'writing solution to file          : CFILE  = '
        CALL FTOUTC(MSHOW,MFILE,MTERM,CMSG,CFILE)
C
        CALL GMVOUT(MUNIT,CFILE,LU)
      ENDIF
C
C=======================================================================
C     Error case
C=======================================================================
C
99997 CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
      CALL FTOUT(MSHOW,MFILE,MTERM,'           FEATFEM DONE!')
      CALL FTOUT(MSHOW,MFILE,MTERM,CLINE)
C
      GOTO 99999
C
99998 IF (MSHOW.GE.2) WRITE(MTERM,*) 'IER', IER
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IN SUBROUTINE ',SUB
C
99999 CLOSE(MFILE)
      WRITE (CFILE,'(I12)') NX
      CFILE = 'results_mgfeat_nx='//TRIM(ADJUSTL(CFILE))//'.log'
      OPEN  (UNIT=66,FILE=CFILE,ACCESS='APPEND')
c$$$      WRITE (66,'(I6,I12,I12,I12,F10.3,I6,F10.3,E10.2)')
c$$$     *            NX,NEL,NVT,NA,TTASM,ITE,TTSOL,DNM
      WRITE (66,'(I6,E16.8,E16.8,E16.8,E16.8,E16.8,E16.8,E16.8,E16.8)')
     *            NX,TTGRID,TTPTR,TTASM,TTRHS,TTBDR,0D0,TTSPMV,TTSOL
      CLOSE (66)
C
C
C
   1  FORMAT(80('-'))
   4  FORMAT(80('%'))
   5  FORMAT(80('$'))
C
      END
C
C
C
************************************************************************
      DOUBLE PRECISION FUNCTION COEF1(X,Y,IA,IB,IBLOC,BFIRST)
************************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION X,Y
      INTEGER IA,IB,IBLOC
      LOGICAL BFIRST
      COEF1 = 1D0
      END
C
C
C
************************************************************************
      DOUBLE PRECISION FUNCTION UE(X,Y)
************************************************************************
      IMPLICIT NONE
      DOUBLE PRECISION X,Y
      UE = 0D0
      END
