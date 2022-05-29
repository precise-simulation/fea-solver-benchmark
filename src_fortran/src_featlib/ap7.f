************************************************************************
*                                                                      *
* XAP7                                                                 *
*                                                                      *
* Purpose  Call AP7                                                    *
*          Allocate KLD and KCOL on DWORK                              *
*                                                                      *
* Subroutines/functions called  AP7, EA00, ZNEW, ZDISP, ZFREE          *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions  *
*                 for dummy call only                                  *
* For the description of the remaining input parameter see AP7         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* LCOL     I*4    Number of KCOL                                       *
* LLD      I*4    Number of KLD                                        *
* NA       I*4    Length of KCOL (and VA (DA))                         *
* NEQ      I*4    Number of equations                                  *
*                                                                      *
************************************************************************
*                                                                      *
* AP7                                                                  *
*                                                                      *
* Purpose  Calculation of the pointer vectors KLD and KCOL             *
*          for a matrix corresponding to a given element type          *
*          Storage technique 7/8                                       *
*                                                                      *
* Subroutines/functions called  NDFL, NDFGL, EA00                      *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KCOL1    I*4    Auxiliary array                                      *
* KIND     I*4    Auxiliary array                                      *
* NA       I*4    Maximal length of KCOL                               *
* NEQ      I*4    Number of equations                                  *
* ELE      SUBR   EXTERNAL Subroutine - evaluation of basis functions- *
*                 for dummy call only                                  *
* ISYMM    I*4    >=1 matrix symmetric - storage technique 8           *
* KVERT    I*4    Arrays describing the triangulation                  *
* KMID     I*4                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KCOL     I*4    Pointer vector containing column indices             *
* KLD      I*4    Pointer vector to diagonal elements                  *
* NA       I*4    Effective number of elements in KCOL                 *
* IER      I*4    Error indicator                                      *
*                 -118  Not enough space for KCOL                      *
*                       Error occured on element IEL                   *
*                                                                      *
************************************************************************
C
      SUBROUTINE XAP7(LCOL,LLD,NA,NEQ,ELE,ISYMM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      DIMENSION VWORK(1),KWORK(1)
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      EXTERNAL ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XAP7'
      IF (ICHECK.GE.997) CALL OTRC('XAP7  ','12/11/89')
C
C *** Determine total number of degrees of freedom
      CALL EA00(ELE,ELE,NVE,IELTYP)
      NEQ=NDFG(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
C *** Allocate KLD on DWORK ***
      CALL ZNEW(NEQ+1,-3,LLD,'KLD   ')
      IF (IER.NE.0) GOTO 99999
C
C *** Determine free space on DWORK
      IWMAX0=IWMAX
      CALL ZFREE(3,IFREE)
      IF (IER.NE.0) GOTO 99999
C
      NA=IFREE/3-3
      CALL ZNEW(NA,-3,LCOL,'KCOL  ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LCOL1,'KCOL1 ')
      IF (IER.NE.0) GOTO 99999
      CALL ZNEW(NA,-3,LIND,'KIND  ')
      IF (IER.NE.0) GOTO 99999
C
      CALL AP7(KWORK(L(LCOL)),KWORK(L(LCOL1)),
     *         KWORK(L(LIND)),KWORK(L(LLD)),
     *         NA,NEQ,ELE,ISYMM, KWORK(L(LVERT)),KWORK(L(LMID)))
      IF (IER.NE.0) GOTO 99999
C
C *** Release space on DWORK not used for KCOL ***
      CALL ZDISP(NA,LIND,'KIND  ')
      CALL ZDISP(NA,LCOL1,'KCOL1 ')
      CALL ZDISP(NA,LCOL,'KCOL  ')
C
      IWMAX=MAX(IWORK,IWMAX0)
C
      CALL ZDISP(0,LIND,'KIND  ')
      CALL ZDISP(0,LCOL1,'KCOL1 ')
99999 END
C
C
C
      SUBROUTINE AP7(KCOL,KCOL1,KIND,KLD,NA,NEQ,ELE,ISYMM,KVERT,KMID)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNVE=4)
      DIMENSION KLD(*),KDFG(NNBAS),KDFL(NNBAS)
      DIMENSION KCOL(*),KCOL1(*),KIND(*),KVERT(NNVE,*),KMID(NNVE,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      EXTERNAL ELE
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='AP7'
      IF (ICHECK.GE.997) CALL OTRC('AP7   ','12/11/89')
C
      IER=0
C
      IF (NA.LT.NEQ) THEN
C *** DWORK exhausted
       WRITE (CPARAM,'(I15)') 1
       CALL WERR(-118,'AP7   ')
       GOTO 99999
      ENDIF
      IFREE=NA-NEQ
      NA=NEQ
C *** Initialization of KIND and KCOL1
      DO 10 IEQ=1,NEQ
      KIND(IEQ)=0
10    KCOL1(IEQ)=IEQ
C
C *** Set element number by dummy call ***
      CALL EA00(ELE,ELE,NVE,IELTYP)
C *** Determine number of degrees of freedom per element ***
      IDFL=NDFL(IELTYP)
      IF (IER.NE.0) GOTO 99999
C
      BSYMM=ISYMM.EQ.1
      IDOFE1=1
C *** Loop over elements
      DO 100 IEL=1,NEL
C
C *** KDFG returns the global degrees of freedom in increasing order
      CALL NDFGL(IEL,0,IELTYP,KVERT,KMID,KDFG,KDFL)
      IF (IER.NE.0) GOTO 99999
C
C *** Loop over local number of degrees of freedom
      DO 110 IDOFE=1,IDFL
      IROW=KDFG(IDOFE)
      IF (BSYMM) THEN
       IF (IROW.EQ.NEQ) GOTO 110
       IDOFE1=IDOFE
      ENDIF
C
C *** Loop over off-diagonal elements
C *** only upper triangular part in symmetric case
      DO 120 JDOFE=IDOFE1,IDFL
      IF (IDOFE.EQ.JDOFE) GOTO 120
      JCOL=KDFG(JDOFE)
C *** JCOL is the global number of d.o.f. to be inserted into row IROW
C
C *** Look whether entry (IROW,JCOL) is already provided
      IPOS=IROW
121   IF (KIND(IPOS).EQ.0) THEN
C *** final entry in row IROW up to now
C
      IF (IFREE.LE.0) THEN
C *** DWORK exhausted
        WRITE (CPARAM,'(I15)') IEL
        CALL WERR(-118,'AP7   ')
        GOTO 99999
       ENDIF
C
C *** insert new entry
       NA=NA+1
       IFREE=IFREE-1
       KCOL1(NA)=JCOL
       KIND(IPOS)=NA
       KIND(NA)=0
       GOTO 120
C
      ELSE
C
       IPOS=KIND(IPOS)
       IF (KCOL1(IPOS).NE.JCOL) GOTO 121
C
      ENDIF
C
120   CONTINUE
110   CONTINUE
100   CONTINUE
C
C
C *** Collect entries on KCOL1 separately for each row
C
      NA=0
      DO 200 IEQ=1,NEQ
      NA=NA+1
      KLD(IEQ)=NA
      KCOL(NA)=IEQ
      IPOS=IEQ
C
201   IF (KIND(IPOS).NE.0) THEN
       NA=NA+1
       IPOS=KIND(IPOS)
       KCOL(NA)=KCOL1(IPOS)
       GOTO 201
      ENDIF
C
200   CONTINUE
      KLD(NEQ+1)=NA+1
C
C
C *** Sort off-diagonal entries on KCOL separately for each row
C
      DO 300 IEQ=1,NEQ
C
301   BSORT=.TRUE.
      DO 302 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-2
      IF (KCOL(ICOL).GT.KCOL(ICOL+1)) THEN
       IHELP=KCOL(ICOL)
       KCOL(ICOL)=KCOL(ICOL+1)
       KCOL(ICOL+1)=IHELP
       BSORT=.FALSE.
      ENDIF
302   CONTINUE
      IF (.NOT.BSORT) GOTO 301
C
300   CONTINUE
C
99999 END
