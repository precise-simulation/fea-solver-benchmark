************************************************************************
*                                                                      *
* ZINIT                                                                *
*                                                                      *
* Purpose  Initialization of COMMON blocks                             *
*                                                                      *
* Subroutines/functions called   None                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NNWORK   I*4    Specified length of DWORK                            *
* CMSG     C*n    File containing messages (Default: FEAT.MSG)         *
* CERR     C*n    File linked to MERR      (Default: FEAT.ERR)         *
* CPRT     C*n    File linked to MPROT     (Default: FEAT.PRT)         *
* CSYS     C*n    File linked to MSYS      (Default: FEAT.SYS)         *
* CTRC     C*n    File linked to MTRC      (Default: FEAT.TRC)         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZINIT(NNWORK,CMSG,CERR,CPRT,CSYS,CTRC)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER*(*) CMSG,CERR,CPRT,CSYS,CTRC,BLANK
      PARAMETER (NNARR=999)
      PARAMETER (BLANK=' ')
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/CHAR/
      SAVE BFIRST
      DATA BFIRST/.TRUE./
C
C *** This step is performed only once
C
      IF (BFIRST) THEN
       NWORK=NNWORK
       IWORK=0
       IWMAX=0
       DO 10 IARR=1,NNARR
10     L(IARR)=0
       DTIMEH=-1D0
       CALL ZTIME(DTIMEH)
C
      ELSE
C
       CLOSE (10)
       CLOSE (MERR)
       CLOSE (MPROT)
       CLOSE (MSYS)
       CLOSE (MTRC)
C
      ENDIF
C
      IF (MTRC.NE.MTERM.AND.MTRC.NE.MERR.AND.MTRC.NE.MPROT
     *    .AND.MTRC.NE.MSYS) THEN
       IF (CTRC(1:1).EQ.BLANK) THEN
        OPEN (MTRC,FILE='FEAT.TRC',ERR=100)
       ELSE
        OPEN (MTRC,FILE=CTRC,ERR=100)
       ENDIF
      ENDIF
C
      IF (MERR.NE.MTERM.AND.MERR.NE.MPROT) THEN
       IF (CERR(1:1).EQ.BLANK) THEN
        OPEN (MERR,FILE='FEAT.ERR',ERR=100)
       ELSE
        OPEN (MERR,FILE=CERR,ERR=100)
       ENDIF
      ENDIF
C
      IF (MSYS.NE.MTERM.AND.MSYS.NE.MERR.AND.MSYS.NE.MPROT) THEN
       IF (CSYS(1:1).EQ.BLANK) THEN
        OPEN (MSYS,FILE='FEAT.SYS',ERR=100)
       ELSE
        OPEN (MSYS,FILE=CSYS,ERR=100)
       ENDIF
      ENDIF
C
      IF (MPROT.NE.MTERM) THEN
       IF (CPRT(1:1).EQ.BLANK) THEN
        OPEN (MPROT,FILE='FEAT.PRT',ERR=100)
       ELSE
        OPEN (MPROT,FILE=CPRT,ERR=100)
       ENDIF
      ENDIF
C
      IF (CMSG(1:1).EQ.BLANK) THEN
       OPEN (10,FILE='FEAT.MSG',STATUS='OLD',ERR=100)
      ELSE
       OPEN (10,FILE=CMSG,STATUS='OLD',ERR=100)
      ENDIF
      WRITE (CPARAM,'(I15)') 9
      CALL OMSG(1,'ZINIT ')
      CALL OERR(1,'ZINIT ')
      CALL OTRC('ZINIT ','T')
C
      GOTO 99999
C
100   CALL OMSG(0,'ZINIT ')
      CALL OERR(0,'ZINIT ')
      CALL OTRC('ZINIT ','F')
      WRITE (MTERM,10000)
C
10000 FORMAT(' *** ZINIT  ***  ERROR DURING INITIALIZATION,',
     *                         ' NO MESSAGES DISPLAYED')
C
99999 BFIRST=.FALSE.
      END
