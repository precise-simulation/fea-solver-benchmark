************************************************************************
*                                                                      *
* ZCTYPE                                                               *
*                                                                      *
* Purpose  Type conversion of array on DWORK                           *
*                                                                      *
* Subroutines/functions called    ZTYPE, ZNEW, ZDISP                   *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ITYPE    I*4    Desired new type                                     *
* LNR      I*4    Number of array                                      *
* ARR      C*6    Name of array                                        *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 -104  Value of LNR invalid                           *
*                 -101  Value of ITYPE invalid                         *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZCTYPE(ITYPE,LNR,ARR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      DIMENSION VWORK(1),KWORK(1)
C
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /TABLE/  KTYPE(NNARR),KLEN(NNARR),KLEN8(NNARR),IFLAG
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TABLE/
C
      IF (ICHECK.GE.998) CALL OTRC('ZCTYPE','01/02/89')
      IER=0
C
      CALL ZTYPE(LNR,JTYPE)
      IF (JTYPE.LE.0) THEN
       WRITE (CPARAM,'(I15,A6)') LNR,ARR
       CALL WERR(-104,'ZCTYPE')
       GOTO 99999
      ENDIF
C
      IF ((ITYPE.LE.0).OR.(ITYPE.GE.4)) THEN
       WRITE (CPARAM,'(I15,A6)') ITYPE,ARR
       CALL WERR(-101,'ZCTYPE')
       GOTO 99999
      ENDIF
C
      ICOUNT=KLEN(LNR)
      ILNR= L(LNR)
      ILNR4=(L(LNR)-1)*2+1
C
C
      GOTO (100,200,300),JTYPE
C
100   GOTO (99997,120,130),ITYPE
C
120   DO 121 I=0,ICOUNT-1
121   VWORK(ILNR4+I)=REAL(DWORK(ILNR+I))
      GOTO 140
130   DO 131 I=0,ICOUNT-1
131   KWORK(ILNR4+I)=INT(DWORK(ILNR+I))
140   L(LNR)=ILNR4
      KTYPE(LNR)=ITYPE
      KLEN(LNR)=2*KLEN(LNR)
      CALL ZDISP(ICOUNT,LNR,ARR)
      GOTO 99998
C
200   GOTO (210,99997,230),ITYPE
C
210   CALL ZNEW(ICOUNT,ITYPE,LNX,'HELP  ')
      IF (IER.NE.0) GOTO 99999
      ILNX=L(LNX)
      DO 211 I=0,ICOUNT-1
211   DWORK(ILNX+I)=DBLE(VWORK(ILNR+I))
      CALL ZDISP(0,LNR,'HELP  ')
      LNR=LNX
      GOTO 99998
230   DO 233 I=0,ICOUNT-1
233   KWORK(ILNR+I)=INT(VWORK(ILNR+I))
      GOTO 99998
C
300   GOTO (310,320,99997),ITYPE
C
310   CALL ZNEW(ICOUNT,ITYPE,LNX,'HELP  ')
      IF (IER.NE.0) GOTO 99999
      ILNX=L(LNX)
      DO 311 I=0,ICOUNT-1
311   DWORK(ILNX+I)=DBLE(KWORK(ILNR+I))
      CALL ZDISP(0,LNR,'HELP  ')
      LNR=LNX
      GOTO 99998
320   DO 322 I=0,ICOUNT-1
322   VWORK(ILNR+I)=REAL(KWORK(ILNR+I))
      GOTO 99998
C
99997 WRITE (CPARAM,'(A6,I15)') ARR,LNR
      CALL OMSG(53,'ZCTYPE')
      GOTO 99999
C
99998 KTYPE(LNR)=ITYPE
      WRITE (CPARAM,'(A6,2I15)') ARR,LNR,ITYPE
      CALL OMSG(6,'ZCTYPE')
C
99999 END
