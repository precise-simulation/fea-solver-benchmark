************************************************************************
*                                                                      *
* ZDISP                                                                *
*                                                                      *
* Purpose  Compress and release arrays on DWORK                        *
*                                                                      *
* Subroutines/functions called   ZTYPE                                 *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* ILONG    I*4    ILONG > 0  Number of elements to be kept             *
*                 ILONG = 0  The whole array is deleted                *
* LNR      I*4    Number of the array to be compressed                 *
*                 as set by ZNEW                                       *
* ARR      C*6    Name of the array - for error messages only          *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 -104  Wrong value of LNR                             *
*                 -105  ILONG > length of array LNR                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE ZDISP(ILONG,LNR,ARR)
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
      DATA LNR1/0/
C
      IF (ICHECK.GE.998) CALL OTRC('ZDISP ','01/02/89')
      BMSG2=M.GE.2.OR.MT.GE.2
      IER=0
      CALL ZTYPE(LNR,ITYPE)
C
      IF (ITYPE.EQ.0) THEN
C ***  Warning  LNR not specified ***
       IF (BMSG2) THEN
        WRITE (CPARAM,'(A6)') ARR
        CALL OMSG(52,'ZDISP ')
       ENDIF
       GOTO 99999
      ENDIF
C
      IF (ITYPE.LT.0) THEN
C ***  Wrong value of LNR ***
       WRITE (CPARAM,'(I15,A6)') LNR,ARR
       CALL WERR(-104,'ZDISP ')
       GOTO 99999
      ENDIF
C
C *** There is nothing to do...
      IF (ILONG.EQ.KLEN(LNR)) GOTO 99999
C
      IF (IFLAG.NE.0) THEN
C
C *** Block 1
C *** Last call of ZNEW using ILONG=0
C
       IF (LNR.NE.IFLAG) THEN
C ***   Warning IFLAG<>LNR ***
        WRITE (CPARAM,'(A6,2I15)') ARR,LNR,IFLAG
        CALL OMSG(51,'ZDISP ')
C ***   Use IFLAG instead of the given LNR
        LNR=IFLAG
       ENDIF
C ***  Array LNR is the last array on DWORK
       L1=L(LNR)
C ***  Determine correct length corresponding to datatype
       JLONG=ILONG
       IF (ITYPE.GT.1) THEN
        JLONG=(JLONG+1)/2
        L1=(L1+1)/2
       ENDIF
       IF (JLONG.GT.KLEN8(LNR)) GOTO 300
       IWORK=L1+JLONG-1
       IWMAX=MAX(IWMAX,IWORK)
       IFLAG=0
C
      ELSE
C
C *** Block 2
C *** Last call of ZNEW using ILONG > 0
C
C ***  Determine correct length corresponding to datatype
       JLONG=ILONG
       IF (ITYPE.GT.1) JLONG=(JLONG+1)/2
C ***  Determine correct offset according to data type ***
       ID=KLEN8(LNR)-JLONG
C
       IF (ID.LT.0) GOTO 300
C
C ***  Adjust starting address corresponding to datatype REAL*8
       DO 210 IARR=1,NNARR
         IF (L(IARR).GE.1.AND.KTYPE(IARR).GT.1) THEN
           L(IARR) = (L(IARR)+1)/2
         ENDIF
210    CONTINUE
C
C ***  Find first element of the next array ***
       J0=L(LNR)+KLEN8(LNR)
C ***  Present array is the last array?
       IF (J0.EQ.IWORK+1) GOTO 240
C ***  Revise L ***
       DO 220 IARR=1,NNARR
       IF (L(IARR).GT.L(LNR)) L(IARR)=L(IARR)-ID
220    CONTINUE
C ***  Compress DWORK copying INTEGER*4 elements ***
       J1=J0*2-1
       CALL LCP3(KWORK(J1),KWORK(J1-2*ID),2*(IWORK-J0+1))
c$$$       DO 230 J1=J0*1-1,2*IWORK
c$$$         KWORK(J1-2*ID) = KWORK(J1)
c$$$230    CONTINUE
C
C ***  Determine new value of IWORK ***
240    IWORK=IWORK-ID
C
C ***  Adjust starting address corresponding to datatype
       DO 250 IARR=1,NNARR
       IF (L(IARR).GE.1.AND.KTYPE(IARR).GT.1) L(IARR)=2*L(IARR)-1
250    CONTINUE
C
      ENDIF
C
      KLEN(LNR)=ILONG
      IF (ILONG.EQ.0) THEN
C ***  Delete array LNR if ILONG=0 ***
       L(LNR)=0
       KTYPE(LNR)=0
       KLEN8(LNR)=0
       LNR1=LNR
       LNR=0
      ELSE
       KLEN8(LNR)=JLONG
      ENDIF
C
      IF (ILONG.EQ.0) THEN
       WRITE (CPARAM,'(A6,I15)') ARR,LNR1
       CALL OMSG(3,'ZDISP ')
      ELSE
       WRITE (CPARAM,'(A6,2I15)') ARR,LNR,ILONG
       CALL OMSG(4,'ZDISP ')
      ENDIF
      GOTO 99999
C
C *** Error *** Length > length of array LNR ***
300   WRITE (CPARAM,'(A6)') ARR
      CALL WERR(-105,'ZDISP ')
C
99999 END
