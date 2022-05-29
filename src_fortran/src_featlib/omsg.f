************************************************************************
*                                                                      *
* OMSG                                                                 *
*                                                                      *
* Purpose  Writing notes during runtime                                *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IMSG     I*4    Number of message                                    *
* SUB0     C*6    Name of calling routine                              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE OMSG(IMSG,SUB0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER SUB0*6,A1*6,A2*6,A3*60,CFMT*140,CFMT0*140,CLAST*140
C
      PARAMETER (NNFMT=9,NNMSG=99)
      DIMENSION CFMT(NNFMT),IFMT(NNFMT),KMT(NNMSG),KM(NNMSG)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/CHAR/
      SAVE IFMT,CFMT,KM,KMT,BMSG,CLAST,LASTF,LASTMS,LASTM,LASTMT,NMSG
      DATA IFMT/NNFMT*0/,BMSG/.FALSE./,LASTMS/0/
      DATA KM/NNMSG*-99/,KMT/NNMSG*-99/,NMSG/0/,LASTF/0/
C
C *** Initialization
      IF (SUB0.EQ.'ZINIT ') THEN
       BMSG=IMSG.EQ.1
       IF (BMSG) THEN
        READ (CPARAM,'(I15)') NFMT
        NMSG=MIN(NFMT,NNFMT)
        DO 2 JFMT=1,NMSG
        CFMT(JFMT)='('' ***'',A7,'' *** '','
        CFMT(JFMT)(140:140)=')'
2       READ (10,'(I5,2I2,I3,A)',ERR=99999)
     *        I,KM(JFMT),KMT(JFMT),IFMT(JFMT),CFMT(JFMT)(20:139)
        CLAST='('' ***'',A7,'' *** '','
        CLAST(140:140)=')'
C
        DO 3 JFMT=NMSG+1,NNMSG
        READ (10,'(I5,2I2)',ERR=99999) IMSG0,M0,MT0
        IF (IMSG0.GT.NNMSG) THEN
         REWIND (10)
         GOTO 99999
        ENDIF
        KM(IMSG0)=M0
3       KMT(IMSG0)=MT0
C
      ENDIF
C
      ELSE IF (BMSG) THEN
C
       IF (KM(IMSG).EQ.-99) GOTO 99999
       IF (IMSG.GT.6) THEN
        MOUT=MPROT
       ELSE
        MOUT=MSYS
       ENDIF
       BTERM=MT.GE.KMT(IMSG)
       BPROT=M.GE.KM(IMSG).AND.MOUT.GE.0
       IF (.NOT.BTERM.AND..NOT.BPROT) GOTO 99999
       IF (MOUT.EQ.MTERM) THEN
        BPROT=.FALSE.
        BTERM=.TRUE.
       ENDIF
C
       IF (IMSG.LE.NMSG) THEN
        IFMT0=IFMT(IMSG)
        CFMT0=CFMT(IMSG)
       ELSE IF (IMSG.EQ.LASTMS) THEN
        IFMT0=LASTF
        CFMT0=CLAST
       ELSE
        IF (LASTMS.GT.IMSG) REWIND (10)
1       READ (10,'(I5)',ERR=99999) LASTMS
        IF (LASTMS.LT.IMSG) GOTO 1
        BACKSPACE (10)
        READ (10,'(I5,2I2,I3,A)',ERR=99999)
     *        LASTMS,LASTM,LASTMT,LASTF,CLAST(20:139)
        IF (LASTMS.GT.IMSG) GOTO 99999
        IFMT0=LASTF
        CFMT0=CLAST
       ENDIF
C
       GOTO (10,20,30,40,50,60,70,80,90,100,110,120,130,140) IFMT0
       GOTO 99999

10     IF (BPROT) WRITE (MOUT,CFMT0) SUB0
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0
       GOTO 99999

20     READ (CPARAM,'(A6,2I15)') A1,I1,I2
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,A1,I1,I2
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,A1,I1,I2
       GOTO 99999
C
30     READ (CPARAM,'(A6,I15)') A1,I1
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,A1,I1
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,A1,I1
       GOTO 99999
C
40     READ (CPARAM,'(A6)') A1
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,A1
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,A1
       GOTO 99999
C
50     READ (CPARAM,'(I15,A6)') I1,A1
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,I1,A1
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,I1,A1
       GOTO 99999
C
60     READ (CPARAM,'(A6,I15,A6,I15)') A1,I1,A2,I2
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,A1,I1,A2,I2
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,A1,I1,A2,I2
       GOTO 99999
C
70     CONTINUE
C       READ (CPARAM,'(A15,I15,A6)') A3,I1,A1
C       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,A3,I1,A1
C       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,A3,I1,A1
C       GOTO 99999
C
80     READ (CPARAM,'(I15)') I1
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,I1
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,I1
       GOTO 99999
C
90     READ (CPARAM,'(A60,2I15)') A3,I1,I2
       DO 91 IPOS=1,140
91     IF (CFMT0(IPOS:IPOS).EQ.'?') GOTO 92
92     WRITE (CFMT0(IPOS:IPOS+1),'(I2)') I2
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,A3(1:I2),I1
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,A3(1:I2),I1
       CFMT0(IPOS:IPOS+1)='??'
       GOTO 99999
C
100    READ (CPARAM,'(A60,2I15)') A3,I1,I2
       DO 101 IPOS=1,140
101    IF (CFMT0(IPOS:IPOS).EQ.'?') GOTO 102
102    WRITE (CFMT0(IPOS:IPOS+1),'(I2)') I2
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,I1,A3
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,I1,A3
       CFMT0(IPOS:IPOS+1)='??'
       GOTO 99999
C
110    READ (CPARAM,'(2I15)') I1,I2
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,I1,I2
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,I1,I2
       GOTO 99999
C
120    write (*,*) 'CPARAM = ',CPARAM
       READ (CPARAM,'(I15,D25.16)') I1,D1
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,I1,D1
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,I1,D1
       GOTO 99999
C
130    READ (CPARAM,'(I15,2D25.16)') I1,D1,D2
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,I1,D1,D2
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,I1,D1,D2
       GOTO 99999
C
140    READ (CPARAM,'(D25.16)') D1
       IF (BPROT) WRITE (MOUT,CFMT0) SUB0,D1
       IF (BTERM) WRITE (MTERM,CFMT0) SUB0,D1
       GOTO 99999
C
      ENDIF
C
99999 END
