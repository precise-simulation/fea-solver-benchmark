************************************************************************
*                                                                      *
* OERR                                                                 *
*                                                                      *
* Purpose  Writing error messages during runtime                       *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IER      I*4    Number of error message                              *
* SUB0     C*6    Name of calling routine                              *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
*                                                                      *
************************************************************************
C
      SUBROUTINE OERR(IER,SUB0)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER SUB0*6,A1*6,A2*6,A3*60,CFMT*140
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/CHAR/
      SAVE BERR
      DATA BERR/.FALSE./
C
C *** Initialization
      IF (SUB0.EQ.'ZINIT ') THEN
       BERR=IER.EQ.1
C
      ELSE IF (BERR) THEN
C
       CFMT='('' ***'',A7,'' *** ERROR , IER = '',I4,2X,'
       CFMT(140:140)=')'
1      READ (10,'(I5)',ERR=99999) IERR
       IF (IERR.NE.-IER) GOTO 1
       BACKSPACE (10)
       READ (10,'(I5,2I2,I3,A)',ERR=99999) IERR,M0,MT0,IFMT,CFMT(40:139)
       REWIND (10)
       IF (IERR.GT.-IER) GOTO 99999
C
       BTERM=MT.GE.MT0
       BERR0=M.GE.M0.AND.MERR.GE.0
       IF (.NOT.BTERM.AND..NOT.BERR0) GOTO 99999
       IF (MERR.EQ.MTERM) THEN
        BERR0=.FALSE.
        BTERM=.TRUE.
       ENDIF
C
C
       GOTO (10,20,30,40,50,60,70,80,90,100,110,120,130,140) IFMT
       GOTO 99999
C
10     IF (BERR0) WRITE (MERR,CFMT) SUB0,IER
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER
       GOTO 99999
C
20     READ (CPARAM,'(A6,2I15)') A1,I1,I2
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,A1,I1,I2
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,A1,I1,I2
       GOTO 99999
C
30     READ (CPARAM,'(A6,I15)') A1,I1
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,A1,I1
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,A1,I1
       GOTO 99999
C
40     READ (CPARAM,'(A6)') A1
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,A1
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,A1
       GOTO 99999
C
50     READ (CPARAM,'(I15,A6)') I1,A1
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,I1,A1
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,I1,A1
       GOTO 99999
C
60     READ (CPARAM,'(A6,I15,A6,I15)') A1,I1,A2,I2
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,A1,I1,A2,I2
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,A1,I1,A2,I2
       GOTO 99999
C
70     READ (CPARAM,'(A15,I15,A6)') A3,I1,A1
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,A3,I1,A1
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,A3,I1,A1
       GOTO 99999
C
80     READ (CPARAM,'(I15)') I1
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,I1
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,I1
       GOTO 99999
C
90     READ (CPARAM,'(A60,2I15)') A3,I1,I2
       DO 91 IPOS=1,140
91     IF (CFMT(IPOS:IPOS).EQ.'?') GOTO 92
92     WRITE (CFMT(IPOS:IPOS+1),'(I2)') I2
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,A3,I1
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,A3,I1
       GOTO 99999
C
100    READ (CPARAM,'(A60,2I15)') A3,I1,I2
       DO 101 IPOS=1,140
101    IF (CFMT(IPOS:IPOS).EQ.'?') GOTO 102
102    WRITE (CFMT(IPOS:IPOS+1),'(I2)') I2
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,I1,A3
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,I1,A3
       GOTO 99999
C
110    READ (CPARAM,'(2I15)') I1,I2
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,I1,I2
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,I1,I2
       GOTO 99999
C
120    READ (CPARAM,'(I15,D25.16)') I1,D1
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,I1,D1
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,I1,D1
       GOTO 99999
C
130    READ (CPARAM,'(I15,2D25.16)') I1,D1,D2
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,I1,D1,D2
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,I1,D1,D2
       GOTO 99999
C
140    READ (CPARAM,'(D25.16)') D1
       IF (BERR0) WRITE (MERR,CFMT) SUB0,IER,D1
       IF (BTERM) WRITE (MTERM,CFMT) SUB0,IER,D1
       GOTO 99999
C
      ENDIF
C
99999 END
