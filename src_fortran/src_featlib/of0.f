************************************************************************
*                                                                      *
* OF0                                                                  *
*                                                                      *
* Purpose  Open I/O-file                                               *
*                                                                      *
* Subroutines/functions called  None                                   *
*                                                                      *
* Version from  12/11/89                                               *
*                                                                      *
* INPUT   TYPE                                                         *
* -----    ----                                                        *
* MFILE   I*4    Unit number of external file                          *
* CFILE   C*(*)  Filename                                              *
*                CFILE.EQ.'SCRATCH' means temporary file               *
* IFMT    I*4    0  Format free input                                  *
*                1  Formatted input                                    *
*                                                                      *
* OUTPUT  TYPE                                                         *
* ------  ----                                                         *
* MFILE   I*4    Unit number used                                      *
* CFILE   C*(*)  Filename used (length not larger than for input)      *
* IER     I*4    Error indicator                                       *
*                -110  MFILE  exceeds range  16,...,80                 *
*                -112  File  CFILE  could not be opened                *
*                                                                      *
************************************************************************
C
      SUBROUTINE OF0(MFILE,CFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER CFILE*(*),FM*11
      PARAMETER (NNARR=999)
      DIMENSION FM(0:1)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,FM
      DATA FM/'UNFORMATTED','FORMATTED  '/
C
      IF (ICHECK.GE.997) CALL OTRC('OF0   ','12/11/89')
C
      IER=0
      BWARN=.FALSE.
C
      ILEN0=LEN(CFILE)
      IF (ILEN0.GT.60) THEN
       CALL WERR(-109,SUB)
       GOTO 99999
      ENDIF
      CPARAM(1:60)=CFILE
      ILEN=ILEN0
C
C *** Valid unit number ?
C
1     IRETRY=0
2     IF (MFILE.EQ.0) THEN
       CALL OMSG(9,SUB)
       READ (MKEYB,*) MFILE
       READ (MKEYB,'(A)') CPARAM(1:60)
       DO 3 ILEN=1,60
3      IF (CPARAM(ILEN:ILEN).EQ.' ') GOTO 4
4      ILEN=ILEN-1
      ENDIF
      ICLOSE=0
C
      IF ((MFILE.LE.15).OR.(MFILE.GT.80)) THEN
       WRITE (CPARAM,'(I15)') MFILE
       CALL OMSG(30,SUB)
       IF (IRETRY.GT.3) THEN
        WRITE (CPARAM,'(I15)') MFILE
        CALL WERR(-110,SUB)
        GOTO 99999
       ELSE
        IRETRY=IRETRY+1
        MFILE=0
        GOTO 2
       ENDIF
      ENDIF
C
      INQUIRE(UNIT=MFILE,OPENED=BOPEN,NAMED=BNAMED,NAME=CPARAM(61:120))
      IF (ICHECK.GT.0.AND.BOPEN) THEN
       IF (BNAMED.AND.CPARAM(1:60).NE.CPARAM(61:120)) THEN
        BWARN=.TRUE.
        DO 5 ILEN1=61,120
5       IF (CPARAM(ILEN1:ILEN1).EQ.' ') GOTO 6
6       ILEN1=ILEN1-61
        CPARAM(1:60)=CPARAM(61:120)
        WRITE (CPARAM(61:120),'(2I15)') MFILE,ILEN1
        CALL OMSG(32,SUB)
       ELSE IF (.NOT.BNAMED.AND.CPARAM(1:7).NE.'SCRATCH') THEN
        BWARN=.TRUE.
        WRITE (CPARAM,'(I15)') MFILE
        CALL OMSG(33,SUB)
       ENDIF
       IF (BWARN) THEN
        BWARN=.FALSE.
        CALL OMSG(34,SUB)
        READ (MKEYB,*) ICORR
        CALL OMSG(35,SUB)
        READ (MKEYB,*) ICLOSE
        IF (ICLOSE.EQ.1) CLOSE (MFILE)
        IF (ICORR.EQ.1) THEN
         MFILE=0
         GOTO 1
        ENDIF
       ENDIF
      ENDIF
C
100   IF (.NOT.BOPEN.OR.ICLOSE.EQ.1) THEN
       IF (CPARAM(1:7).EQ.'SCRATCH') THEN
        OPEN(UNIT=MFILE,STATUS='SCRATCH',FORM=FM(IFMT),IOSTAT=ICHK)
        IF (ICHK.NE.0) THEN
         WRITE (CPARAM(61:120),'(2I15)') MFILE,ILEN
         CALL WERR(-112,SUB)
         GOTO 99999
        ENDIF
       ELSE
        IF (CPARAM(1:1).EQ.' ') GOTO 99998
        OPEN(UNIT=MFILE,FILE=CFILE,FORM=FM(IFMT),IOSTAT=ICHK)
        IF (ICHK.NE.0) GOTO 99998
       ENDIF
      ENDIF
C
      GOTO 99999
C
C
99998 WRITE (CPARAM(61:120),'(2I15)') MFILE,ILEN
      CALL OMSG(31,SUB)
      IF (IRETRY.GT.3) THEN
       WRITE (CPARAM(61:120),'(2I15)') MFILE,ILEN
       CALL WERR(-112,SUB)
       GOTO 99999
      ENDIF
      IRETRY=IRETRY+1
      CALL OMSG(9,SUB)
      READ(MKEYB,'(I15)') MFILE
      READ(MKEYB,'(A)') CPARAM(1:60)
      GOTO 100
C
99999 IF (ILEN.GT.ILEN0) CALL OMSG(36,SUB)
      CFILE=CPARAM(1:ILEN)
C
      END
