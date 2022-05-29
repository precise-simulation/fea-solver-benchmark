************************************************************************
*     Subroutine to output a vector/array to file                      *
************************************************************************
      SUBROUTINE MWRITEV(DVECT,ILEN,CFILE)
      IMPLICIT NONE
      CHARACTER SUB*6,FMT*15,CPARAM*1,CFILE*(*)
      COMMON /CHAR/   SUB,FMT(3),CPARAM(120)
      INTEGER ILEN,MUNIT,IFORMAT,I
      DOUBLE PRECISION DVECT(*)
C
      MUNIT = 71
C
      OPEN   (UNIT=MUNIT,FILE=CFILE)
      REWIND (MUNIT)
C
      IFORMAT = 0
      IF (IFORMAT.EQ.0) THEN
        DO 10 I=1,ILEN
          WRITE (MUNIT,1000) DVECT(I)
10      CONTINUE
      ENDIF
C
      REWIND (MUNIT)
      CLOSE(MUNIT)
C
1000  FORMAT(3D25.16)
      END
C
************************************************************************
*     Subroutine to output a full matrix/array to file                 *
************************************************************************
      SUBROUTINE MWRITEM(DA,NA,NEQ,KCOLA,KLDA,CFILE)
      IMPLICIT NONE
      INTEGER NA,NEQ,KCOLA(NA),KLDA(NEQ+1),MUNIT,IROW,JCOL,IA
      DOUBLE PRECISION DA(NA),DROW(NEQ)
      CHARACTER CFILE*(*),FMTS*25
C
      MUNIT = 71
      WRITE(FMTS,'(A,I12,A)') '(',NEQ,'E16.8)'
C
C     
      OPEN   (UNIT=MUNIT,FILE=CFILE)
      REWIND (MUNIT)
C
      DO 10 IROW=1,NEQ
C
C       Zero entries in DROW
        DO 20 JCOL=1,NEQ
          DROW(JCOL) = 0D0       
20      CONTINUE
C
        DO 30 IA=KLDA(IROW),KLDA(IROW+1)-1
          JCOL = KCOLA(IA)
          DROW(JCOL) = DA(IA)
30      CONTINUE
        WRITE (MUNIT,FMT=FMTS) (DROW(JCOL),JCOL=1,NEQ)
10    CONTINUE      
C
      REWIND (MUNIT)
      CLOSE  (MUNIT)
C
      END
C
************************************************************************
*     Subroutines to output to string, real, integer to log file       *
************************************************************************
      SUBROUTINE FTOUT(MSHOW,MFILE,MTERM,CMSG)
      IMPLICIT NONE
      INTEGER MSHOW,MTERM,MFILE
      CHARACTER CMSG*(*)
C
      IF (MSHOW.GE.0) WRITE(MTERM,*) CMSG
      IF (MSHOW.GE.2) WRITE(MFILE,*) CMSG
C
      END
C
      SUBROUTINE FTOUTC(MSHOW,MFILE,MTERM,CMSG,CIN)
      IMPLICIT NONE
      INTEGER MSHOW,MTERM,MFILE
      CHARACTER CMSG*(*),CIN*(*)
C
      IF (MSHOW.GE.0) WRITE(MTERM,'(A,A)') CMSG,CIN
      IF (MSHOW.GE.2) WRITE(MFILE,*) CMSG//CIN
C
      END
C
      SUBROUTINE FTOUTI(MSHOW,MFILE,MTERM,CMSG,IIN)
      IMPLICIT NONE
      INTEGER MSHOW,MTERM,MFILE,IIN
      CHARACTER CMSG*(*)
C
      IF (MSHOW.GE.0) WRITE(MTERM,'(A,I10)') CMSG,IIN
      IF (MSHOW.GE.2) WRITE(MFILE,*) CMSG,IIN
C
      END
C
      SUBROUTINE FTOUTD(MSHOW,MFILE,MTERM,CMSG,DIN)
      IMPLICIT NONE
      INTEGER MSHOW,MTERM,MFILE
      DOUBLE PRECISION DIN
      CHARACTER CMSG*(*)
C
      IF (MSHOW.GE.0) WRITE(MTERM,'(A,D10.2)') CMSG,DIN
      IF (MSHOW.GE.2) WRITE(MFILE,*) CMSG,DIN
C
      END
C
      SUBROUTINE FTOUTF(MSHOW,MFILE,MTERM,CMSG,DIN)
      IMPLICIT NONE
      INTEGER MSHOW,MTERM,MFILE
      DOUBLE PRECISION DIN
      CHARACTER CMSG*(*)
C
      IF (MSHOW.GE.0) WRITE(MTERM,'(A,F10.2)') CMSG,DIN
      IF (MSHOW.GE.2) WRITE(MFILE,*) CMSG,DIN
C
      END
C
      SUBROUTINE FTOUTP(MSHOW,MFILE,MTERM,CMSG,DIN1,DIN2)
      IMPLICIT NONE
      INTEGER MSHOW,MTERM,MFILE
      DOUBLE PRECISION DIN1,DIN2
      CHARACTER CMSG*(*)
C
      IF (MSHOW.GE.0) WRITE(MTERM,'(A,D10.2,A,F4.1,A)') 
     * CMSG,DIN1,' ( ',DIN2,'% )'
      IF (MSHOW.GE.2) WRITE(MFILE,*)
     * CMSG,DIN1,' ( ',DIN2,'% )'
C
      END
C
C
C 
      SUBROUTINE WVECREF(DX,NEQ,A)
      IMPLICIT NONE
      INTEGER NEQ,I
      DOUBLE PRECISION DX(NEQ),A
C
      DO 1 I=1,NEQ
        DX(I) = A
1     CONTINUE
C
      END
************************************************************************
*                                                                      *
* IC01n                                                                *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using SOR-method      *
*          Double/double precision version                             *
*                                                                      *
* Subroutines/functions called  LLI1, LCL1                             *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DA       R*8    Matrix                                               *
* KCOL     I*4    Pointer vectors for the matrix DA corresponding to   *
* KLD      I*4    the storage technique n                              *
* KOP      I*4                                                         *
* DX       R*8    Starting vector                                      *
* DB       R*8    Right hand side                                      *
* NEQ      I*4    Number of equations                                  *
* NIT      I*4    Maximum number of iterations                         *
* EPS      R*8    Desired precision                                    *
* OMEGA    R*8    Relaxation parameter                                 *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector                                      *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
*                 +1  Precision EPS not achieved after NIT iterations  *
*                                                                      *
************************************************************************
C
      SUBROUTINE IC017_4C(DA,KCOL,KLD,DX,DB,NEQ,NIT,ITE,EPS,OMEGA,N)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      DIMENSION DA(*),KCOL(*),KLD(*),DX(*),DB(*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/
C
      SUB='IC017'
      IF (ICHECK.GE.997) CALL OTRC('IC017 ','01/02/89')
C
      BMSG2=M.GE.2.OR.MT.GE.2
      DM1=0D0
      DM2=1D0
C
      IF (ICHECK.GT.0) THEN
       CALL LLI1(DB,NEQ,RBNORM,IND)
       IF (RBNORM.EQ.0D0) THEN
        CALL LCL1(DX,NEQ)
        IF (BMSG2) CALL OMSG(70,'IC017 ')
        GOTO 99999
       ENDIF
      ENDIF
C
      NSEC   = NEQ/N   ! Number of sections/horizontal rows with N nodes
                       ! equal to the number of rows of vertices in grid
C
      DO 2 ITE=1,NIT
      DM1=0D0
      DM2=0D0

      DO 3 IEQ=1,NEQ
      AUX=DB(IEQ)
      DO 4 ICOL=KLD(IEQ)+1,KLD(IEQ+1)-1
4     AUX=AUX-DA(ICOL)*DX(KCOL(ICOL))
      AUX=OMEGA*(AUX/DA(KLD(IEQ))-DX(IEQ))+DX(IEQ)
      DM1=MAX(DM1,ABS(AUX-DX(IEQ)))
      DM2=MAX(DM2,ABS(AUX))
3     DX(IEQ)=AUX
C
C
C     Loop over odd grid rows
      DO 20 IROW=1,NSEC,2
C
C       Loop over odd grid columns
        DO 25 ICOL=1,N,2
C
          IEQ = ICOL+N*(IROW-1)
C
          AUX = DB(IEQ)
          DO ICOLA=KLD(IEQ)+1,KLD(IEQ+1)-1
            AUX = AUX - DA(ICOLA)*DX(KCOL(ICOLA))
          ENDDO
          AUX     = OMEGA*(AUX/DA(KLD(IEQ)) - DX(IEQ))+DX(IEQ)
          DM1=MAX(DM1,ABS(AUX-DX(IEQ)))
          DM2=MAX(DM2,ABS(AUX))
          DX(IEQ) = AUX
25      CONTINUE
C
C       Loop over even grid columns
        DO 28 ICOL=2,N,2
C
          IEQ = ICOL+N*(IROW-1)
C
          AUX = DB(IEQ)
          DO ICOLA=KLD(IEQ)+1,KLD(IEQ+1)-1
            AUX = AUX - DA(ICOLA)*DX(KCOL(ICOLA))
          ENDDO
          AUX     = OMEGA*(AUX/DA(KLD(IEQ)) - DX(IEQ))+DX(IEQ)
          DM1=MAX(DM1,ABS(AUX-DX(IEQ)))
          DM2=MAX(DM2,ABS(AUX))
          DX(IEQ) = AUX
28      CONTINUE
C
20    CONTINUE
C
C     Loop over even grid rows
      DO 30 IROW=2,NSEC,2
C
C       Loop over odd grid columns
        DO 35 ICOL=1,N,2
C
          IEQ = ICOL+N*(IROW-1)
C
          AUX = DB(IEQ)
          DO ICOLA=KLD(IEQ)+1,KLD(IEQ+1)-1
            AUX = AUX - DA(ICOLA)*DX(KCOL(ICOLA))
          ENDDO
          AUX     = OMEGA*(AUX/DA(KLD(IEQ)) - DX(IEQ))+DX(IEQ)
          DM1=MAX(DM1,ABS(AUX-DX(IEQ)))
          DM2=MAX(DM2,ABS(AUX))
          DX(IEQ) = AUX
35      CONTINUE
C
C       Loop over even grid columns
        DO 38 ICOL=2,N,2
C
          IEQ = ICOL+N*(IROW-1)
C
          AUX = DB(IEQ)
          DO ICOLA=KLD(IEQ)+1,KLD(IEQ+1)-1
            AUX = AUX - DA(ICOLA)*DX(KCOL(ICOLA))
          ENDDO
          AUX     = OMEGA*(AUX/DA(KLD(IEQ)) - DX(IEQ))+DX(IEQ)
          DM1=MAX(DM1,ABS(AUX-DX(IEQ)))
          DM2=MAX(DM2,ABS(AUX))
          DX(IEQ) = AUX
38      CONTINUE
C
30    CONTINUE
C
C
      IF (DM1.LE.EPS*DM2) GOTO 99998
      IF (BMSG2) THEN
       WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
       CALL OMSG(74,'IC017 ')
      ENDIF
2     CONTINUE
C
      IER=1
      WRITE (CPARAM,'(I15,D25.16)') NIT,DM1/DM2
      CALL OMSG(71,'IC017 ')
      CALL OMSG(75,'IC017 ')
      GOTO 99999
C
99998 IER=0
      WRITE (CPARAM,'(I15,D25.16)') ITE,DM1/DM2
      CALL OMSG(75,'IC017 ')
C
99999 END
C
C
C
