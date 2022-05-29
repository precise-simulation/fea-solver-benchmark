************************************************************************
*                                                                      *
* NDFL                                                                 *
*                                                                      *
* Purpose  Determination of the number of degrees of freedom           *
*          per element                                                 *
*                                                                      *
* Subroutines/functions called   none                                  *
*                                                                      *
* Version from  10/27/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IELTYP   I*4    Type of element corresponding to  Ennn               *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* NDFL     I*4    Number of degrees of freedom per element             *
* IER      I*4    Error indicator                                      *
*                 -120 Wrong value of IELTYP                           *
*                                                                      *
************************************************************************
C
      FUNCTION NDFL(IELTYP)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='NDFL'
      IF (ICHECK.EQ.999) CALL OTRC('NDFL  ','10/27/89')
      IER=0
C
      GO TO (  1, 10, 20, 30, 40,999,999,999,999,999,
     *       100,110,120,130,140,999,999,999,999,999,
     *       200,210,220,230,999,999,999,999,999,999,
     *       300,310,320,330,340,999,999,999,999,999,
     *       400,999,999,999,999,999,999,999,999,999,
     *       500,510,520,530,999,999,999,999,999,999,
     *       600,610,999,999,999,999,999,999,999,999), IELTYP+1
C
999   CONTINUE
      NDFL=0
      WRITE (CPARAM,'(I15)') IELTYP
      CALL WERR(-120,'NDFL  ')
      GOTO 99999
C
1     CONTINUE
100   CONTINUE
400   CONTINUE
      NDFL=1
      GOTO 99999
C
10    CONTINUE
110   CONTINUE
      NDFL=NVE
      GOTO 99999
C
20    CONTINUE
120   CONTINUE
220   CONTINUE
500   CONTINUE
      NDFL=2*NVE
      GOTO 99999
C
30    CONTINUE
40    CONTINUE
      NDFL=10
      GOTO 99999
C
130   CONTINUE
230   CONTINUE
330   NDFL=2*NVE+1
      GOTO 99999
C
140   NDFL=16
      GOTO 99999
C
200   CONTINUE
300   CONTINUE
310   CONTINUE
      NDFL=NVE
      GOTO 99999
C
210   CONTINUE
340   NDFL=4
      GOTO 99999
320   NDFL=5
      GOTO 99999
C
510   CONTINUE
      NDFL=3*NVE
      GOTO 99999
C
520   NDFL=21
      GOTO 99999
530   NDFL=18
      GOTO 99999
C
600   CONTINUE
610   CONTINUE
      NDFL=3
C
99999 END
