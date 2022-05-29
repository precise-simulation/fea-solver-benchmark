************************************************************************
*                                                                      *
* NDFG                                                                 *
*                                                                      *
* Purpose  Determination of the total number of degrees of freedom     *
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
* NDFG     I*4    Total number of degrees of freedom                   *
* IER      I*4    Error indicator                                      *
*                 -120 Wrong value of IELTYP                           *
*                                                                      *
************************************************************************
C
      FUNCTION NDFG(IELTYP)
C
      CHARACTER SUB*6,FMT*15,CPARAM*120
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='NDFG'
      IF (ICHECK.EQ.999) CALL OTRC('NDFG  ','10/27/89')
C
      IER=0
      GOTO  (  1, 10, 20, 30, 40,999,999,999,999,999,
     *       100,110,120,130,140,999,999,999,999,999,
     *       200,210,220,230,999,999,999,999,999,999,
     *       300,310,320,330,340,999,999,999,999,999,
     *       400,999,999,999,999,999,999,999,999,999,
     *       500,510,520,530,999,999,999,999,999,999,
     *       600,610,999,999,999,999,999,999,999,999), IELTYP+1
C
999   CONTINUE
      NDFG=0
      WRITE (CPARAM,'(I15)') IELTYP
      CALL WERR(-120,'NDFG  ')
      GOTO 99999
C
1     CONTINUE
100   CONTINUE
400   CONTINUE
      NDFG=NEL
      GOTO 99999
C
10    CONTINUE
110   CONTINUE
      NDFG=NVT
      GOTO 99999
C
20    CONTINUE
120   CONTINUE
220   CONTINUE
500   CONTINUE
      NDFG=NVT+NMT
      GOTO 99999
C
30    NDFG=NVT*3+NEL
      GOTO 99999
40    NDFG=NVT+2*NMT+NEL
      GOTO 99999
C
130   CONTINUE
230   CONTINUE
330   NDFG=NVT+NMT+NEL
      GOTO 99999
C
140   NDFG=NVT+2*NMT+4*NEL
      GOTO 99999
C
200   CONTINUE
300   CONTINUE
310   CONTINUE
      NDFG=NMT
      GOTO 99999
C
210   NDFG=NVT+NEL
      GOTO 99999
320   NDFG=NMT+NEL
      GOTO 99999
340   NDFG=4*NEL
      GOTO 99999
510   NDFG=NVT*3
      GOTO 99999
520   NDFG=NVT*6+NMT
      GOTO 99999
530   NDFG=NVT*6
      GOTO 99999
C
600   CONTINUE
610   CONTINUE
      NDFG=3*NEL
C
99999 END
