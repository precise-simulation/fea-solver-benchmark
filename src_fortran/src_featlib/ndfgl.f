************************************************************************
*                                                                      *
* NDFGL                                                                *
*                                                                      *
* Purpose  Determination of the global degrees of freedom              *
*          on a given element                                          *
*                                                                      *
* Subroutines/functions called  NGLS                                   *
*                                                                      *
* Version from  10/27/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IEL      I*4    Number of element                                    *
* IPAR     I*4    Controls sorting of the resulting vectors            *
*                 (see below)                                          *
* IELTYP   I*4    Type of element corresponding to  Ennn               *
* KVERT    I*4    Array of vertices of the elements                    *
* KMID     I*4    Array of midpoints (or edges) of the elements        *
*                 if NMT > 0                                           *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KDFG     I*4    Global degrees of freedom on element IEL             *
*                  KDFG is sorted if  IPAR >= 0                        *
* KDFL     I*4    Local degrees od freedom corresponding to KDFG       *
*                  KDFL is determined only if  IPAR = 1                *
* IER      I*4    Error indicator                                      *
*                 -120 Wrong value of IELTYP                           *
*                                                                      *
************************************************************************
C
      SUBROUTINE NDFGL(IEL,IPAR,IELTYP,KVERT,KMID,KDFG,KDFL)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNVE=4,NNVME=8)
      DIMENSION KDFG(*),KDFL(*),JVG(NNVME),JVL(NNVME)
      DIMENSION KVERT(NNVE,*),KMID(NNVE,*)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/
C
      SUB='NDFGL'
      IF (ICHECK.GE.999) CALL OTRC('NDFGL ','10/27/89')
C
      IER=0
C
C *** Constant elements with 1 dof per element
C
      IF ((IELTYP.EQ.0).OR.(IELTYP.EQ.10).OR.(IELTYP.EQ.40)) THEN
       KDFG(1)=IEL
       IF (IPAR.EQ.1) KDFL(1)=1
       GOTO 99999
      ENDIF
C
      DO  1 IVE=1,NVE
      JVG(IVE)=KVERT(IVE,IEL)
1     JVL(IVE)=IVE
      NKE=NVE
      IF (NMT.GT.0) THEN
       DO 2 IME=1,NVE
       JVG(NVE+IME)=KMID(IME,IEL)
2      JVL(NVE+IME)=NVE+IME
       NKE=NKE+NVE
      ENDIF
      IF (IPAR.GE.0) CALL NGLS(JVG,JVL,NKE)
C
      GO TO ( 10, 20, 30, 40,999,999,999,999,999,
     *       999,110,120,130,140,999,999,999,999,999,
     *       200,210,220,230,999,999,999,999,999,999,
     *       300,310,320,330,340,999,999,999,999,999,
     *       999,999,999,999,999,999,999,999,999,999,
     *       500,510,520,530,999,999,999,999,999,999,
     *       600,610,999,999,999,999,999,999,999,999), IELTYP
C
999   CONTINUE
      WRITE (CPARAM,'(I15)') IELTYP
      CALL WERR(-120,'NDFGL ')
      GOTO 99999
C
C *** Lagrange type elements and Morley element
C
10    CONTINUE
20    CONTINUE
110   CONTINUE
120   CONTINUE
220   CONTINUE
500   CONTINUE
      DO 11 IKE=1,NKE
11    KDFG(IKE)=JVG(IKE)
      IF (IPAR.EQ.1) THEN
       DO 12 IKE=1,NKE
12     KDFL(IKE)=JVL(IKE)
      ENDIF
      GOTO 99999
C
C *** Cubic Hermite element
C
30    DO 31 IDFL=1,7,3
      J1=(IDFL+2)/3
      KDFG(IDFL)=JVG(J1)*3-2
      KDFG(IDFL+1)=KDFG(IDFL)+1
31    KDFG(IDFL+2)=KDFG(IDFL)+2
C *** Dof corresponding to the center of gravity
      KDFG(10)=3*NVT+IEL
      IF (IPAR.EQ.1) THEN
       DO 32 IDFL=1,7,3
       J1=(IDFL+2)/3
       KDFL(IDFL)=JVL(J1)*3-2
       KDFL(IDFL+1)=KDFL(IDFL)+1
32     KDFL(IDFL+2)=KDFL(IDFL)+2
C ***  Dof corresponding to the center of gravity
       KDFL(10)=10
      ENDIF
      GOTO 99999
C
C *** Cubic Lagrange element
C
40    DO 41 IDFL=1,3
41    KDFG(IDFL)=JVG(IDFL)
      DO 42 IDFL=4,6
      J1=2*IDFL-4
      KDFG(J1)=2*JVG(IDFL)-NVT-1
42    KDFG(J1+1)=2*JVG(IDFL)-NVT
C *** Dof corresponding to the center of gravity
      KDFG(10)=NVT+2*NMT+IEL
      IF (IPAR.EQ.1) THEN
       DO 43 IDFL=1,3
43     KDFL(IDFL)=JVL(IDFL)
       DO 44 IDFL=4,6
       J1=2*IDFL-4
       IVT1=KVERT(JVL(IDFL)-3,IEL)
       IDFL1=JVL(IDFL)-2
       IF (IDFL1.EQ.4) IDFL1=1
       IVT2=KVERT(IDFL1,IEL)
       IF (IVT1.LT.IVT2) THEN
        KDFL(J1)=2*JVL(IDFL)-4
        KDFL(J1+1)=KDFL(J1)+1
       ELSE
        KDFL(J1+1)=2*JVL(IDFL)-4
        KDFL(J1)=KDFL(J1+1)+1
       ENDIF
44     CONTINUE
C ***  Dof corresponding to the center of gravity
       KDFL(10)=10
      ENDIF
      GOTO 99999
C
C *** Biquadratic and piecewise bilinear elements with 9 dof
C *** Augmented quadratic element (P2+bulb)
C
130   CONTINUE
230   CONTINUE
330   CONTINUE
      DO 331 IKE=1,NKE
331   KDFG(IKE)=JVG(IKE)
C *** Dof corresponding to the center of gravity
      KDFG(2*NVE+1)=NVT+NMT+IEL
      IF (IPAR.EQ.1) THEN
       DO 332 IKE=1,NKE
332    KDFL(IKE)=JVL(IKE)
C ***  Dof corresponding to the center of gravity
       KDFL(2*NVE+1)=2*NVE+1
      ENDIF
      GOTO 99999
C
C *** Bicubic Lagrange element
C
140   DO 141 IDFL=1,4
141   KDFG(IDFL)=JVG(IDFL)
      DO 142 IDFL=5,8
      J1=2*IDFL-5
      KDFG(J1)=2*JVG(IDFL)-NVT-1
142   KDFG(J1+1)=2*JVG(IDFL)-NVT
C *** Dofs corresponding to interior nodes
      DO 143 IDFL=1,4
143   KDFG(12+IDFL)=NVT+2*NMT+4*(IEL-1)+IDFL
      IF (IPAR.EQ.1) THEN
       DO 144 IDFL=1,4
144    KDFL(IDFL)=JVL(IDFL)
       DO 145 IDFL=5,8
       J1=2*IDFL-5
       IVT1=KVERT(JVL(IDFL)-4,IEL)
       IDFL1=JVL(IDFL)-3
       IF (IDFL1.EQ.5) IDFL1=1
       IVT2=KVERT(IDFL1,IEL)
       IF (IVT1.LT.IVT2) THEN
        KDFL(J1)=2*JVL(IDFL)-5
        KDFL(J1+1)=KDFL(J1)+1
       ELSE
        KDFL(J1+1)=2*JVL(IDFL)-5
        KDFL(J1)=KDFL(J1+1)+1
       ENDIF
145    CONTINUE
C ***  Dofs corresponding to the interior nodes
       DO 146 IDFL=13,16
146    KDFL(IDFL)=IDFL
      ENDIF
      GOTO 99999
C
C *** Nonconforming elements with NVE dof
C
200   CONTINUE
300   CONTINUE
310   CONTINUE
      DO 201 IME=1,NVE
201   KDFG(IME)=JVG(IME+NVE)-NVT
      IF (IPAR.EQ.1) THEN
       DO 202 IME=1,NVE
202    KDFL(IME)=JVL(IME+NVE)-NVE
      ENDIF
      GOTO 99999
C
C *** Mini element
C
210   CONTINUE
      DO 211 IVE=1,NVE
211   KDFG(IVE)=JVG(IVE)
      KDFG(4)=NVT+IEL
      IF (IPAR.EQ.1) THEN
       DO 212 IVE=1,NVE
212    KDFL(IVE)=JVL(IVE)
      KDFL(4)=4
      ENDIF
      GOTO 99999
C
C *** Nonconforming Han element with 5 dof
C
320   DO 321 IME=1,NVE
321   KDFG(IME)=JVG(IME+NVE)-NVT
C *** Dof corresponding to the center of gravity
      KDFG(5)=NMT+IEL
      IF (IPAR.EQ.1) THEN
       DO 322 IME=1,NVE
322    KDFL(IME)=JVL(IME+NVE)-NVE
C ***  Dof corresponding to the center of gravity
       KDFL(5)=5
      ENDIF
      GOTO 99999
C
C *** 4Q0-element
C
340   CONTINUE
      IEL4=4*(IEL-1)
      DO 341 IDFL=1,4
341   KDFG(IDFL)=IEL4+IDFL
      IF (IPAR.EQ.1) THEN
       DO 342 IDFL=1,4
342    KDFL(IDFL)=IDFL
      ENDIF
      GOTO 99999
C
C *** Cubic Zienkiewicz element
C
510   DO 511 IDFL=1,7,3
      J1=(IDFL+2)/3
      KDFG(IDFL)=JVG(J1)*3-2
      KDFG(IDFL+1)=KDFG(IDFL)+1
511   KDFG(IDFL+2)=KDFG(IDFL)+2
      IF (IPAR.EQ.1) THEN
       DO 512 IDFL=1,7,3
       J1=(IDFL+2)/3
       KDFL(IDFL)=JVL(J1)*3-2
       KDFL(IDFL+1)=KDFL(IDFL)+1
512    KDFL(IDFL+2)=KDFL(IDFL)+2
      ENDIF
      GOTO 99999
C
C *** Argyris element
C
520   DO 521 IDFL=1,16,6
      J1=(IDFL+5)/6
      KDFG(IDFL)=JVG(J1)*6-5
      DO 522 J=1,5
522   KDFG(IDFL+J)=KDFG(IDFL)+J
C *** Dof corresponding to midpoints of edges
521   KDFG(J1+18)=JVG(NVE+J1)+5*NVT
      IF (IPAR.EQ.1) THEN
       DO 523 IDFL=1,16,6
       J1=(IDFL+5)/6
       KDFL(IDFL)=JVL(J1)*6-5
       DO 524 J=1,5
524    KDFL(IDFL+J)=KDFL(IDFL)+J
523    KDFL(J1+18)=JVL(NVE+J1)+18
      ENDIF
      GOTO 99999
C
C *** Reduced Argyris element (Bell element)
C
530   DO 531 IDFL=1,16,6
      L1=(IDFL+5)/6
      KDFG(IDFL)=JVG(L1)*6-5
      DO 532 J=1,5
532   KDFG(IDFL+J)=KDFG(IDFL)+J
531   CONTINUE
      IF (IPAR.EQ.1) THEN
       DO 533 I=1,16,6
       J1=(IDFL+5)/6
       KDFL(IDFL)=JVL(J1)*6-5
       DO 534 J=1,5
534    KDFL(IDFL+J)=KDFL(IDFL)+J
533    CONTINUE
      ENDIF
      GOTO 99999
C
600   CONTINUE
610   CONTINUE
      KDFG(1)=3*IEL-2
      KDFG(2)=KDFG(1)+1
      KDFG(3)=KDFG(1)+2
      IF (IPAR.EQ.1) THEN
       KDFL(1)=1
       KDFL(2)=2
       KDFL(3)=3
      ENDIF
C
99999 END
C
C
C
      SUBROUTINE NGLS(KV1,KV2,IDIM)
C
C *** Bubble sort of the arrays KV1 and KV2
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KV1(*),KV2(*)
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/
C
      IF (ICHECK.EQ.999) CALL OTRC('NGLS  ','03/21/89')
C
      BMORE=.TRUE.
C *** WHILE (BMORE) DO ***
5     IF (.NOT.BMORE) GOTO 99999
       BMORE=.FALSE.
       DO 10 ICOMP=1,IDIM-1
       IF (KV1(ICOMP).GT.KV1(ICOMP+1)) THEN
       JAUX1=KV1(ICOMP)
       JAUX2=KV2(ICOMP)
       KV1(ICOMP)=KV1(ICOMP+1)
       KV2(ICOMP)=KV2(ICOMP+1)
       KV1(ICOMP+1)=JAUX1
       KV2(ICOMP+1)=JAUX2
       BMORE=.TRUE.
       ENDIF
10     CONTINUE
C *** ENDWHILE ***
      GOTO 5
C
99999 END
