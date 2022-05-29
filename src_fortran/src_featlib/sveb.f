************************************************************************
*                                                                      *
* XSVEB                                                                *
*                                                                      *
* Purpose  Allocation of the vectors  KVBD, KEBD, KBCT                 *
*                                                                      *
* Subroutines/functions called  SNVB, SVEB, SVEBS                      *
*                                                                      *
* Version from  04/12/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* IPAR     I*4    -1 Determination of NVBD and KBCT only               *
*                    If NVBD>0 this step is assumed to be performed    *
*                    before                                            *
*                  0 Determination of boundary vertices                *
*                  1 Sorting of boundary vertices and determination    *
*                    of element numbers corresponding to the edges     *
* TMAX     R*8     EXTERNAL FUNCTION - maximum parameter of each       *
*                  boundary component                                  *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* The output quantities are contained in COMMON block /TRIA/           *
*                                                                      *
************************************************************************
*                                                                      *
* SVEB                                                                 *
*                                                                      *
* Purpose  Determination of the numbers of boundary vertices           *
*          and element numbers corresponding to edges if desired       *
*                                                                      *
*                                                                      *
* Subroutines/functions called  SVEBS                                  *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KNPR     I*4    Vector of nodal properties                           *
* KVERT    I*4                                                         *
* KSORT    I*4    Help vector - needed for sorting the boundary nodes  *
* DCORVG   R*8    Cartesian coordinates                                *
* TMAX     R*8    EXTERNAL FUNCTION                                    *
* IPAR     I*4    Vertices are sorted with respect to the parameter    *
*                 if IPAR=1 . Element numbers corresponding to edges   *
*                 are determined                                       *
*                                                                      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KVBD     I*4    Numbers of boundary vertices sorted with respect     *
*                 to the number of boundary components                 *
* KEBD     I*4    Numbers of elements corresponding to edges           *
*                                                                      *
************************************************************************
*                                                                      *
* SVEBS                                                                *
*                                                                      *
* Purpose                                                              *
*                                                                      *
* Subroutines/functions called  None                                   *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KVBD     I*4    Numbers of boundary vertices                         *
* NBD      I*4    Number of points on the boundary                     *
* KSORT    I*4    Help vector - needed for sorting the boundary nodes  *
* DCORVG   R*8    Cartesian coordinates                                *
* H        R*8                                                         *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* KVBD     I*4    Numbers of boundary vertices sorted with respect     *
*                 to the number of boundary components                 *
* KEBD     I*4    Numbers of corresponding elements                    *
*                                                                      *
************************************************************************
*                                                                      *
* SNVB                                                                 *
*                                                                      *
* Purpose  Determination of the total number of boundary vertices      *
*                                                                      *
* Subroutines/functions called  None                                   *
*                                                                      *
* Version from  01/02/89                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* KNPR     I*4    Field of nodal properties                            *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* NVBD     I*4    Total number of boundary vertices                    *
* KBCT     I*4    KBCT(1)=1                                            *
*                 KBCT(IBCT+1)-KBCT(IBCT) = Number of boundary vertices*
*                                           on boundary component IBCT *
*                                                                      *
************************************************************************
C
      SUBROUTINE XSVEB(IPAR,TMAX)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      PARAMETER (NNARR=999)
      INCLUDE 'nnwork.inc'
      DIMENSION VWORK(1),KWORK(1)
      EXTERNAL TMAX
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/TRIAD/,/TRIAA/
C
      SUB='XSVEB'
      IF (ICHECK.GE.997) CALL OTRC('XSVEB ','04/12/91')
C
C *** Determination of the number of boundary vertices
C
      IF (LBCT.EQ.0) THEN
       CALL ZNEW(NBCT+1,-3,LBCT,'KBCT  ')
       IF (IER.NE.0)  GOTO 99999
      ENDIF
      CALL SNVB(KWORK(L(LNPR)),KWORK(L(LBCT)))
      IF (IPAR.EQ.-1.OR.NVBD.EQ.0) GOTO 99999
C
C *** Allocation of the arrays
C
      CALL ZNEW(NVBD,-3,LVBD,'KVBD  ')
      IF (IER.NE.0) GOTO 99999
C
C *** Allocate vector for element numbers if needed
      IF (IPAR.EQ.1) THEN
       CALL ZNEW(NVBD,-3,LEBD,'KEBD  ')
       IF (IER.NE.0) GOTO 99999
      ENDIF
C
      CALL ZNEW(2*NVBD,-3,LSORT,'KSORT ')
      IF (IER.NE.0) GOTO 99999
      CALL SVEB(KWORK(L(LNPR)),KWORK(L(LVERT)),KWORK(L(LADJ)),
     *          KWORK(L(LVBD)),KWORK(L(LEBD)),KWORK(L(LBCT)),
     *          KWORK(L(LSORT)),DWORK(L(LCORVG)),TMAX,IPAR)
C
      CALL ZDISP(0,LSORT,'KSORT ')
C
99999 END
C
C
C
      SUBROUTINE SVEB(KNPR,KVERT,KADJ,KVBD,KEBD,KBCT,KSORT,DCORVG,
     *                TMAX,IPAR)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNVE=4)
      DIMENSION KNPR(*),KVERT(NNVE,*),KADJ(NNVE,*),KVBD(*),KEBD(*),
     *          KSORT(2,*),KBCT(*)
      DIMENSION DCORVG(2,*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /ERRCTL/,/TRIAD/
C
      IF (ICHECK.GE.997) CALL OTRC('SVEB  ','01/02/89')
C
C *** Determination of boundary vertices
C
      DO 10 IBCT=1,NBCT
10    KSORT(1,IBCT)=KBCT(IBCT)
C
      DO 20 IVT=1,NVT
      IBCT=KNPR(IVT)
      IF (IBCT.GT.0) THEN
       ISORT=KSORT(1,IBCT)
       KVBD(ISORT)=IVT
       KSORT(1,IBCT)=ISORT+1
      ENDIF
20    CONTINUE
C
C *** Sorting of vertices
C
      IF (IPAR.GE.1) THEN
C
       DO 50 IBCT=1,NBCT
       IBCT1=KBCT(IBCT)
       IBCT2=KBCT(IBCT+1)-1
       NBD=IBCT2-IBCT1+1
       H=TMAX(IBCT)/NBD
       CALL SVEBS(KVBD(IBCT1),NBD,KSORT,DCORVG,H)
C
C *** Check if KEBD should be generated
C
       IF (IPAR.EQ.2) GOTO 50
C
       DO 51 IVBD=IBCT1,IBCT2
       IVBD1=KVBD(IVBD)
       IF (IVBD.NE.IBCT2) THEN
        IVBD2=KVBD(IVBD+1)
       ELSE
        IVBD2=KVBD(IBCT1)
       ENDIF
C
C *** Find first element of boundary component via extensive search
C
       IF (IVBD.EQ.IBCT1) THEN
C
        DO 52 IEL=1,NEL
        IF (KVERT(4,IEL).EQ.0) THEN
         NVE=3
        ELSE
         NVE=4
        ENDIF
        DO 52 IVE=1,NVE
        IVE1=KVERT(IVE,IEL)
        IVE2=KVERT(MOD(IVE,NVE)+1,IEL)
        IF ((IVBD1.EQ.IVE1).AND.(IVBD2.EQ.IVE2)) THEN
         KEBD(IVBD)=IEL
         GOTO 51
        ENDIF
52      CONTINUE
C
       ELSE
C
C *** Find all other elements using KADJ
C
        DO 53 IVE=1,NVE
        IF (KVERT(IVE,IEL).EQ.IVBD1 .AND.
     *      KVERT(MOD(IVE,NVE)+1,IEL).EQ.IVBD2) THEN
         KEBD(IVBD)=IEL
         GOTO 51
        ENDIF
53      CONTINUE
C
60      DO 54 IVE=1,NVE
        IF (KVERT(IVE,IEL).EQ.IVBD1) THEN
         IEL=KADJ(IVE,IEL)
         IF (KVERT(4,IEL).EQ.0) THEN
          NVE1=3
         ELSE
          NVE1=4
         ENDIF
C
         DO 55 IVE1=1,NVE1
         IF (KVERT(IVE1,IEL).EQ.IVBD1 .AND.
     *       KVERT(MOD(IVE1,NVE1)+1,IEL).EQ.IVBD2) THEN
          KEBD(IVBD)=IEL
          NVE=NVE1
          GOTO 51
         ENDIF
55       CONTINUE
         GOTO 60
        ENDIF
54      CONTINUE
       ENDIF
C
51     CONTINUE
C
50     CONTINUE
      ENDIF
C
      END
C
C
C
      SUBROUTINE SVEBS(KVBD,NBD,KSORT,DCORVG,H)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KVBD(*),KSORT(2,*),DCORVG(2,*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /TRIAD/,/ERRCTL/
C
      IF (ICHECK.GE.998) CALL OTRC('SVEBS ','01/02/89')
C
      DO 1 IBD=1,NBD
1     KSORT(1,IBD)=0
C
C *** Step 1 : How many elements in the interval ((J-1)*H,J*H)
C
      DO 10 IBD=1,NBD
      IVT=KVBD(IBD)
      J=DCORVG(1,IVT)/H+1
10    KSORT(1,J)=KSORT(1,J)+1
      J=0
      DO 20 IBD=1,NBD
      J=J+KSORT(1,IBD)
20    KSORT(1,IBD)=J
C
C *** Step 2 : Insert numbers of elements into KSORT(2,.)
C
      DO 30 IBD=1,NBD
      IVT=KVBD(IBD)
      J=DCORVG(1,IVT)/H+1
      J1=KSORT(1,J)
      KSORT(1,J)=J1-1
30    KSORT(2,J1)=IVT
C
C *** Step 3 : Bubble sort
C
      DO 40 IBD=1,NBD
40    KVBD(IBD)=KSORT(2,IBD)
49    BEND=.TRUE.
      DO 50 IBD=1,NBD-1
      J1=KVBD(IBD)
      J2=KVBD(IBD+1)
      IF (DCORVG(1,J1).GT.DCORVG(1,J2)) THEN
       KVBD(IBD+1)=J1
       KVBD(IBD)=J2
       BEND=.FALSE.
      ENDIF
50    CONTINUE
      IF (.NOT.BEND) GOTO 49
C
      END
C
C
C
      SUBROUTINE SNVB(KNPR,KBCT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      DIMENSION KNPR(*),KBCT(*)
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /ERRCTL/ IER,ICHECK
      SAVE /TRIAD/,/ERRCTL/
C
      IF (ICHECK.GE.997) CALL OTRC('SNVB  ','01/02/89')
C
      KBCT(1)=1
      DO 10 IBCT=2,NBCT+1
10    KBCT(IBCT)=0
C
      DO 20 IVT=1,NVT
      INPR=KNPR(IVT)
      IF (INPR.GT.0) KBCT(INPR+1)=KBCT(INPR+1)+1
20    CONTINUE
C
      DO 30 IBCT=1,NBCT
30    KBCT(IBCT+1)=KBCT(IBCT)+KBCT(IBCT+1)
      NVBD=KBCT(NBCT+1)-1
C
      END
