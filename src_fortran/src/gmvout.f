C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C     Output results in GMV ASCII format
C
      SUBROUTINE GMVOUT(MUNIT,CFILE,LU)
C=======================================================================
C     Declarations
C=======================================================================
      IMPLICIT DOUBLE PRECISION(A,C-H,O-U,W-Z),LOGICAL(B)
      PARAMETER (NNARR=999,NNLEV=19)
      INCLUDE 'nnwork.inc'
      CHARACTER CFILE*(*)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
C-----------------------------------------------------------------------
C     C O M M O N S
C-----------------------------------------------------------------------
C
C *** Standard COMMON blocks
      COMMON          NWORK,IWORK,IWMAX,L(NNARR),DWORK(NNWORK)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /TRIAA/  LCORVG,LCORMG,LVERT,LMID,LADJ,LVEL,LMEL,LNPR,LMM,
     *                LVBD,LEBD,LBCT,LVBDP,LMBDP
C-----------------------------------------------------------------------
C *** Standard dimensioning for workspace concept
      DIMENSION VWORK(1),KWORK(1)
      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C
      SUB = 'GMVOUT'
C
      OPEN (UNIT=MUNIT,FILE=CFILE)
C
      WRITE(MUNIT,'(A)')'gmvinput ascii'
C
C
C *** Write node coordinates
      WRITE(MUNIT,*)'nodes ',NVT
C
      DO 100 IVT=1,NVT
100   WRITE(MUNIT,1000) REAL(DWORK(L(LCORVG)+2*IVT-2))
      DO 101 IVT=1,NVT
101   WRITE(MUNIT,1000) REAL(DWORK(L(LCORVG)+2*IVT-1))
      DO 102 IVT=1,NVT
102   WRITE(MUNIT,1000) 0.E0
C
C
C *** Write cell information
      WRITE(MUNIT,*)'cells ',NEL
      DO 110 IEL=1,NEL
      WRITE(MUNIT,*)'quad 4'
110   WRITE(MUNIT,1100) KWORK(L(LVERT)+4*IEL-4),KWORK(L(LVERT)+4*IEL-3),
     *                  KWORK(L(LVERT)+4*IEL-2),KWORK(L(LVERT)+4*IEL-1)
C
C
C *** Write variable information
      WRITE(MUNIT,*) 'variable'
C
      WRITE(MUNIT,*) 'solution 1'
      DO 200 IVT=1,NVT
200     WRITE(MUNIT,1000) REAL(DWORK(L(LU)+IVT-1))
C
      WRITE(MUNIT,*)  'endvars'
      WRITE(MUNIT,*)  'probtime ',0.E0
      WRITE(MUNIT,*)  'endgmv'
C
      REWIND (MUNIT)
      CLOSE  (MUNIT)
C
1000  FORMAT(E15.8)
1100  FORMAT(8I8)
C
      GOTO 99999
C=======================================================================
C     Error case
C=======================================================================
99998 WRITE(*,*) 'IER', IER
      WRITE(*,*) 'IN SUBROUTINE ',SUB
C
99999 END
