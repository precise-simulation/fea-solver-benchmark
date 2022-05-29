      SUBROUTINE ASMMAT(A,KCOLA,KLDA,NA,NEQ,KVERT,KMID,DCORVG,
     *                  ELE,COEFF,BCON,ICUB,DENTOUT)
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNBAS=21,NNDER=6,NNCUBP=36,NNVE=4)
      INTEGER KCOLA(*),KLDA(*),NA,NEQ,KVERT(NNVE,*),KMID(NNVE,*),ICUB
      DOUBLE PRECISION A(*),DCORVG(2,*),COEFF,DENTOUT(4,4)
      LOGICAL BCON
C
      INTEGER          KENTRY(NNBAS,NNBAS)
      DOUBLE PRECISION COECON(3,3),DENTRY(NNBAS,NNBAS)
C
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /ELEM/   DX(NNVE),DY(NNVE),DJAC(2,2),DETJ,
     1                DBAS(NNBAS,NNDER),BDER(NNDER),KVE(NNVE),IEL
      COMMON /TRIAD/  NEL,NVT,NMT,NVE,NVEL,NBCT,NVBD
      COMMON /CUB/    DXI(NNCUBP,3),DOMEGA(NNCUBP),NCUBP,ICUBP
      COMMON /COAUX1/ KDFG(NNBAS),KDFL(NNBAS),IDFL
      SAVE /OUTPUT/,/ERRCTL/,/CHAR/,/ELEM/,/TRIAD/,/CUB/,/COAUX1/
C
      SUB='ASMMAT'
      IF (ICHECK.GE.997) CALL OTRC('ASMMAT','14/01/17')
C
C     Preparation - evaluation of parameters
      IER = 0
C     Which derivatives of basis functions are needed?
      DO 1 I=1,NNDER
      BDER(I) = .TRUE.
      BDER(2) = .TRUE.
1     BDER(3) = .TRUE.
C
C     For non sorting of KDFG and KDFL
      KDFL(1) = 1
      KDFL(2) = 2
      KDFL(3) = 3
      KDFL(4) = 4
C
C     Dummy call of ELE sets number of element
      IELTYP = -1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL = NDFL(IELTYP)
      CALL CB2Q(ICUB)
      IF (IER.NE.0) GOTO 99999
C
C     Dummy call of COEFF for nonlinear problems
C     COEFF must set BDER(IDER)=.TRUE. if derivative IDER is needed
      AUX = COEFF(0D0,0D0,-1,-1,0,BFIRST)
      IF (BCON) THEN
        IA = 2
        IB = 2
        COECON(IA,IB)=COEFF(0D0,0D0,IA,IB,1,.TRUE.)
        IA = 3
        IB = 3
        COECON(IA,IB)=COEFF(0D0,0D0,IA,IB,1,.TRUE.)
      ENDIF
C
C     Dummy call - ELE may save arithmetic operations
      ICUBP = ICUB
      CALL ELE(0D0,0D0,-2)
C
C *** Loop over all cells
      DO 100 IEL=1,NEL
C
C       Get local to global dof mapping
        CALL NDFGL(IEL,1,IELTYP,KVERT,KMID,KDFG,KDFL)   ! Don't sort for global vertex numbers
        IF (IER.LT.0) GOTO 99999
C
C       Determine entry positions in matrix
        DO 110 IROWL=1,IDFL
          IROW = KDFG(IROWL)
          IA   = KLDA(IROW)                  ! Pointer to diagonal entry in A(IROW,IROW), and start of row IROW
          KENTRY(IROWL,IROWL) = IA
          DENTRY(IROWL,IROWL) = 0D0
          DO 111 JCOLL=1,IDFL                ! For unsymmetric version
            IF ( JCOLL.EQ.IROWL ) GOTO 111   ! Diagonal already processed
C
            JCOL = KDFG(JCOLL)
            DO 112 IAJ=IA+1,NA   ! Loop over row entries
              IF ( KCOLA(IAJ).EQ.JCOL ) GOTO 113
112         CONTINUE
113         IA = IA+1   ! Not necessary if KDFL is sorted???
            KENTRY(IROWL,JCOLL) = IAJ
            DENTRY(IROWL,JCOLL) = 0D0
111       CONTINUE
110     CONTINUE
C
C       Extract vertex coordinates
        DO 120 IVE=1,NVE
          JP       = KVERT(IVE,IEL)
          KVE(IVE) = JP
          DX(IVE)  = DCORVG(1,JP)
          DY(IVE)  = DCORVG(2,JP)
120     CONTINUE
C
        DJ1 = 0.5D0*(-DX(1)-DX(2)+DX(3)+DX(4))
        DJ2 = 0.5D0*( DX(1)-DX(2)+DX(3)-DX(4))
        DJ3 = 0.5D0*(-DY(1)+DY(2)-DY(3)+DY(4))
        DJ4 = 0.5D0*(-DY(1)+DY(2)+DY(3)-DY(4))
C
C ***   Loop over all cubature points
        DO 200 ICUBP=1,NCUBP
C
          XI1 = DXI(ICUBP,1)
          XI2 = DXI(ICUBP,2)
C
C         Jacobian of the bilinear mapping onto the reference element
          DJAC(1,1) = 0.5D0*(DX(2)-DX(1)+DJ2)+0.5D0*DJ2*XI2
          DJAC(1,2) = 0.5D0*DJ1+0.5D0*DJ2*XI1
          DJAC(2,1) = 0.5D0*DJ4-0.5D0*DJ3*XI2
          DJAC(2,2) = 0.5D0*(DY(3)-DY(1)-DJ4)-0.5D0*DJ3*XI1
          DETJ      = DJAC(1,1)*DJAC(2,2)-DJAC(1,2)*DJAC(2,1)
          OM        = DOMEGA(ICUBP)*DETJ
C
C         ELE needs the information ICUBP because of preceeding
C         dummy call using IPAR = -2
          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999
          IF (BCON) THEN
            XX  = 0.5D0*(DX(1)+DX(2)+DJ1)+0.5D0*(DX(2)-DX(1)+DJ2)*XI1
     *           +0.5D0*DJ1*XI2+0.5D0*DJ2*XI1*XI2
            YY  = 0.5D0*(DY(1)+DY(3)+DJ3)+0.5D0*DJ4*XI1+0.5D0*
     *            (DY(3)-DY(1)-DJ4)*XI2-0.5D0*DJ3*XI1*XI2
            AUX = COEFF(XX,YY,IA,IB,IBLOC,BFIRST)
          ELSE
            AUX = COECON(IA,IB)
          ENDIF
C
C ***     Summing up over all pairs of multiindices
          DO 231 JDOFE=1,IDFL
            JDOFEH = KDFL(JDOFE)
            HBASJ2 = DBAS(JDOFEH,2)
            HBASJ3 = DBAS(JDOFEH,3)
C
            DO 241 IDOFE=1,IDFL
              IF ( IDOFE.EQ.JDOFE ) THEN
                AH     = AUX*(HBASJ2**2+HBASJ3**2)
!                AH     = AUX*DBAS(JDOFEH,1)**2
              ELSE
                IDOFEH = KDFL(IDOFE)
                HBASI2 = DBAS(IDOFEH,2)
                HBASI3 = DBAS(IDOFEH,3)
                AH     = AUX*(HBASJ2*HBASI2+HBASJ3*HBASI3)
!                AH     = AUX*DBAS(JDOFEH,1)*DBAS(IDOFEH,1)
              ENDIF
              DENTRY(IDOFE,JDOFE) = DENTRY(IDOFE,JDOFE)+OM*AH
241         CONTINUE
231       CONTINUE
C
200     CONTINUE   ! End loop over cubature points
C
        DO 400 JDOFE=1,IDFL
          DO 401 IDOFE=1,IDFL
            IA    = KENTRY(JDOFE,IDOFE)
            A(IA) = A(IA)+DENTRY(JDOFE,IDOFE)
C
            DENTOUT(KDFL(JDOFE),KDFL(IDOFE)) = DENTRY(JDOFE,IDOFE)
401       CONTINUE
400     CONTINUE
C
100   CONTINUE   ! End loop over cells
C
99999 END
