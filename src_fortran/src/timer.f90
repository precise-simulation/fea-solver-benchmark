!***********************************************************************
!                                                                      *
! ZTIME                                                                *
!                                                                      *
! Purpose  Return cpu time  (machine dependent)                        *
!                                                                      *
! INPUT    TYPE                                                        *
! -----    ----                                                        *
! T        R*8    -1D0: Initialization                                 *
!                                                                      *
! OUTPUT   TYPE                                                        *
! ------   ----                                                        *
! T        R*8    Absolute cpu time                                    *
!                                                                      *
!***********************************************************************
SUBROUTINE ZTIME(T)

IMPLICIT NONE

REAL ( KIND = 8 ) T

!$    DOUBLE PRECISION omp_get_wtime
!$    T = omp_get_wtime()
!$    RETURN

!!$ REAL ( KIND = 4 ) VTA(2)
!!$ T = ETIME(VTA)

CALL CPU_TIME ( T )

END SUBROUTINE
